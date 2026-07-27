// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::ZeroEqProblem
 */
#ifndef DUMUX_RANS_ZEROEQ_PROBLEM_HH
#define DUMUX_RANS_ZEROEQ_PROBLEM_HH

#include <cmath>
#include <vector>
#include <algorithm>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/walldistance.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/fcstaggered/velocitygradients.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds a zero-equation (algebraic mixing-length) RANS turbulence
 *        closure to a face-centered-staggered Navier-Stokes momentum problem.
 *
 * This model does not solve any extra transport equation and does not require a
 * separate coupled turbulence sub-model. The eddy viscosity at a given element is
 * computed purely algebraically from the local mean-strain-rate magnitude and the
 * distance to the nearest wall, following Prandtl's mixing-length hypothesis with
 * an optional Van Driest near-wall damping correction. See turbulenceequations.md
 * (repository root) for the underlying physics and exact formulas, and
 * whatisimplemented.md for how this maps onto (and deviates from) the turbulence
 * models that existed in DuMux at releases/3.10 before the mass/momentum-split
 * discretization refactor.
 *
 * Usage: derive your test problem from this class instead of directly from
 * Dumux::NavierStokesMomentumProblem<TypeTag>, e.g.
 * \code{.cc}
 * struct Problem<TypeTag, TTag::MyTestMomentum>
 * { using type = MyTestProblem<TypeTag, Dumux::ZeroEqProblem<TypeTag>>; };
 * \endcode
 * Two extension points must be implemented by the most derived problem class:
 * - `isOnWallAtPos(globalPos)` (or `isOnWall(scvf)`), returning whether a given
 *   boundary position/scvf is a solid (no-slip) wall relevant for wall-distance and
 *   near-wall damping computations (as opposed to an inlet/outlet/symmetry boundary).
 * Additionally, the two update routines below must be called explicitly from `main.cc`:
 * `updateStaticWallProperties()` once after the problem has been fully constructed,
 * and `updateDynamicWallProperties(sol, gridVariables)` after every (converged, or
 * Picard-iteration) solve of the momentum balance, before the next assembly that
 * should use the updated eddy viscosity.
 */
template<class TypeTag>
class ZeroEqProblem : public NavierStokesMomentumProblem<TypeTag>
{
    using ParentType = NavierStokesMomentumProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static_assert(GridGeometry::discMethod == DiscretizationMethods::fcstaggered,
                  "ZeroEqProblem is only implemented for the face-centered staggered discretization");

    using WallDistanceHelper = WallDistance<GridGeometry>;
    using WallData = typename WallDistanceHelper::WallData;
    using VelocityGradients = FCStaggeredVelocityGradients;

public:
    using ParentType::ParentType;

    //! The von Kármán constant, used throughout the mixing-length/near-wall formulas.
    Scalar karmanConstant() const
    { return 0.41; }

    /*!
     * \brief Computes the distance from every element to the nearest wall, and records
     *        which element/scvf that nearest wall belongs to (needed to evaluate the
     *        near-wall friction velocity used by the Van Driest damping function).
     * \note Purely geometric (does not depend on the solution) - call this exactly once,
     *       after the problem has been fully constructed (it accesses the CRTP-derived
     *       implementation via isOnWall/isOnWallAtPos, which requires complete construction).
     */
    void updateStaticWallProperties()
    {
        const WallDistanceHelper wallDistance(
            this->gridGeometry(), WallDistanceHelper::atElementCenters,
            [&](const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
            { return asImp_().isOnWall(scvf); }
        );

        wallDistance_ = wallDistance.wallDistance();
        wallData_ = wallDistance.wallData();

        strainRateMagnitude_.assign(wallDistance_.size(), 0.0);
        eddyViscosity_.assign(wallDistance_.size(), 0.0);
    }

    /*!
     * \brief Updates the local mean-strain-rate magnitude and, from it, the eddy viscosity
     *        at every element. Must be called with the current momentum solution/grid
     *        variables after (each iteration of) solving the momentum balance, before the
     *        next assembly that is meant to use the updated eddy viscosity.
     */
    template<class SolutionVector, class GridVariables>
    void updateDynamicWallProperties(const SolutionVector& sol, const GridVariables& gridVariables)
    {
        const auto& gridGeometry = this->gridGeometry();
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());

        // first pass: local mean-strain-rate magnitude |S| = sqrt(2 S_ij S_ij) at every element
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            strainRateMagnitude_[eIdx] = elementStrainRateMagnitude_(fvGeometry, elemVolVars);
        }

        // second pass: eddy viscosity, using the wall-adjacent element's strain rate
        // to estimate the friction velocity needed for the Van Driest damping function
        const std::string eddyViscosityModel = getParamFromGroup<std::string>(this->paramGroup(), "RANS.EddyViscosityModel", "vanDriest");

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);

            const auto& scv = *scvs(fvGeometry).begin();
            const Scalar density = this->density(element, scv);
            const Scalar molecularViscosity = ParentType::effectiveViscosity(element, fvGeometry, scv);
            const Scalar kinematicViscosity = molecularViscosity/density;

            const Scalar y = wallDistance_[eIdx];
            const Scalar wallStrainRateMagnitude = strainRateMagnitude_[wallData_[eIdx].eIdx];
            const Scalar frictionVelocity = std::max(std::sqrt(kinematicViscosity*wallStrainRateMagnitude), 1e-10);
            const Scalar yPlus = y*frictionVelocity/kinematicViscosity;

            eddyViscosity_[eIdx] = density*kinematicEddyViscosity_(eddyViscosityModel, y, yPlus, strainRateMagnitude_[eIdx]);
        }
    }

    //! Adds the (spatially constant per element, precomputed) eddy viscosity to the molecular one.
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scvf) + dynamicEddyViscosity_(fvGeometry.elementIndex()); }

    //! \copydoc effectiveViscosity(const Element&,const FVElementGeometry&,const SubControlVolumeFace&) const
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scv) + dynamicEddyViscosity_(scv.elementIndex()); }

    /*!
     * \brief Returns whether the given boundary scvf lies on a solid (no-slip) wall.
     * \note Overload this in the most derived problem, or overload isOnWallAtPos directly.
     */
    bool isOnWall(const SubControlVolumeFace& scvf) const
    { return asImp_().isOnWallAtPos(scvf.center()); }

    //! \copydoc isOnWall(const SubControlVolumeFace&) const
    bool isOnWallAtPos(const GlobalPosition&) const
    { DUNE_THROW(Dune::NotImplemented, "isOnWallAtPos or isOnWall not implemented for this problem."); }

private:
    Scalar dynamicEddyViscosity_(std::size_t eIdx) const
    { return eddyViscosity_[eIdx]; }

    //! Prandtl mixing-length closure, with an optional Van Driest near-wall damping factor.
    Scalar kinematicEddyViscosity_(const std::string& model, Scalar y, Scalar yPlus, Scalar strainRateMagnitude) const
    {
        using std::exp;
        using std::sqrt;

        if (model == "none")
            return 0.0;

        Scalar mixingLength = karmanConstant()*y;

        if (model == "prandtl")
        {
            // undamped mixing length, mixingLength = kappa*y, used as is
        }
        else if (model == "vanDriest")
            mixingLength *= (1.0 - exp(-yPlus/26.0)) / sqrt(1.0 - exp(-0.26*yPlus));
        else
            DUNE_THROW(Dune::NotImplemented, "RANS.EddyViscosityModel " << model << " (use \"none\", \"prandtl\" or \"vanDriest\")");

        return mixingLength*mixingLength*strainRateMagnitude;
    }

    //! Computes |S| = sqrt(2 S_ij S_ij) for the element, using the full local velocity
    //! gradient tensor reconstructed from the surrounding staggered velocity dofs.
    //! Falls back to zero if the element has no interior frontal scvf to seed the
    //! full-gradient reconstruction from (only possible in a mesh one cell wide).
    template<class ElementVolumeVariables>
    Scalar elementStrainRateMagnitude_(const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars) const
    {
        for (const auto& scv : scvs(fvGeometry))
        {
            for (const auto& scvf : scvfs(fvGeometry, scv))
            {
                if (scvf.isFrontal() && !scvf.boundary())
                {
                    const auto gradV = VelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars, /*fullGradient=*/true);

                    Scalar strainRateSquared = 0.0;
                    for (int i = 0; i < GridView::dimension; ++i)
                        for (int j = 0; j < GridView::dimension; ++j)
                        {
                            const Scalar sIJ = 0.5*(gradV[i][j] + gradV[j][i]);
                            strainRateSquared += sIJ*sIJ;
                        }

                    using std::sqrt;
                    return sqrt(2.0*strainRateSquared);
                }
            }
        }

        return 0.0;
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::vector<Scalar> wallDistance_;
    std::vector<WallData> wallData_;
    std::vector<Scalar> strainRateMagnitude_;
    std::vector<Scalar> eddyViscosity_;
};

} // end namespace Dumux

#endif
