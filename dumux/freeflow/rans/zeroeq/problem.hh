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
#include <string>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/rans/common/ransmomentumproblem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds a zero-equation (algebraic mixing-length) RANS turbulence
 *        closure to a face-centered-staggered Navier-Stokes momentum problem.
 *
 * This model does not solve any extra transport equation and does not require a
 * separate coupled turbulence sub-model. The eddy viscosity at a given element is
 * computed algebraically from the wall-normal derivative of the flow-direction
 * velocity component and the distance to the nearest wall, following Prandtl's
 * mixing-length hypothesis with an optional Van Driest near-wall damping correction
 * - matching the exact formula used by releases/3.10:dumux/freeflow/rans/zeroeq/volumevariables.hh
 * (velGrad = |velocityGradients()[flowDirectionAxis][wallNormalAxis]|).
 *
 * Wall distance and the velocity gradients are taken from Dumux::RANSMomentumProblem, i.e.
 * with the same (releases/3.10) reconstruction as every other RANS closure. The axes are fixed,
 * runtime-configurable constants (RANS.WallNormalAxis, defaulting to the last coordinate
 * direction, with the flow direction assumed to be the complementary axis in 2D) rather than
 * the old code's per-element dynamic axis detection - a deliberate simplification for
 * straight-channel-type test geometries.
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
class ZeroEqProblem : public RANSMomentumProblem<TypeTag>
{
    using ParentType = RANSMomentumProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static_assert(GridGeometry::discMethod == DiscretizationMethods::fcstaggered,
                  "ZeroEqProblem is only implemented for the face-centered staggered discretization");

public:
    using ParentType::ParentType;

    //! \copydoc RANSMomentumProblem::updateStaticWallProperties
    void updateStaticWallProperties()
    {
        ParentType::updateStaticWallProperties();
        eddyViscosity_.assign(this->gridGeometry().elementMapper().size(), 0.0);
    }

    /*!
     * \brief Updates the wall properties of the shared base and, from them, the eddy viscosity
     *        at every element. Must be called with the current momentum solution/grid
     *        variables after (each iteration of) solving the momentum balance, before the
     *        next assembly that is meant to use the updated eddy viscosity.
     */
    template<class SolutionVector, class GridVariables>
    void updateDynamicWallProperties(const SolutionVector& sol, const GridVariables& gridVariables)
    {
        ParentType::updateDynamicWallProperties(sol, gridVariables);

        const std::string eddyViscosityModel = getParamFromGroup<std::string>(this->paramGroup(), "RANS.EddyViscosityModel", "vanDriest");
        const auto flowAxis = this->flowDirectionAxis();
        const auto wallAxis = this->wallNormalAxis();

        for (std::size_t eIdx = 0; eIdx < eddyViscosity_.size(); ++eIdx)
        {
            using std::abs;
            using std::max;
            using std::sqrt;

            // the friction velocity (for the Van Driest damping) is estimated from the velocity
            // gradient at the wall-adjacent element: uStar = sqrt(nu * |velGrad|_wall)
            const Scalar kinematicViscosity = this->kinematicViscosity(eIdx);
            const Scalar wallVelocityGradient = abs(this->velocityGradient(this->wallElementIndex(eIdx), flowAxis, wallAxis));
            const Scalar frictionVelocity = max(sqrt(kinematicViscosity*wallVelocityGradient), 1e-10);

            const Scalar y = this->wallDistance(eIdx);
            const Scalar yPlus = y*frictionVelocity/kinematicViscosity;
            const Scalar velocityGradient = abs(this->velocityGradient(eIdx, flowAxis, wallAxis));

            eddyViscosity_[eIdx] = this->storedDensity(eIdx)*kinematicEddyViscosity_(eddyViscosityModel, y, yPlus, velocityGradient);
        }
    }

    //! Adds the (constant per element, precomputed) eddy viscosity to the molecular one.
    //! On interior faces the inside and outside values are averaged, exactly like the molecular
    //! viscosity is in the coupling manager: otherwise the two staggered half-cells sharing a
    //! lateral face would see different viscosities, so the shear flux leaving one would not
    //! equal the one entering the other (a spurious body force where the eddy viscosity varies).
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    {
        Scalar eddyViscosity = dynamicEddyViscosity_(fvGeometry.scv(scvf.insideScvIdx()).elementIndex());
        if (!scvf.boundary())
            eddyViscosity = 0.5*(eddyViscosity + dynamicEddyViscosity_(fvGeometry.scv(scvf.outsideScvIdx()).elementIndex()));
        return ParentType::effectiveViscosity(element, fvGeometry, scvf) + eddyViscosity;
    }

    //! \copydoc effectiveViscosity(const Element&,const FVElementGeometry&,const SubControlVolumeFace&) const
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scv) + dynamicEddyViscosity_(scv.elementIndex()); }

    //! The (spatially constant per element, precomputed) eddy viscosity - public so a
    //! nonisothermal mass-side problem (Dumux::ZeroEqMassProblem) can forward it into the
    //! energy equation's eddy thermal conductivity; zero-eq is the only RANS model with no
    //! mass-domain representation of its own, since its eddy viscosity is computed and added
    //! purely on the momentum side (see effectiveViscosity() above).
    Scalar dynamicEddyViscosity(std::size_t eIdx) const
    { return eddyViscosity_[eIdx]; }

private:
    Scalar dynamicEddyViscosity_(std::size_t eIdx) const
    { return eddyViscosity_[eIdx]; }

    //! Prandtl mixing-length closure, with an optional Van Driest near-wall damping factor.
    Scalar kinematicEddyViscosity_(const std::string& model, Scalar y, Scalar yPlus, Scalar velocityGradient) const
    {
        using std::exp;
        using std::sqrt;

        if (model == "none")
            return 0.0;

        Scalar mixingLength = this->karmanConstant()*y;

        if (model == "prandtl")
        {
            // undamped mixing length, mixingLength = kappa*y, used as is
        }
        else if (model == "vanDriest")
            mixingLength *= (1.0 - exp(-yPlus/26.0)) / sqrt(1.0 - exp(-0.26*yPlus));
        else
            DUNE_THROW(Dune::NotImplemented, "RANS.EddyViscosityModel " << model << " (use \"none\", \"prandtl\" or \"vanDriest\")");

        return mixingLength*mixingLength*velocityGradient;
    }

    std::vector<Scalar> eddyViscosity_;
};

} // end namespace Dumux

#endif
