// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::RANSMomentumProblem
 */
#ifndef DUMUX_RANS_COMMON_MOMENTUM_PROBLEM_HH
#define DUMUX_RANS_COMMON_MOMENTUM_PROBLEM_HH

#include <vector>
#include <array>
#include <cmath>

#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/walldistance.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/fcstaggered/velocitygradients.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin providing the momentum-side wall-distance/velocity-gradient/vorticity
 *        bookkeeping shared by transport-equation-based RANS turbulence closures
 *        (one-equation Spalart-Allmaras, and, in the future, two-equation models).
 *
 * This is a from-scratch reimplementation of
 * releases/3.10:dumux/freeflow/rans/problem.hh's RANSProblemBase, restricted to what the
 * one- and two-equation models actually need (no flat-wall-bounded y+/u+ bookkeeping - see
 * whatisimplemented.md for the rationale). Like the old class, all
 * quantities here are computed explicitly from main.cc, once per time step
 * (updateDynamicWallProperties), from the *previous* converged solution - a deliberate,
 * lagged/Picard-in-time treatment matching releases/3.10's actual numerics (as opposed to
 * differentiating these terms through Newton), see whatisimplemented.md and
 * proposedimplementation.md for why this replaces the fully-implicit 3-domain design that
 * was attempted and parked before this.
 *
 * Unlike the old code (which reconstructed velocity gradients via finite differences over
 * a hand-rolled structured-neighbor search), gradients are obtained by directly reusing the
 * existing, unmodified FCStaggeredVelocityGradients helper (fullGradient=true) - the same
 * mechanism already used and proven working by Dumux::ZeroEqProblem
 * (dumux/freeflow/rans/zeroeq/problem.hh) - rather than re-deriving gradient reconstruction.
 *
 * Usage: derive your test problem from this class (or, more commonly, from a further mixin
 * like Dumux::OneEqMomentumProblem that itself extends this one), analogous to
 * Dumux::ZeroEqProblem. Two extension points must be implemented by the most derived problem
 * class: `isOnWall(scvf)`/`isOnWallAtPos(globalPos)`, exactly as for Dumux::ZeroEqProblem.
 */
template<class TypeTag>
class RANSMomentumProblem : public NavierStokesMomentumProblem<TypeTag>
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
                  "RANSMomentumProblem is only implemented for the face-centered staggered discretization");

    static constexpr int dim = GridView::dimension;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using WallDistanceHelper = WallDistance<GridGeometry>;
    using VelocityGradients = FCStaggeredVelocityGradients;

public:
    using ParentType::ParentType;

    //! The von Kármán constant.
    Scalar karmanConstant() const
    { return 0.41; }

    /*!
     * \brief Computes the (solution-independent) distance from every element to the
     *        nearest wall. Call exactly once, after the problem has been fully constructed.
     */
    void updateStaticWallProperties()
    {
        const WallDistanceHelper wallDistance(
            this->gridGeometry(), WallDistanceHelper::atElementCenters,
            [&](const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
            { return asImp_().isOnWall(scvf); }
        );

        wallDistance_ = wallDistance.wallDistance();
        wallElementIdx_.resize(wallDistance_.size());
        const auto& wallData = wallDistance.wallData();
        for (std::size_t eIdx = 0; eIdx < wallDistance_.size(); ++eIdx)
            wallElementIdx_[eIdx] = wallData[eIdx].eIdx;

        velocityGradientTensor_.assign(wallDistance_.size(), DimMatrix(0.0));
        vorticityTensorScalarProduct_.assign(wallDistance_.size(), 0.0);
        stressTensorScalarProduct_.assign(wallDistance_.size(), 0.0);
        storedDensity_.assign(wallDistance_.size(), 0.0);
        storedViscosity_.assign(wallDistance_.size(), 0.0);

        findNeighborIndices_();
    }

    /*!
     * \brief Updates the velocity gradient tensor, the vorticity tensor scalar product, and
     *        the stored (molecular) density/viscosity at every element, from the given
     *        (typically: the last converged) momentum solution. Called once per time step
     *        from main.cc, before the next assembly - not re-evaluated within Newton.
     */
    template<class SolutionVector, class GridVariables>
    void updateDynamicWallProperties(const SolutionVector& sol, const GridVariables& gridVariables)
    {
        const auto& gridGeometry = this->gridGeometry();
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);

            const auto gradV = elementVelocityGradientTensor_(fvGeometry, elemVolVars);
            velocityGradientTensor_[eIdx] = gradV;

            Scalar vorticitySquaredSum = 0.0;
            Scalar stressSquaredSum = 0.0;
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                {
                    const Scalar vorticity_ij = 0.5*(gradV[i][j] - gradV[j][i]);
                    vorticitySquaredSum += vorticity_ij*vorticity_ij;

                    const Scalar stress_ij = 0.5*(gradV[i][j] + gradV[j][i]);
                    stressSquaredSum += stress_ij*stress_ij;
                }
            vorticityTensorScalarProduct_[eIdx] = vorticitySquaredSum;
            stressTensorScalarProduct_[eIdx] = stressSquaredSum;

            const auto& scv = *scvs(fvGeometry).begin();
            storedDensity_[eIdx] = this->density(element, scv);
            // ParentType::effectiveViscosity() forwards to the coupling manager and returns
            // the molecular (mass-domain) viscosity only - exactly the same trick
            // Dumux::ZeroEqProblem already uses to obtain it without duplicating the lookup.
            storedViscosity_[eIdx] = ParentType::effectiveViscosity(element, fvGeometry, scv);
        }
    }

    Scalar wallDistance(std::size_t eIdx) const
    { return wallDistance_[eIdx]; }

    //! The element index of the wall-adjacent cell in eIdx's wall-normal column - needed by
    //! two-equation wall-function models (k-epsilon); unused by zero-/one-equation/k-omega.
    std::size_t wallElementIndex(std::size_t eIdx) const
    { return wallElementIdx_[eIdx]; }

    //! The neighboring element index along the given axis/side (0 = "lower", 1 = "upper"),
    //! ported from releases/3.10:dumux/freeflow/rans/problem.hh's RANSProblemBase - needed for
    //! k-epsilon's per-column near-wall-region/matching-point search on structured,
    //! axis-aligned (flat-wall-bounded) meshes.
    std::size_t neighborIndex(std::size_t eIdx, int axisIdx, int sideIdx) const
    { return neighborIdx_[eIdx][axisIdx][sideIdx]; }

    GlobalPosition cellCenter(std::size_t eIdx) const
    { return this->gridGeometry().element(eIdx).geometry().center(); }

    //! One component of the full local velocity gradient tensor already reconstructed in
    //! updateDynamicWallProperties() - needed by k-epsilon's zero-equation (Van Driest)
    //! blended eddy viscosity, mirroring ZeroEqProblem's own gradient usage.
    Scalar velocityGradient(std::size_t eIdx, int i, int j) const
    { return velocityGradientTensor_[eIdx][i][j]; }

    //! The wall-normal axis used by wall-function/mixing-length formulas (RANS.WallNormalAxis,
    //! defaulting to the last coordinate direction) - same fixed-axis convention already
    //! established by Dumux::ZeroEqProblem for straight, axis-aligned channels.
    int wallNormalAxis() const
    {
        static const int axis = getParamFromGroup<int>(this->paramGroup(), "RANS.WallNormalAxis", GridView::dimension - 1);
        return axis;
    }

    //! The flow-direction axis: the complementary axis to wallNormalAxis() in 2D.
    int flowDirectionAxis() const
    {
        static_assert(GridView::dimension == 2, "flowDirectionAxis() as the complementary axis is only well-defined in 2D");
        return 1 - wallNormalAxis();
    }

    Scalar vorticityTensorScalarProduct(std::size_t eIdx) const
    { return vorticityTensorScalarProduct_[eIdx]; }

    //! S_ij*S_ij with S_ij = 1/2*(dv_i/dx_j + dv_j/dx_i), needed by two-equation models'
    //! production term (P = 2*mu_t*S:S) - not used by the one-equation model.
    Scalar stressTensorScalarProduct(std::size_t eIdx) const
    { return stressTensorScalarProduct_[eIdx]; }

    Scalar kinematicViscosity(std::size_t eIdx) const
    { return storedViscosity_[eIdx]/storedDensity_[eIdx]; }

    Scalar storedDensity(std::size_t eIdx) const
    { return storedDensity_[eIdx]; }

    Scalar storedViscosity(std::size_t eIdx) const
    { return storedViscosity_[eIdx]; }

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
    //! Full local velocity gradient tensor (du_i/dx_j), reusing the existing, unmodified
    //! FCStaggeredVelocityGradients helper - the same one Dumux::ZeroEqProblem already uses.
    //! Falls back to zero if the element has no interior frontal scvf to seed the
    //! full-gradient reconstruction from (only possible in a mesh one cell wide).
    template<class ElementVolumeVariables>
    DimMatrix elementVelocityGradientTensor_(const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars) const
    {
        for (const auto& scv : scvs(fvGeometry))
            for (const auto& scvf : scvfs(fvGeometry, scv))
                if (scvf.isFrontal() && !scvf.boundary())
                    return VelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars, /*fullGradient=*/true);

        return DimMatrix(0.0);
    }

    //! Axis-aligned nearest-neighbor search via cell-center comparison, ported from
    //! releases/3.10:dumux/freeflow/rans/problem.hh's RANSProblemBase::findNeighborIndices_() -
    //! only meaningful/used for structured, axis-aligned (flat-wall-bounded) meshes, exactly as
    //! in the old code.
    void findNeighborIndices_()
    {
        using std::abs;
        const auto& gridGeometry = this->gridGeometry();
        neighborIdx_.resize(gridGeometry.elementMapper().size());

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            for (int axisIdx = 0; axisIdx < dim; ++axisIdx)
            {
                neighborIdx_[eIdx][axisIdx][0] = eIdx;
                neighborIdx_[eIdx][axisIdx][1] = eIdx;
            }

            for (const auto& intersection : intersections(gridGeometry.gridView(), element))
            {
                if (intersection.boundary())
                    continue;

                const auto neighborIdx = gridGeometry.elementMapper().index(intersection.outside());
                for (int axisIdx = 0; axisIdx < dim; ++axisIdx)
                {
                    if (abs(cellCenter(eIdx)[axisIdx] - cellCenter(neighborIdx)[axisIdx]) > 1e-8)
                    {
                        if (cellCenter(eIdx)[axisIdx] > cellCenter(neighborIdx)[axisIdx])
                            neighborIdx_[eIdx][axisIdx][0] = neighborIdx;
                        if (cellCenter(eIdx)[axisIdx] < cellCenter(neighborIdx)[axisIdx])
                            neighborIdx_[eIdx][axisIdx][1] = neighborIdx;
                    }
                }
            }
        }
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::vector<Scalar> wallDistance_;
    std::vector<std::size_t> wallElementIdx_;
    std::vector<std::array<std::array<std::size_t, 2>, dim>> neighborIdx_;
    std::vector<DimMatrix> velocityGradientTensor_;
    std::vector<Scalar> vorticityTensorScalarProduct_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> storedDensity_;
    std::vector<Scalar> storedViscosity_;
};

} // end namespace Dumux

#endif
