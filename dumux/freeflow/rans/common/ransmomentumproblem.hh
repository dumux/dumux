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

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin providing the momentum-side wall-distance/velocity-gradient/vorticity
 *        bookkeeping shared by transport-equation-based RANS turbulence closures
 *        (one-equation Spalart-Allmaras, and, in the future, two-equation models).
 *
 * This provides the subset of wall-distance/near-wall bookkeeping the one- and two-equation
 * models actually need (no flat-wall-bounded y+/u+ bookkeeping). All quantities here are
 * computed explicitly from main.cc, once per time step (updateDynamicWallProperties), from
 * the *previous* converged solution - a deliberate, lagged/Picard-in-time treatment rather
 * than differentiating these terms through Newton.
 *
 * Velocity gradients are reconstructed by a cell-center-averaged central difference over each
 * element's two neighbors along the differentiated axis (releases/3.10's own
 * calculateCCVelocities_()/calculateCCVelocityGradients_() convention, reproduced verbatim -
 * including its Dirichlet-velocity-boundary handling, which differences the boundary value
 * against the next interior element rather than the boundary-adjacent element's own value),
 * rather than the generic FCStaggeredVelocityGradients helper: the latter's shorter effective
 * distance at boundary-adjacent cells was found to bias every model's near-wall turbulence
 * level, and for the low-Re k-epsilon closure to destabilize the wall-adjacent cell outright.
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
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using WallDistanceHelper = WallDistance<GridGeometry>;

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
        const auto numElements = gridGeometry.elementMapper().size();

        // First pass: each element's own per-axis dof value and dof position - the latter only
        // to determine, below, which side of the cell that dof sits on (needed to pair it with
        // the correct neighbor when averaging). Everything here is read from the currently-bound
        // element alone, no neighbor access required yet.
        std::vector<DimVector> ownVelocity(numElements);
        std::vector<DimVector> ownDofPosition(numElements);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            for (const auto& scv : scvs(fvGeometry))
            {
                ownVelocity[eIdx][scv.dofAxis()] = elemVolVars[scv].velocity();
                ownDofPosition[eIdx][scv.dofAxis()] = scv.dofPosition()[scv.dofAxis()];
            }
        }

        // Second pass: the cell-center-averaged velocity at every element, reproducing
        // releases/3.10:dumux/freeflow/rans/problem.hh's calculateCCVelocities_() ("faces are
        // equidistant to cell center"). This element's own dof sits on one side of the cell; the
        // dof on the other side is the same-axis dof of the neighboring element on that side.
        std::vector<DimVector> ccVelocity(numElements, DimVector(0.0));
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            for (int i = 0; i < dim; ++i)
            {
                const bool ownDofOnUpperSide = ownDofPosition[eIdx][i] > cellCenter(eIdx)[i];
                const auto otherIdx = neighborIndex(eIdx, i, ownDofOnUpperSide ? 0 : 1);
                ccVelocity[eIdx][i] = 0.5*(ownVelocity[eIdx][i] + ownVelocity[otherIdx][i]);
            }
        }

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);

            // Central-difference gradient over the *two* neighboring elements' cell-center
            // velocities (skipping this element's own value entirely), reproducing
            // releases/3.10:dumux/freeflow/rans/problem.hh's calculateCCVelocityGradients_().
            DimMatrix gradV(0.0);
            for (int j = 0; j < dim; ++j)
            {
                const auto neighborIdx0 = neighborIndex(eIdx, j, 0);
                const auto neighborIdx1 = neighborIndex(eIdx, j, 1);
                const auto distance = cellCenter(neighborIdx1)[j] - cellCenter(neighborIdx0)[j];

                if (std::abs(distance) < 1e-8)
                    continue;

                for (int i = 0; i < dim; ++i)
                    gradV[i][j] = (ccVelocity[neighborIdx1][i] - ccVelocity[neighborIdx0][i])/distance;
            }

            // At every Dirichlet-velocity boundary face (walls, and - on these axis-aligned
            // channel meshes - the inlet too), releases/3.10 instead differenced the boundary
            // value directly against the *next* interior element - one element further away
            // than the central-difference stencil above would otherwise reach.
            for (const auto& scv : scvs(fvGeometry))
            {
                for (const auto& scvf : scvfs(fvGeometry, scv))
                {
                    if (!scvf.boundary())
                        continue;

                    const auto axisIdx = scvf.normalAxis();
                    const auto bcTypes = asImp_().boundaryTypes(element, scvf);
                    bool anyDirichletVelocity = false;
                    for (int i = 0; i < dim; ++i)
                        anyDirichletVelocity = anyDirichletVelocity || bcTypes.isDirichlet(Indices::velocity(i));
                    if (!anyDirichletVelocity)
                        continue;

                    auto neighborIdx = neighborIndex(eIdx, axisIdx, 0);
                    if (scvf.center()[axisIdx] < cellCenter(eIdx)[axisIdx])
                        neighborIdx = neighborIndex(eIdx, axisIdx, 1);
                    const auto distance = cellCenter(neighborIdx)[axisIdx] - scvf.center()[axisIdx];
                    const auto dirichletValues = asImp_().dirichlet(element, scvf);

                    for (int i = 0; i < dim; ++i)
                        if (bcTypes.isDirichlet(Indices::velocity(i)))
                            gradV[i][axisIdx] = (ccVelocity[neighborIdx][i] - dirichletValues[Indices::velocity(i)])/distance;
                }
            }

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
