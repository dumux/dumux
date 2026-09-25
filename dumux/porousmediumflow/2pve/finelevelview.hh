// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Provides access to all fine-level concepts and containers
 */

#ifndef DUMUX_TWOPVE_FINE_LEVEL_VIEW_HH
#define DUMUX_TWOPVE_FINE_LEVEL_VIEW_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/porousmediumflow/2pve/quantityreconstruction.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>
#include <dumux/porousmediumflow/2pve/elementstatefine.hh>
#include <dumux/porousmediumflow/2pve/fieldstoragefine.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief The gas plume distance of a coarse-level column at the last time step and its minimum over all previous time steps
 *
 * The minimum gas plume distance bounds the region with residually trapped gas.
 */
template<class Scalar>
struct TwoPVEColumnHistory
{
    Scalar gasPlumeDistance{};
    Scalar minimumGasPlumeDistance{};
};

/*!
 * \ingroup TwoPVEModel
 * \brief Connects the coarse level of the two-phase VE model to its fine level
 *
 * Provides the mapping between coarse-level columns and fine-level elements, the history of each column,
 * and the reconstruction of the fine-level solution from the coarse-level solution.
 *
 * \tparam GridGeometry the grid geometry of both levels
 * \tparam Scalar the scalar type
 * \tparam FluidSystem the immiscible two-phase fluid system
 * \tparam Indices the primary variable indices of the model
 * \tparam SolutionVector the type of the coarse-level and the fine-level solution vectors
 * \tparam SpatialParamsFine the fine-level spatial parameters, providing `permeabilityAtElement(fineElement)`,
 *         `porosityAtElement(fineElement)` and `gridGeometry()`
 */
template<class GridGeometry, class Scalar, class FluidSystem, class Indices, class SolutionVector, class SpatialParamsFine>
class TwoPVEFineLevelView
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using QuantityReconstructor = TwoPVEQuantityReconstruction<Scalar, FluidSystem>;
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int saturationIdx = Indices::saturationIdx;
    static constexpr int wettingPhaseIdx = FluidSystem::phase0Idx;
    static constexpr int nonwettingPhaseIdx = FluidSystem::phase1Idx;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == 2 || dim == 3, "TwoPVEFineLevelView only supports two- and three-dimensional grids");
    using ColumnState = TwoPVEColumnState<Scalar>;
    using ColumnMapping = TwoPVEColumnMapping<GridGeometry, Scalar>;

    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:
    using FineLevelElementState = TwoPVEFineLevelElementState<GridGeometry, Scalar, FluidSystem>;
    using FineLevelFields = TwoPVEFineLevelFieldStorage<GridGeometry, Scalar, FluidSystem>;

    /*!
     * \brief Builds the mapping between the levels and checks that the grids fulfill the requirements of the model
     *
     * \param gridGeometry       fine-level grid geometry
     * \param gridGeometryCoarse coarse-level grid geometry
     * \param spatialParamsFine  fine-level spatial parameters
     */
    TwoPVEFineLevelView(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<const GridGeometry> gridGeometryCoarse,
                        std::shared_ptr<const SpatialParamsFine> spatialParamsFine)
    : gridGeometryFine_(gridGeometry),
      gridGeometryCoarse_(gridGeometryCoarse),
      spatialParamsFine_(spatialParamsFine),
      quantityReconstructor_{},
      columnMapping_(gridGeometryCoarse, gridGeometry),
      solution_(gridGeometryFine_->numDofs()),
      fineLevelFields_(gridGeometryFine_->numDofs()),
      fineCellHeight_(computeFineCellHeight_()),
      columnHistory_(gridGeometryCoarse_->numDofs())
    {
        const Scalar domainHeight = gridGeometryFine_->bBoxMax()[dim-1] - gridGeometryFine_->bBoxMin()[dim-1];
        checkColumns_(domainHeight);
        for (auto& history : columnHistory_)
        {
            // initially there is no gas in the columns
            history.gasPlumeDistance = domainHeight;
            history.minimumGasPlumeDistance = domainHeight;
        }
    }

    /*!
     * \brief Returns the height of a fine-level cell (for uniform grid)
     */
    Scalar fineCellHeight() const
    {
        return fineCellHeight_;
    }

    /*!
     * \brief Getter function for fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParams() const
    {
        return *spatialParamsFine_;
    }

    /*!
     * \brief Getter function for pointer to fine-level spatial parameters
     */
    std::shared_ptr<const SpatialParamsFine> spatialParamsPtr() const
    {
        return spatialParamsFine_;
    }

    /*!
     * \brief Getter function for mapping between coarse and fine level
     */
    const ColumnMapping& columnMap() const
    {
        return columnMapping_;
    }

    /*!
     * \brief Getter function for fine-level grid geometry
     */
    const GridGeometry& gridGeometry() const
    {
        return *gridGeometryFine_;
    }

    /*!
     * \brief Getter function for coarse-level grid geometry
     */
    const GridGeometry& gridGeometryCoarse() const
    {
        return *gridGeometryCoarse_;
    }

    /*!
     * \brief Getter function for fine-level solution vector
     */
    const SolutionVector& solution() const
    {
        return solution_;
    }

    /*!
     * \brief Getter function for fine-level solution vector
     */
    SolutionVector& solution()
    {
        return solution_;
    }

    /*!
     * \brief Getter function for pointer to quantity reconstructor member
     */
    const QuantityReconstructor& quantityReconstructor() const
    {
        return quantityReconstructor_;
    }

    /*!
     * \brief Getter function for fine-level solution fields
     */
    const auto& fields() const
    {
        return fineLevelFields_;
    }

    /*!
     * \brief Computes and returns all quantities that define the state of a coarse-level column
     *
     * \param coarseElement       coarse-level element
     * \param coarsePriVars       coarse-level primary variables
     * \param coarseSpatialParams coarse-level spatial parameters
     */
    template<class PrimaryVariables, class SpatialParamsCoarse>
    ColumnState makeColumnState(const Element& coarseElement,
                                const PrimaryVariables& coarsePriVars,
                                const SpatialParamsCoarse& coarseSpatialParams) const
    {
        ColumnState state;

        const auto coarsePosition = coarseElement.geometry().center();
        const auto columnIdx = gridGeometryCoarse_->elementMapper().index(coarseElement);

        state.pwCoarse = coarsePriVars[pressureIdx];
        state.swCoarse = 1.0 - coarsePriVars[saturationIdx];
        state.temperature = coarseSpatialParams.temperatureAtPos(coarsePosition);
        const auto& history = columnHistory_[columnIdx];

        state.domainHeight = gridGeometryFine_->bBoxMax()[dim - 1] - gridGeometryFine_->bBoxMin()[dim - 1];

        const auto fluidMatrixInteraction = coarseSpatialParams.fluidMatrixInteractionAtPos(coarsePosition);
        using BasicParams = std::decay_t<decltype(fluidMatrixInteraction.pcSwCurve().basicParams())>;
        static_assert(std::is_same_v<BasicParams, FluidMatrix::BrooksCorey::Params<Scalar>>, "The VE reconstruction requires a Brooks-Corey material law");
        const auto& brooksCoreyParams = fluidMatrixInteraction.pcSwCurve().basicParams();
        const auto& absoluteSaturationParams = fluidMatrixInteraction.pcSwCurve().effToAbsParams();
        state.swr = absoluteSaturationParams.swr();
        state.snr = absoluteSaturationParams.snr();
        state.brooksCoreyLambda = brooksCoreyParams.lambda();
        state.entryPressure = brooksCoreyParams.pcEntry();
        const auto& gravity = coarseSpatialParams.gravity(coarsePosition);
        for (int dirIdx = 0; dirIdx < dimWorld - 1; ++dirIdx)
            if (gravity[dirIdx] != 0.0)
                DUNE_THROW(Dune::InvalidStateException, "The VE model requires gravity to be aligned with the vertical coordinate axis");
        state.gravityNorm = gravity.two_norm();

        // the fluid properties of both phases are evaluated at the coarse-level wetting-phase pressure
        ImmiscibleFluidState<Scalar, FluidSystem> fluidState;
        fluidState.setTemperature(state.temperature);
        fluidState.setPressure(wettingPhaseIdx, state.pwCoarse);
        fluidState.setPressure(nonwettingPhaseIdx, state.pwCoarse);
        state.densityW = FluidSystem::density(fluidState, wettingPhaseIdx);
        state.viscosityW = FluidSystem::viscosity(fluidState, wettingPhaseIdx);
        state.densityNw = FluidSystem::density(fluidState, nonwettingPhaseIdx);
        state.viscosityNw = FluidSystem::viscosity(fluidState, nonwettingPhaseIdx);

        state.gasPlumeDistance = quantityReconstructor_.computeGasPlumeDist(
                 PhaseDensities{state.densityW, state.densityNw},
                 ResidualSaturations{state.swr, state.snr},
                 state.gravityNorm,
                 state.domainHeight,
                 state.swCoarse,
                 history.minimumGasPlumeDistance,
                 BrooksCoreyParameters{state.brooksCoreyLambda, state.entryPressure});
        using std::min;
        state.minimumGasPlumeDistance = min(history.minimumGasPlumeDistance, state.gasPlumeDistance);

        return state;
    }

    /*!
     * \brief Updates fine-level solution and solution fields given the coarse-level solution
     *
     * Also updates the history of each column, on which the coarse-level volume variables depend.
     * Cached coarse-level volume variables have to be updated afterwards.
     *
     * \param coarseProblem  coarse-level problem
     * \param coarseSolution coarse-level solution vector
     */
    template<class ProblemCoarse>
    void updateSol(const ProblemCoarse& coarseProblem,
                   const SolutionVector& coarseSolution)
    {
        const auto& coarseSpatialParams = coarseProblem.spatialParams();

        for (const auto& coarseElement : elements(gridGeometryCoarse_->gridView()))
        {
            const auto coarseIdx = gridGeometryCoarse_->elementMapper().index(coarseElement);
            const auto& coarsePriVars = coarseSolution[coarseIdx];
            const auto columnState = makeColumnState(coarseElement, coarsePriVars, coarseSpatialParams);

            auto& history = columnHistory_[coarseIdx];
            history.gasPlumeDistance = columnState.gasPlumeDistance;
            history.minimumGasPlumeDistance = columnState.minimumGasPlumeDistance;

            const auto& column = columnMapping_.column(coarseIdx);
            for (const auto& fineElement : column)
            {
                const auto fineIdx = gridGeometryFine_->elementMapper().index(fineElement);

                FineLevelElementState fineElementState;
                fineElementState.update(fineElement, columnState, *spatialParamsFine_, quantityReconstructor_, fineCellHeight_);
                solution_[fineIdx][pressureIdx] = fineElementState.pressure(wettingPhaseIdx);
                solution_[fineIdx][saturationIdx] = fineElementState.saturation(nonwettingPhaseIdx);
                fineLevelFields_.set(fineIdx, fineElementState);
            }
        }
    }

private:

    /*!
     * \brief Helper function for computing the extent of a fine-level cell along a major axis
     *
     * \param geometry  geometry of a fine-level element
     * \param direction index that represents a major axis
     */
    template<class ElementGeometry>
    Scalar elementExtent_(const ElementGeometry& geometry, const int direction) const
    {
        using std::max; using std::min;
        Scalar lower = geometry.corner(0)[direction];
        Scalar upper = lower;

        for (int cornerIdx = 1; cornerIdx < geometry.corners(); ++cornerIdx)
        {
            const Scalar coordinate = geometry.corner(cornerIdx)[direction];
            lower = min(lower, coordinate);
            upper = max(upper, coordinate);
        }

        return upper - lower;
    }

    /*!
     * \brief Checks that the fine-level elements have a uniform height and that each column spans the domain height
     *
     * \param domainHeight height of the domain
     */
    void checkColumns_(Scalar domainHeight) const
    {
        using std::abs;
        const Scalar tolerance = 1e-10*domainHeight;
        for (const auto& element : elements(gridGeometryFine_->gridView()))
            if (abs(elementExtent_(element.geometry(), dim - 1) - fineCellHeight_) > tolerance)
                DUNE_THROW(Dune::InvalidStateException, "The fine grid of the VE model has to be uniform in the vertical direction");

        for (std::size_t columnIdx = 0; columnIdx < columnMapping_.numberOfColumns(); ++columnIdx)
            if (abs(columnMapping_.column(columnIdx).size()*fineCellHeight_ - domainHeight) > tolerance)
                DUNE_THROW(Dune::InvalidStateException, "Column " << columnIdx << " does not span the domain height, the coarse grid of the VE model has to consist of a single layer of elements");
    }

    /*!
     * \brief Helper function for computing the height of uniform, fine-level elements
     */
    Scalar computeFineCellHeight_() const
    {
        const auto& gridView = gridGeometryFine_->gridView();
        const auto element = *gridView.template begin<0>();
        return elementExtent_(element.geometry(), dim - 1);
    }

    std::shared_ptr<const GridGeometry> gridGeometryFine_;
    std::shared_ptr<const GridGeometry> gridGeometryCoarse_;
    std::shared_ptr<const SpatialParamsFine> spatialParamsFine_;
    QuantityReconstructor quantityReconstructor_;
    ColumnMapping columnMapping_;
    SolutionVector solution_;
    FineLevelFields fineLevelFields_;
    Scalar fineCellHeight_;
    std::vector<TwoPVEColumnHistory<Scalar>> columnHistory_;
};

} // end namespace Dumux

#endif
