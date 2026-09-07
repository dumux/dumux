// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief Provides access to all fine-level concepts and containers
 */

#ifndef DUMUX_TWOPVE_FINE_LEVEL_VIEW_HH
#define DUMUX_TWOPVE_FINE_LEVEL_VIEW_HH

#include <algorithm>

#include <dumux/porousmediumflow/2pve/columnmapping.hh>
#include <dumux/porousmediumflow/2pve/elementstatefine.hh>
#include <dumux/porousmediumflow/2pve/fieldstoragefine.hh>

namespace Dumux {

// history is required for hysteresis
template<class Scalar>
struct TwoPVEColumnHistory
{
    Scalar gasPlumeDistance{};
    Scalar minimumGasPlumeDistance{};
};


template<class TypeTag, class FineProblemType>
class TwoPVEFineLevelView
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using QuantityReconstructor = TwoPVEQuantityReconst<TypeTag>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationGasIdx = Indices::saturationIdx,
        wettingPhaseIdx = FluidSystem::phase0Idx,
        nonwettingPhaseIdx = FluidSystem::phase1Idx
    };
    enum {
        dim = GridView::dimension,
    };
    using WettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::WettingPhase;
    using NonwettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::NonwettingPhase;
    using ColumnState = TwoPVEColumnState<Scalar>;
    using ColumnMapping = VEColumnMapping<GridGeometry, Scalar>;
    using SpatialParamsCoarse = GetPropType<TypeTag, Properties::SpatialParams>;
    using ProblemCoarse = GetPropType<TypeTag, Properties::Problem>;

    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:
    using FineLevelElementState = TwoPVEFineLevelElementState<TypeTag>;
    using FineLevelFields = TwoPVEFineLevelFieldStorage<TypeTag>;
    using FineProblem = FineProblemType;
    using SpatialParamsFine = typename FineProblem::SpatialParamsFine;

    TwoPVEFineLevelView(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<const GridGeometry> gridGeometryCoarse)
    : gridGeometryFine_(gridGeometry),
      gridGeometryCoarse_(gridGeometryCoarse),
      quantityReconstructor_{},
      columnMapping_(gridGeometryCoarse, gridGeometry),
      solution_(gridGeometryFine_->numDofs()),
      fineLevelFields_(gridGeometryFine_->numDofs()),
      fineCellHeight_(computeFineCellHeight_()),
      fineCellDepth_(computeFineCellDepth_()),
      columnHistory_(gridGeometryCoarse_->numDofs())
    {
        problemFine_ = std::make_unique<FineProblemType>(gridGeometryFine_, fineCellHeight_, fineCellDepth_);

        const Scalar domainHeight = gridGeometryFine_->bBoxMax()[dim-1] - gridGeometryFine_->bBoxMin()[dim-1];
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
    const Scalar& fineCellHeight() const
    {
        return fineCellHeight_;
    }

    /*!
     * \brief Returns the depth of a fine-level cell (for uniform grid)
     */
    const Scalar& fineCellDepth() const
    {
        return fineCellDepth_;
    }

    /*!
     * \brief Getter function for fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParams() const
    {
        return problemFine_->spatialParams();
    }

    /*!
     * \brief Getter function for pointer to fine-level spatial parameters
     */
    std::shared_ptr<const SpatialParamsFine> spatialParamsPtr() const
    {
        return problemFine_->spatialParamsPtr();
    }

    /*!
     * \brief Getter function for fine-level problem
     */
    const FineProblemType& problem() const
    {
        return *problemFine_;
    }

    /*!
     * \brief Getter function for fine-level problem
     */
    FineProblemType& problem()
    {
        return *problemFine_;
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
    ColumnState makeColumnState(const Element& coarseElement,
                                const PrimaryVariables& coarsePriVars,
                                const SpatialParamsCoarse& coarseSpatialParams) const
    {
        ColumnState state;

        const auto coarsePosition = coarseElement.geometry().center();
        const auto columnIdx = gridGeometryCoarse_->elementMapper().index(coarseElement);

        state.pwCoarse = coarsePriVars[pressureH2OIdx];
        state.swCoarse = 1.0 - coarsePriVars[saturationGasIdx];
        state.temperature = coarseSpatialParams.temperatureAtPos(coarsePosition);
        const auto& history = columnHistory_[columnIdx];

        state.domainHeight = gridGeometryFine_->bBoxMax()[dim - 1] - gridGeometryFine_->bBoxMin()[dim - 1];

        const auto fluidMatrixInteraction = coarseSpatialParams.fluidMatrixInteractionAtPos(coarsePosition);
        const auto& brooksCoreyParams = fluidMatrixInteraction.pcSwCurve().basicParams();
        const auto& absoluteSaturationParams = fluidMatrixInteraction.pcSwCurve().effToAbsParams();
        state.swr = absoluteSaturationParams.swr();
        state.snr = absoluteSaturationParams.snr();
        state.brooksCoreyLambda = brooksCoreyParams.lambda();
        state.entryPressure = brooksCoreyParams.pcEntry();
        state.gravityNorm = coarseSpatialParams.gravity(coarsePosition).two_norm();
        state.densityW = WettingPhase::density(state.temperature, state.pwCoarse);
        state.viscosityW = WettingPhase::viscosity(state.temperature, state.pwCoarse);
        state.densityNw = NonwettingPhase::density(state.temperature, state.pwCoarse);
        state.viscosityNw = NonwettingPhase::viscosity(state.temperature, state.pwCoarse);

        state.gasPlumeDistance = quantityReconstructor_.computeGasPlumeDist(
                 PhaseDensities{state.densityW, state.densityNw},
                 ResidualSaturations{state.swr, state.snr},
                 state.gravityNorm,
                 state.domainHeight,
                 state.swCoarse,
                 BrooksCoreyParameters{state.brooksCoreyLambda, state.entryPressure});
        state.minimumGasPlumeDistance = std::min(history.minimumGasPlumeDistance, state.gasPlumeDistance);

        return state;
    }

    /*!
     * \brief Updates fine-level solution and solution fields given the coarse-level solution
     *
     * \param coarseProblem  coarse-level problem
     * \param coarseSolution coarse-level solution vector
     */
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
            for(const auto& fineElement : column)
            {
                const auto fineIdx = gridGeometryFine_->elementMapper().index(fineElement);

                FineLevelElementState fineElementState;
                fineElementState.update(fineElement, columnState, problemFine_->spatialParams(), quantityReconstructor_, fineCellHeight_);
                solution_[fineIdx][pressureH2OIdx] = fineElementState.pressure(wettingPhaseIdx);
                solution_[fineIdx][saturationGasIdx] = fineElementState.saturation(nonwettingPhaseIdx);
                fineLevelFields_.set(fineIdx,fineElementState);
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
        Scalar lower = geometry.corner(0)[direction];
        Scalar upper = lower;

        for (int cornerIdx = 1; cornerIdx < geometry.corners(); ++cornerIdx)
        {
            const Scalar coordinate = geometry.corner(cornerIdx)[direction];
            lower = std::min(lower, coordinate);
            upper = std::max(upper, coordinate);
        }

        return upper - lower;
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

    /*!
     * \brief Helper function for computing the depth of uniform, fine-level elements
     */
    Scalar computeFineCellDepth_() const
    {
        if constexpr (dim == 2)
            return 1.0; // assume same depth for all cells, otherwise use extrusionFactor per element
        else
        {
            static_assert(dim == 3, "TwoPVEFineLevelView currently supports only 2D and 3D grids");
            const auto& gridView = gridGeometryFine_->gridView();
            const auto element = *gridView.template begin<0>();
            return elementExtent_(element.geometry(), dim - 2);
        }
    }

    std::shared_ptr<const GridGeometry> gridGeometryFine_;
    std::shared_ptr<const GridGeometry> gridGeometryCoarse_;
    QuantityReconstructor quantityReconstructor_;
    ColumnMapping columnMapping_;
    SolutionVector solution_;
    FineLevelFields fineLevelFields_;
    Scalar fineCellHeight_, fineCellDepth_;
    std::vector<TwoPVEColumnHistory<Scalar>> columnHistory_;
    std::unique_ptr<FineProblemType> problemFine_;
};

} // end namespace Dumux

#endif
