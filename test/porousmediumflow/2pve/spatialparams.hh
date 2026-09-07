// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief The spatial params for the vertical equilibrium Darcy test.
 */

#ifndef DUMUX_TEST_TWOPVE_SPATIAL_PARAMS_HH
#define DUMUX_TEST_TWOPVE_SPATIAL_PARAMS_HH

#include <memory>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>

#include "spatialparams_fine.hh"

namespace Dumux {

template<class GridGeometry, class Scalar>
class TwoPTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TwoPTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ThisType = TwoPTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;
    using SpatialParamsFine = TwoPTestFineSpatialParams<GridGeometry, Scalar>;
    using ColumnMapping = VEColumnMapping<GridGeometry, Scalar>;

public:
    using PermeabilityType = Scalar;

    TwoPTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                          const ColumnMapping& columnMapping,
                          std::shared_ptr<const SpatialParamsFine> spatialParamsFine,
                          const Scalar fineCellHeight)
        : ParentType(gridGeometry),
          gridGeometry_(gridGeometry),
          columnMapping_(columnMapping),
          spatialParamsFine_(spatialParamsFine),
          myPcKrSwCurve_("SpatialParams"),
          deltaZ_(fineCellHeight)
    {}

    /*!
     * \brief Return a reference to the fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParamsFine() const
    { return *spatialParamsFine_; }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param element coarse-level element
     * \param scv     respective subcontrol volume
     * \param elemSol local solution for this element
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        return permeabilityCoarseAtElement(element);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor for coarse-level element
     *
     * \param element coarse-level element
     */
    decltype(auto) permeabilityCoarseAtElement(const Element& element) const
    {
        const unsigned int elementIdx = gridGeometry_->elementMapper().index(element);
        return permeabilityCoarsePreComputed_[elementIdx];
    }

    /*!
     * \brief Returns the porosity
     *
     * \param element coarse-level element
     * \param scv     respective subcontrol volume
     * \param elemSol local solution for this element
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosityCoarseAtElement(element);
    }

    /*!
     * \brief Returns the porosity for coarse-level element
     *
     * \param element coarse-level element
     */
    Scalar porosityCoarseAtElement(const Element& element) const
    {
        const unsigned int elementIdx = gridGeometry_->elementMapper().index(element);
        return porosityCoarsePreComputed_[elementIdx];
    }


    /*!
     * \brief Computes the coarse-level permeability and porosity
     */
    void calcSpatialParamsPermeabilityAndPorosityCoarse()
    {
        permeabilityCoarsePreComputed_ = preCalculatePermeabilityCoarse(gridGeometry_);
        porosityCoarsePreComputed_ = preCalculatePorosityCoarse(gridGeometry_);
    }

    /*!
     * \brief Returns a vector of permeabilities for the coarse-grid \f$[m^2]\f$ elements
     *
     * \param gridGeometry gridGeometry belonging to the coarse-level grid
     */
    std::vector<PermeabilityType> preCalculatePermeabilityCoarse(std::shared_ptr<const GridGeometry> gridGeometry)
    {
        std::vector<PermeabilityType> storeCoarsePermeabilities(0.0);
        int numberCoarseElements = gridGeometry->elementMapper().size();
        storeCoarsePermeabilities.resize(numberCoarseElements);

        for (const auto& element : Dune::elements(gridGeometry->gridView()))
        {
            const int coarseElementIdx = gridGeometry->elementMapper().index(element);
            storeCoarsePermeabilities[coarseElementIdx] = calculatePermeabilityCoarseAtElement(coarseElementIdx);
        }

        return storeCoarsePermeabilities;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor for the coarse-grid \f$[m^2]\f$ at the given index
     *
     * \param coarseElementIdx index of the coarse-level element
     */
    PermeabilityType calculatePermeabilityCoarseAtElement(const unsigned int& coarseElementIdx) const
    {
        const auto& fineElementsInColumn = columnMapping_.column(coarseElementIdx);
        Scalar integratedValue = 0.0;

        for(const auto& fineElement : fineElementsInColumn)
        {
            integratedValue += spatialParamsFine_->permeabilityAtElement(fineElement)*deltaZ_;
        }

        return integratedValue;
    }

    /*!
     * \brief Returns a vector of porosities for the coarse-grid \f$[-]\f$ elements
     *
     * \param gridGeometry gridGeometry belonging to the coarse-level grid
     */
    std::vector<Scalar> preCalculatePorosityCoarse(std::shared_ptr<const GridGeometry> gridGeometry)
    {
        std::vector<Scalar> storeCoarsePorosities(0.0);
        int numberCoarseElements = gridGeometry->elementMapper().size();
        storeCoarsePorosities.resize(numberCoarseElements);

        for (const auto& element : Dune::elements(gridGeometry->gridView()))
        {
            const int coarseElementIdx = gridGeometry->elementMapper().index(element);
            storeCoarsePorosities[coarseElementIdx] = calculatePorosityCoarseAtElement(coarseElementIdx);
        }

        return storeCoarsePorosities;
    }


    /*!
     * \brief Returns the porosity for the coarse-grid \f$[-]\f$ at the given index
     *
     * \param coarseElementIdx index of the coarse-level element
     */
    Scalar calculatePorosityCoarseAtElement(const unsigned int& coarseElementIdx) const
    {
        const auto& fineElementsInColumn = columnMapping_.column(coarseElementIdx);
        Scalar integratedValue = 0.0;

        for(const auto& fineElement : fineElementsInColumn)
        {
            integratedValue += spatialParamsFine_->porosityAtElement(fineElement)*deltaZ_;
        }

        return integratedValue;
    }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/saturation material law
     *
     * \param globalPos the coordinates
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(myPcKrSwCurve_);
    }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/saturation material law
     */
    const auto fluidMatrixInteraction() const
    {
        return makeFluidMatrixInteraction(myPcKrSwCurve_);
    }

    /*
     * \brief Returns the object which contains the capillary pressure - saturation and relative permeability - saturation curves
     */
    const PcKrSwCurve getPcKrSwCurve() const
    {
        return myPcKrSwCurve_;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPosition) const
    {
        return 326.0; // 53°C
    }

    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        return extrusionFactorAtElement(element);
    }

    Scalar extrusionFactorAtElement(const Element& element) const
    {
        return 1.0;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }


private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    ColumnMapping columnMapping_;
    std::shared_ptr<const SpatialParamsFine> spatialParamsFine_;
    const PcKrSwCurve myPcKrSwCurve_;
    Scalar deltaZ_;
    std::vector<PermeabilityType> permeabilityCoarsePreComputed_;
    std::vector<Scalar> porosityCoarsePreComputed_;
};

} // end namespace Dumux

#endif
