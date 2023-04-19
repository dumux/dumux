// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CO2Tests
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */

#ifndef DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH
#define DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH

#include <dumux/io/grid/griddata.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dumux/porousmediumflow/co2/model.hh>

namespace Dumux {

/*!
 * \ingroup CO2Tests
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */
template<class GridGeometry, class Scalar>
class HeterogeneousSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry,
                                       Scalar,
                                       HeterogeneousSpatialParams<GridGeometry, Scalar>>
{
    using Grid = typename GridGeometry::Grid;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, HeterogeneousSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using PermeabilityType = Scalar;

    HeterogeneousSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<const GridData<Grid>> gridData)
    : ParentType(gridGeometry)
    , gridData_(gridData)
    , pcKrSwCurve_("SpatialParams")
    {
        // Set the permeability for the layers
        barrierTopK_ = 1e-17; //sqm
        barrierMiddleK_ = 1e-15; //sqm
        reservoirK_ = 1e-14; //sqm

        // Set the effective porosity of the layers
        barrierTopPorosity_ = 0.001;
        barrierMiddlePorosity_ = 0.05;
        reservoirPorosity_ = 0.2;

        depthBOR_ = getParam<Scalar>("Problem.DepthBOR");
    }

    /*!
     * \brief Reads layer information from the grid
     */
    void getParamsFromGrid()
    {
        const auto& gridView = this->gridGeometry().gridView();
        paramIdx_.resize(gridView.size(0));

        for (const auto& element : elements(gridView))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            paramIdx_[eIdx] = gridData_->parameters(element)[0];
        }
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     *
     * \return The intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        // Get the global index of the element
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return permeability(eIdx);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param eIdx The element index
     */
    PermeabilityType permeability(std::size_t eIdx) const
    {
        if (paramIdx_[eIdx] == barrierTop_)
            return barrierTopK_;
        else if (paramIdx_[eIdx] == barrierMiddle_)
            return barrierMiddleK_;
        else
            return reservoirK_;
    }

    /*!
     * \brief Returns the volume fraction of the inert component with index compIdx \f$[-]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The element solution
     * \param compIdx The solid component index
     * \return The solid volume fraction
     */
    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        // Get the global index of the element
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return inertVolumeFraction(eIdx);
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param eIdx The element index
     */
    Scalar inertVolumeFraction(std::size_t eIdx) const
    {
        if (paramIdx_[eIdx] == barrierTop_)
            return 1- barrierTopPorosity_;
        else if (paramIdx_[eIdx] == barrierMiddle_)
            return 1- barrierMiddlePorosity_;
        else
            return 1- reservoirPorosity_;

    }


    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {   return makeFluidMatrixInteraction(pcKrSwCurve_); }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::BrineIdx; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The global position
     *
     * This problem assumes a geothermal gradient with
     * a surface temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return 283.0 + (depthBOR_ - globalPos[dimWorld-1])*0.03; }

private:

    std::shared_ptr<const GridData<Grid>> gridData_;

    int barrierTop_ = 1;
    int barrierMiddle_ = 2;
    int reservoir_ = 3;

    Scalar barrierTopPorosity_;
    Scalar barrierMiddlePorosity_;
    Scalar reservoirPorosity_;

    Scalar barrierTopK_;
    Scalar barrierMiddleK_;
    Scalar reservoirK_;

    Scalar depthBOR_;

    std::vector<int> paramIdx_;
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
