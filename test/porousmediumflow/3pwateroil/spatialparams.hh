// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePWaterOilTests
 * \brief Definition of the spatial parameters for the steam-assisted gravity
 *        drainage (SAGD) problem.
 */

#ifndef DUMUX_SAGD_SPATIAL_PARAMS_HH
#define DUMUX_SAGD_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>

#include <dumux/material/fluidmatrixinteractions/3p/parkervangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilTests
 * \brief Definition of the spatial parameters for the steam-assisted gravity
 *        drainage (SAGD) problem
 */
template<class GridGeometry, class Scalar>
class SagdSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           SagdSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    enum { dimWorld=GridView::dimensionworld };

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                         SagdSpatialParams<GridGeometry, Scalar>>;

    using ThreePhasePcKrSw = FluidMatrix::ParkerVanGenuchten3PDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    SagdSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    , eps_(1e-6)
    {
        layerBottom_ = 35.0;

        // intrinsic permeabilities
        fineK_ = 1e-16;
        coarseK_ = 4e-14;

        // porosities
        finePorosity_ = 0.10;
        coarsePorosity_ = 0.1;
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {  return permeabilityAtPos(scv.dofPosition());}

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

private:

    bool isFineMaterial_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] > layerBottom_ - eps_; }

    Scalar layerBottom_;

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    const ThreePhasePcKrSw pcKrSwCurve_;

    Scalar eps_;
};

}

#endif
