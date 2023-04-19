// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem
 *        which uses the isothermal two-phase two component fully implicit model.
 */

#ifndef DUMUX_INFILTRATION_THREEPTHREEC_SPATIAL_PARAMETERS_HH
#define DUMUX_INFILTRATION_THREEPTHREEC_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class GridGeometry, class Scalar>
class InfiltrationThreePThreeCSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           InfiltrationThreePThreeCSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                     InfiltrationThreePThreeCSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using ThreePhasePcKrSw = FluidMatrix::ParkerVanGenuchten3PDefault<Scalar>;

public:
   using PermeabilityType = Scalar;

    InfiltrationThreePThreeCSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        // intrinsic permeabilities
        fineK_ = 1.e-11;
        coarseK_ = 1.e-11;

        // porosities
        porosity_ = 0.40;
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
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The position for which the porosity is evaluated
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
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

    /*!
     * \brief Returns the temperature at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10.0; // -> 10 degrees Celsius
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return
            70.0 <= globalPos[0] && globalPos[0] <= 85.0 &&
            7.0 <= globalPos[1] && globalPos[1] <= 7.50;
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    const ThreePhasePcKrSw pcKrSwCurve_;
};

} // end namespace Dumux

#endif
