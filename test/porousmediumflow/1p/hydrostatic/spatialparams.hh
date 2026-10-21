// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Spatial parameters for the hydrostatic one-phase test (element-wise constant permeability)
 */

#ifndef DUMUX_TEST_ONEP_HYDROSTATIC_SPATIAL_PARAMS_HH
#define DUMUX_TEST_ONEP_HYDROSTATIC_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Spatial parameters for the hydrostatic one-phase test.
 *        The permeability jumps between the elements left and right of x = 1.
 */
template<class GridGeometry, class Scalar>
class OnePHydrostaticSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             OnePHydrostaticSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePHydrostaticSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    OnePHydrostaticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeabilityLeft_ = getParam<Scalar>("SpatialParams.PermeabilityLeft");
        permeabilityRight_ = getParam<Scalar>("SpatialParams.PermeabilityRight");
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return element.geometry().center()[0] < 1.0 ? permeabilityLeft_ : permeabilityRight_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

private:
    Scalar permeabilityLeft_, permeabilityRight_;
};

} // end namespace Dumux

#endif
