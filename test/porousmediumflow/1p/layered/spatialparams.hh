// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial params of the single-phase test with two layers of different permeability
 */
#ifndef DUMUX_ONEP_LAYERED_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_LAYERED_TEST_SPATIALPARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Isotropic permeability that is constant per element and jumps at a plane x = const
 *
 * The permeability is defined per element, so that the local dofs on the plane, which the elements
 * on both sides share, take the permeability of the element the variables belong to.
 */
template<class GridGeometry, class Scalar>
class OnePLayeredSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, OnePLayeredSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePLayeredSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using ElementDiscretization = typename GridGeometry::LocalView;

public:
    using PermeabilityType = Scalar;

    OnePLayeredSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , interfacePosition_(getParam<Scalar>("SpatialParams.InterfacePosition"))
    , permeabilityLeft_(getParam<Scalar>("SpatialParams.PermeabilityLeft"))
    , permeabilityRight_(getParam<Scalar>("SpatialParams.PermeabilityRight"))
    {}

    template<class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization& elemDisc,
                                  const IpData& ipData,
                                  const ElementSolution& elemSol) const
    { return elemDisc.element().geometry().center()[0] < interfacePosition_ ? permeabilityLeft_ : permeabilityRight_; }

    template<class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization& elemDisc,
                    const IpData& ipData,
                    const ElementSolution& elemSol) const
    { return 0.4; }

    Scalar interfacePosition() const
    { return interfacePosition_; }

    Scalar permeabilityLeft() const
    { return permeabilityLeft_; }

    Scalar permeabilityRight() const
    { return permeabilityRight_; }

private:
    Scalar interfacePosition_;
    Scalar permeabilityLeft_;
    Scalar permeabilityRight_;
};

} // end namespace Dumux

#endif
