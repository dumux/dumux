// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial params of the non-isothermal single-phase test with two layers of different thermal conductivity
 */
#ifndef DUMUX_ONEPNI_LAYERED_SPATIALPARAMS_HH
#define DUMUX_ONEPNI_LAYERED_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Porosity that jumps at a plane x = const, so that the effective thermal conductivity
 *        \f$ \phi \lambda_f + (1 - \phi) \lambda_s \f$ jumps as well
 *
 * The porosity of a local dof is the one at the position of the local dof.
 */
template<class GridGeometry, class Scalar>
class OnePNILayeredSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             OnePNILayeredSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePNILayeredSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    OnePNILayeredSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , interfacePosition_(getParam<Scalar>("SpatialParams.InterfacePosition"))
    , porosityLeft_(getParam<Scalar>("SpatialParams.PorosityLeft"))
    , porosityRight_(getParam<Scalar>("SpatialParams.PorosityRight"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return globalPos[0] < interfacePosition_ ? porosityLeft_ : porosityRight_; }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return permeabilityAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return porosityAtPos(ipData.global()); }

    using ParentType::permeability;
    using ParentType::porosity;

    Scalar interfacePosition() const
    { return interfacePosition_; }

    Scalar porosityLeft() const
    { return porosityLeft_; }

    Scalar porosityRight() const
    { return porosityRight_; }

private:
    Scalar permeability_;
    Scalar interfacePosition_;
    Scalar porosityLeft_;
    Scalar porosityRight_;
};

} // end namespace Dumux

#endif
