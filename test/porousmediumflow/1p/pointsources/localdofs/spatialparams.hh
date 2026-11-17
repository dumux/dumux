// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Spatial parameters of the point source tests with unit permeability
 */
#ifndef DUMUX_ONEP_POINTSOURCE_LOCALDOFS_SPATIALPARAMS_HH
#define DUMUX_ONEP_POINTSOURCE_LOCALDOFS_SPATIALPARAMS_HH

#include <type_traits>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/porousmediumflow/spatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class OnePPointSourceSpatialParams;

namespace Detail {
//! The finite volume base class for schemes with sub-control volumes, the general one otherwise
template<class GridGeometry, class Scalar>
using OnePPointSourceSpatialParamsBase = std::conditional_t<
    requires { typename GridGeometry::SubControlVolume; },
    FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, OnePPointSourceSpatialParams<GridGeometry, Scalar>>,
    Experimental::PorousMediumFlowSpatialParams<GridGeometry, Scalar, OnePPointSourceSpatialParams<GridGeometry, Scalar>>
>;
} // end namespace Detail

//! Unit permeability and constant porosity
template<class GridGeometry, class Scalar>
class OnePPointSourceSpatialParams
: public Detail::OnePPointSourceSpatialParamsBase<GridGeometry, Scalar>
{
    using ParentType = Detail::OnePPointSourceSpatialParamsBase<GridGeometry, Scalar>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    using ParentType::permeability;
    using ParentType::porosity;

    template<class ElementDiscretization, class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization& elemDisc,
                                  const IpData& ipData,
                                  const ElementSolution& elemSol) const
    { return permeabilityAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization& elemDisc,
                    const IpData& ipData,
                    const ElementSolution& elemSol) const
    { return porosityAtPos(ipData.global()); }
};

} // end namespace Dumux

#endif
