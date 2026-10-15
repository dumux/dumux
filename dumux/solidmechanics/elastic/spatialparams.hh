// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \ingroup Elastic
 * \brief The base class for spatial parameters of linear elastic problems assembled with
 *        local degrees of freedom
 */
#ifndef DUMUX_SOLIDMECHANICS_ELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_SOLIDMECHANICS_ELASTIC_SPATIAL_PARAMS_HH

#include <memory>

#include <dune/common/exceptions.hh>

#include <dumux/common/spatialparams.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup SpatialParameters
 * \ingroup Elastic
 * \brief The base class for spatial parameters of linear elastic problems assembled with
 *        local degrees of freedom (`Experimental::Assembler`). The implementation provides the
 *        Lamé parameters, `lameParams(elemDisc, ipData)` or `lameParamsAtPos(globalPos)`, and for
 *        problems with inertia the solid density, `solidDensity(elemDisc, ipData)` or
 *        `solidDensityAtPos(globalPos)`.
 */
template<class GridDiscretization, class Scalar, class Implementation>
class ElasticSpatialParams
: public SpatialParams<GridDiscretization, Scalar, Implementation>
{
    using ParentType = SpatialParams<GridDiscretization, Scalar, Implementation>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GlobalPosition = typename GridDiscretization::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

    //! The Lamé parameters at an interpolation point
    template<class IpData>
    decltype(auto) lameParams(const ElementDiscretization&, const IpData& ipData) const
    { return this->asImp_().lameParamsAtPos(ipData.global()); }

    //! The solid density at an interpolation point (required for problems with inertia)
    template<class IpData>
    Scalar solidDensity(const ElementDiscretization&, const IpData& ipData) const
    { return this->asImp_().solidDensityAtPos(ipData.global()); }

    //! The solid density at a position (required for problems with inertia)
    Scalar solidDensityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::InvalidStateException,
            "The spatial parameters of problems with inertia have to provide solidDensity(elemDisc, ipData) or solidDensityAtPos(globalPos)");
    }
};

} // end namespace Dumux::Experimental

#endif
