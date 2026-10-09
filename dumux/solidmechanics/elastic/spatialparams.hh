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

#include <dumux/common/spatialparams.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup SpatialParameters
 * \ingroup Elastic
 * \brief The base class for spatial parameters of linear elastic problems assembled with
 *        local degrees of freedom (`Experimental::Assembler`). The implementation provides the
 *        Lamé parameters, `lameParams(elemDisc, ipData)` or `lameParamsAtPos(globalPos)`.
 */
template<class GridDiscretization, class Scalar, class Implementation>
class ElasticSpatialParams
: public SpatialParams<GridDiscretization, Scalar, Implementation>
{
    using ParentType = SpatialParams<GridDiscretization, Scalar, Implementation>;
    using ElementDiscretization = typename GridDiscretization::LocalView;

public:
    using ParentType::ParentType;

    //! The Lamé parameters at an interpolation point
    template<class IpData>
    decltype(auto) lameParams(const ElementDiscretization&, const IpData& ipData) const
    { return this->asImp_().lameParamsAtPos(ipData.global()); }
};

} // end namespace Dumux::Experimental

#endif
