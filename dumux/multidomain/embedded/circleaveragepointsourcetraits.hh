// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Point source traits for average-based coupling modes
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLE_AVERAGE_POINT_SOURCE_TRAITS_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLE_AVERAGE_POINT_SOURCE_TRAITS_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>

namespace Dumux {

//! point source traits for the circle average coupling mode
template<class MDTraits>
struct CircleAveragePointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using GridGeometry = GetPropType<SubDomainTypeTag<i>, Properties::GridGeometry>;
    template<std::size_t i> using NumEqVector = Dumux::NumEqVector<GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>>;
public:
    //! export the point source type for domain i
    template<std::size_t i>
    using PointSource = IntegrationPointSource<typename GridGeometry<i>::GlobalCoordinate, NumEqVector<i>>;

    //! export the point source helper type  for domain i
    template<std::size_t i>
    using PointSourceHelper = IntegrationPointSourceHelper;

    //! export the point source data type
    using PointSourceData = PointSourceDataCircleAverage<MDTraits>;
};

} // end namespace Dumux

#endif
