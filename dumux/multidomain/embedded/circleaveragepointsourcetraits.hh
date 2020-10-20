// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Point source traits for average-based coupling modes
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLE_AVERAGE_POINT_SOURCE_TRAITS_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLE_AVERAGE_POINT_SOURCE_TRAITS_HH

#include <dumux/common/properties.hh>
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
    template<std::size_t i> using NumEqVector = GetPropType<SubDomainTypeTag<i>, Properties::NumEqVector>;
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
