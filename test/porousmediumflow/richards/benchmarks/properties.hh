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
 * \ingroup RichardsTests
 * \brief Common properties for the Richards benchmarks
 */
#ifndef DUMUX_RICHARDS_BENCHMARKS_PROPERTIES_HH
#define DUMUX_RICHARDS_BENCHMARKS_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsBenchmark { using InheritsFrom = std::tuple<Richards>; };
struct RichardsBenchmarkCCTpfa { using InheritsFrom = std::tuple<RichardsBenchmark, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsBenchmark>
{ using type = Dune::YaspGrid<1, Dune::TensorProductCoordinates<double, 1>>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsBenchmark>
{ using type = RichardsBenchmarkProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsBenchmark>
{
    using type = RichardsBenchmarkSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>
    >;
};

} // end namespace Dumux::Properties

#endif
