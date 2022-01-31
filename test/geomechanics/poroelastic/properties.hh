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
 * \ingroup GeomechanicsTests
 * \brief The properties of a test problem for the poro-elastic model.
 */
#ifndef DUMUX_POROELASTIC_PROPERTIES_HH
#define DUMUX_POROELASTIC_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct TestPoroElastic { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestPoroElastic> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestPoroElastic> { using type = Dumux::PoroElasticProblem<TypeTag>; };

// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestPoroElastic>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};

// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestPoroElastic>
{
    using FS = GetPropType<TypeTag, Properties::FluidSystem>;
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry>,
                                           FS, PV, Indices>;
};

} // end namespace Dumux::Properties

#endif
