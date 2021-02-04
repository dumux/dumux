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
 * \brief Test for the Mandel problem.
 */
#ifndef MANDEL_PROPERTIES_HH
#define MANDEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/poroelastic/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/h2o.hh>

#include "problem.hh"
namespace Dumux {
namespace Properties{

////////////////////////////
//////Fluid Problem/////////
////////////////////////////
namespace TTag {
struct TestMandelFluid { using InheritsFrom = std::tuple<> };
}
////////////////////////////
/////Solid Problem//////////
////////////////////////////
namespace TTag {
// Create new type tags
struct TestMandelSolid { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestMandelSolid> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestMandelSolid> { using type = Dumux::MandelSolidProblem<TypeTag>; };

// The fluid phase consists of water
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestPoroElastic>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::OnePLiquid< Scalar,
                                                  Dumux::Components::h2O<Scalar>
                                                  >;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestPoroElastic>
{
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry> >;
};
}// end properties

}// end namespace
#endif
