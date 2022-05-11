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

#ifndef DUMUX_RICHARDS_ANNULUS_PROPERTIES_HH
#define DUMUX_RICHARDS_ANNULUS_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/discretization/extrusion.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct RichardsAnnulus { using InheritsFrom = std::tuple<Richards, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsAnnulus>
{ using type =  Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsAnnulus>
{
    using type = RichardsAnnulusProblem<
        TypeTag, GetPropType<TypeTag, Properties::FluidSystem>
    >;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsAnnulus>
{
    using type = RichardsAnnulusSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>
    >;
};

} // end namespace Dumux::Properties

#endif
