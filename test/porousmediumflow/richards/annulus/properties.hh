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

#include <dune/grid/yaspgrid.hh>

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
{ using type = Dune::YaspGrid<1, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 1>>; };

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

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::RichardsAnnulus>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    // We take the default traits as basis and exchange the extrusion type
    // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

public:
    // Pass the above traits to the box grid geometry such that it uses the
    // rotation-symmetric sub-control volumes and faces.
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

} // end namespace Dumux::Properties

#endif
