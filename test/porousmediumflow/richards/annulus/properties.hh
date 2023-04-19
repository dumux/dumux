// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

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

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsAnnulus>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} // end namespace Dumux::Properties

#endif
