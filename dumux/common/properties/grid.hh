// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Properties
 * \brief Defines a type tags and some fundamental grid-related properties
 */
#ifndef DUMUX_GRID_PROPERTIES_HH
#define DUMUX_GRID_PROPERTIES_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/pointsource.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {
namespace Properties {

namespace TTag {
//! Type tag for numeric models.
struct GridProperties {};
}

//! Use the minimal point source implementation as default
template<class TypeTag>
struct PointSource<TypeTag, TTag::GridProperties>
{
private:
    using SourceValues = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = typename Dune::FieldVector<typename GridView::ctype, GridView::dimensionworld>;
public:
    using type = Dumux::PointSource<GlobalPosition, SourceValues>;
};

//! Use the point source helper using the bounding box tree as a default
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::GridProperties> { using type = BoundingBoxTreePointSourceHelper; };

} // namespace Properties
} // namespace Dumux

#endif
