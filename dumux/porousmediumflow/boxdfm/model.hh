// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \brief Defines a type tag and some properties for porous medium
 *        flow models using the box scheme extended to discrete fractures.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_MODEL_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_MODEL_HH

#include <dumux/discretization/box.hh>

#include "fvgridgeometry.hh"
#include "fluxvariablescache.hh"

namespace Dumux {
namespace Properties {

//! Type tag for the box scheme.
// Create new type tags
namespace TTag {
struct BoxDfmModel { using InheritsFrom = std::tuple<BoxModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::BoxDfmModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = BoxDfmFVGridGeometry<Scalar, GridView, enableCache>;
};

//! The flux variables cache class specific to box-dfm porous medium flow models
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::BoxDfmModel> { using type = BoxDfmFluxVariablesCache<TypeTag>; };

} // namespace Properties
} // namespace Dumux

#endif
