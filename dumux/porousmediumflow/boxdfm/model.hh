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
