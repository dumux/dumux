// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Properties
 * \ingroup CCTpfaProperties
 * \ingroup CCTpfaModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_CCTPFA_PROPERTY_DEFAULTS_HH
#define DUMUX_CCTPFA_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/implicit/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/implicit/cellcentered/tpfa/fvelementgeometryvector.hh>
#include <dumux/implicit/cellcentered/properties.hh>

namespace Dumux {

// forward declarations
template<class TypeTag> class CCElementBoundaryTypes;

namespace Properties {
//! Set the corresponding discretization method property
SET_INT_PROP(CCTpfaModel, DiscretizationMethod, GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa);

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(CCTpfaModel, FVElementGeometryVector, CCTpfaFVElementGeometryVector<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVElementGeometryCache)>);

//! The sub control volume
SET_PROP(CCTpfaModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::SubControlVolume<ScvGeometry, IndexType, /*isBox=*/false> type;
};

SET_PROP(CCTpfaModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::SubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(CCTpfaModel, FluxVariablesCache, Dumux::CCTpfaPorousMediumFluxVariablesCache<TypeTag>);

} // namespace Properties

} // namespace Dumux

#endif
