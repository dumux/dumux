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
 * \ingroup ImplicitProperties
 * \ingroup MimeticModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_MIMETIC_PROPERTY_DEFAULTS_HH
#define DUMUX_MIMETIC_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include <dumux/discretization/methods.hh>

#include <dumux/porousmediumflow/implicit/mimetic/fluxvariables.hh>
#include <dumux/porousmediumflow/implicit/mimetic/velocityoutput.hh>
#include <dumux/discretization/staggered/mimetic/subcontrolvolume.hh>
#include <dumux/discretization/staggered/mimetic/globalfvgeometry.hh>

#include <dumux/discretization/staggered/mimetic/globalfluxvariablescache.hh>
#include <dumux/discretization/staggered/mimetic/mimeticgeometryhelper.hh>
#include <dumux/discretization/staggered/mimetic/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/mimetic/facevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/localresidual.hh>

namespace Dumux {

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(MimeticModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Mimetic;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(MimeticModel, GlobalFVGeometry, MimeticGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The geometry helper required for the stencils, etc.
SET_PROP(MimeticModel, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = MimeticGeometryHelper<GridView>;
};

//! The sub control volume
SET_PROP(MimeticModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::MimeticSubControlVolume<ScvGeometry, IndexType> type;
};

//! The sub-controlvolume face
SET_PROP(MimeticModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::MimeticSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! The variables living on the faces
SET_TYPE_PROP(MimeticModel, FaceVariables, MimeticFaceVariables<TypeTag>);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(MimeticModel, FluxVariables, PorousMediumFluxVariablesMimetic<TypeTag>);

//! The global flux variables cache vector class
SET_TYPE_PROP(MimeticModel, GlobalFluxVariablesCache, Dumux::MimeticGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

SET_TYPE_PROP(MimeticModel, VelocityOutput, MimeticVelocityOutput<TypeTag>);

//! Per default we have assume isothermal problems. Set this to true to solve an energy equation
SET_BOOL_PROP(MimeticModel, EnableEnergyBalance, false);

SET_TYPE_PROP(MimeticModel, EnergyLocalResidual, EnergyLocalResidual<TypeTag> );

SET_BOOL_PROP(MimeticModel, VtkWriteFaceData, false);
} // namespace Properties

} // namespace Dumux

#endif
