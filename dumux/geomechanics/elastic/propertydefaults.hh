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
 * \ingroup ElasticBoxModel
 * \file
 *
 * \brief Defines some default values for the properties of the
 * linear elasticity model.
 */


#ifndef DUMUX_ELASTIC_PROPERTIES_DEFAULTS_HH
#define DUMUX_ELASTIC_PROPERTIES_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"



namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//!< set the number of equations to the space dimension of the problem
SET_PROP(BoxElastic, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    static const int value = dim;
};

//! Use the linear elasticity local residual function for the elasticity model
SET_TYPE_PROP(BoxElastic,
              LocalResidual,
              ElasticLocalResidual<TypeTag>);

//! define the model
SET_TYPE_PROP(BoxElastic, Model, ElasticModel<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(BoxElastic, VolumeVariables, ElasticVolumeVariablesBase<TypeTag>);

//! define the FluxVariables
SET_TYPE_PROP(BoxElastic, FluxVariables, ElasticFluxVariablesBase<TypeTag>);

//! Set the indices used by the linear elasticity model
SET_TYPE_PROP(BoxElastic, Indices, ElasticIndices<>);

//! enable gravity by default
SET_BOOL_PROP(BoxElastic, ProblemEnableGravity, true);
}
}

#endif
