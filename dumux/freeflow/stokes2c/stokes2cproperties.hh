/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxStokes2cModel
 *
 * \file
 *
 * \brief Defines the properties required for the compositional
 * stokes BOX model.
 */
#ifndef DUMUX_STOKES2CPROPERTIES_HH
#define DUMUX_STOKES2CPROPERTIES_HH

#include "stokes2cindices.hh"

#include <dumux/freeflow/stokes/stokesproperties.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "stokes2clocalresidual.hh"

namespace Dumux
{
/*!
 * \addtogroup BoxStokes2cModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class Stokes2cModel;

template<class TypeTag>
class Stokes2cLocalResidual;

template <class TypeTag>
class Stokes2cVolumeVariables;

template <class TypeTag>
class Stokes2cFluxVariables;

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the compositional stokes problems
NEW_TYPE_TAG(BoxStokes2c, INHERITS_FROM(BoxStokes));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Stokes2cIndices); //!< Enumerations for the compositional stokes models
NEW_PROP_TAG(NumComponents); //!< Number of components

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(BoxStokes2c, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
 public:
    static constexpr int value = 2 + dim;
};

//! Use the stokes2c local jacobian operator for the compositional stokes model
SET_TYPE_PROP(BoxStokes2c,
              LocalResidual,
              Stokes2cLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxStokes2c, Model, Stokes2cModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokes2c, VolumeVariables, Stokes2cVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokes2c, FluxVariables, Stokes2cFluxVariables<TypeTag>);

//! the Indices for the Stokes2c model
SET_TYPE_PROP(BoxStokes2c,
              Stokes2cIndices,
              Stokes2cCommonIndices<TypeTag>);

//! Set the number of components to 2
SET_INT_PROP(BoxStokes2c, NumComponents, 2);

//! Choose the type of the employed fluid state
SET_PROP(BoxStokes2c, FluidState)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;

};

//! Choose the considered phase
SET_PROP(BoxStokes2c, PhaseIndex)
{
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;
public:
    static constexpr int value = Indices::gPhaseIdx;
};

}

}
#endif
