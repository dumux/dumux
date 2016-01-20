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
 * \ingroup BoxStokesncModel
 *
 * \file
 *
 * \brief Defines the properties required for the compositional
 * Stokes box model.
 */
#ifndef DUMUX_STOKESNC_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKESNC_PROPERTY_DEFAULTS_HH

#include <dumux/freeflow/stokes/propertydefaults.hh>

#include "properties.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "model.hh"
#include "volumevariables.hh"

#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Define that mole fractions are used in the balance equations
SET_BOOL_PROP(BoxStokesnc, UseMoles, true);

//! set the number of equations
SET_PROP(BoxStokesnc, NumEq)
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static const int dim = Grid::dimension;
public:
    static constexpr int value = FluidSystem::numComponents + dim;
};

//! Use the Stokes nc local residual function for the Stokes model
SET_TYPE_PROP(BoxStokesnc, LocalResidual, StokesncLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxStokesnc, Model, StokesncModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokesnc, VolumeVariables, StokesncVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokesnc, FluxVariables, StokesncFluxVariables<TypeTag>);

//! Set the Indices for the Stokes nc model.
SET_TYPE_PROP(BoxStokesnc, Indices, StokesncCommonIndices<TypeTag>);

//! Choose the type of the employed fluid state
SET_PROP(BoxStokesnc, FluidState)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Choose the considered phase (single-phase system); the gas phase is used
SET_PROP(BoxStokesnc, PhaseIdx)
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    static constexpr int value = FluidSystem::nPhaseIdx;
};

} // namespace Properties
} // namespace Dumux
#endif // DUMUX_STOKESNC_PROPERTY_DEFAULTS_HH
