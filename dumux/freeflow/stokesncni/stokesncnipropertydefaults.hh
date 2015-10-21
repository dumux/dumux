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
 * \ingroup BoxStokesncniModel
 *
 * \file
 *
 * \brief Sets default properties for the non-isothermal compositional
 *        n-component Stokes box model.
 */
#ifndef DUMUX_STOKESNCNI_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKESNCNI_PROPERTY_DEFAULTS_HH


#include "stokesncnifluxvariables.hh"
#include "stokesncniindices.hh"
#include "stokesncnilocalresidual.hh"
#include "stokesncnimodel.hh"
#include "stokesncnivolumevariables.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(BoxStokesncni, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    static const int dim = Grid::dimension;
 public:
//! for salinization applications only major components exist in the free-flow subdomain
#if SALINIZATION
    static constexpr int value = dim + FluidSystem::numMajorComponents +/*energyequation*/1;
#else
    static constexpr int value = dim + FluidSystem::numComponents +/*energyequation*/1;
#endif
};

//! Use the stokesncni local jacobian operator for the compositional stokes model
SET_TYPE_PROP(BoxStokesncni, LocalResidual, StokesncniLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxStokesncni, Model, StokesncniModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokesncni, VolumeVariables, StokesncniVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokesncni, FluxVariables, StokesncniFluxVariables<TypeTag>);

//! Set the indices for the Stokesncni model
SET_TYPE_PROP(BoxStokesncni, Indices,  StokesncniCommonIndices<TypeTag>);
}
}
#endif
