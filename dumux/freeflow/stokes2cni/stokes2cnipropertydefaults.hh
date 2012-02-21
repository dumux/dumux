// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \ingroup BoxStokes2cniModel
 *
 * \file
 *
 * \brief Sets default properties for the non-isothermal compositional
 *        Stokes box model.
 */
#ifndef DUMUX_STOKES2CNI_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKES2CNI_PROPERTY_DEFAULTS_HH


#include "stokes2cnifluxvariables.hh"
#include "stokes2cniindices.hh"
#include "stokes2cnilocalresidual.hh"
#include "stokes2cnimodel.hh"
#include "stokes2cnivolumevariables.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(BoxStokes2cni, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
 public:
    static constexpr int value = 3 + dim;
};

//! Use the stokes2cni local jacobian operator for the compositional stokes model
SET_TYPE_PROP(BoxStokes2cni,
              LocalResidual,
              Stokes2cniLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxStokes2cni, Model, Stokes2cniModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokes2cni, VolumeVariables, Stokes2cniVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokes2cni, FluxVariables, Stokes2cniFluxVariables<TypeTag>);

// the indices for the Stokes2cni model
SET_TYPE_PROP(BoxStokes2cni,
              Stokes2cniIndices,
              Stokes2cniCommonIndices<TypeTag>);
}
}
#endif
