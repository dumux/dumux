// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \ingroup BoxStokesModel
 *
 * \file
 *
 * \brief Defines the properties required for the Stokes BOX model.
 */

#ifndef DUMUX_STOKES_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKES_PROPERTY_DEFAULTS_HH

#include "stokesproperties.hh"
#include "stokesindices.hh"
#include "stokeslocaljacobian.hh"
#include "stokesnewtoncontroller.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class StokesModel;

template<class TypeTag>
class StokesLocalResidual;

template <class TypeTag>
class StokesVolumeVariables;

template <class TypeTag>
class StokesFluxVariables;

template <class TypeTag>
class StokesNewtonController;

// \}

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

/*!
 * \addtogroup BoxStokesModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! The local jacobian operator for the stokes box scheme
SET_TYPE_PROP(BoxStokes, LocalJacobian, Dumux::StokesLocalJacobian<TypeTag>);

SET_PROP(BoxStokes, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
 public:
    static constexpr int value = 1 + dim;
};

SET_SCALAR_PROP(BoxStokes, Scaling, 1); //!< set scaling to 1 by default

//! Use the Stokes local residual function for the Stokes model
SET_TYPE_PROP(BoxStokes,
              LocalResidual,
              StokesLocalResidual<TypeTag>);

// Use the Stokes specific newton controller for the Stokes model
SET_PROP(BoxStokes, NewtonController)
{public:
    typedef StokesNewtonController<TypeTag> type;
};

#if HAVE_SUPERLU
SET_TYPE_PROP(BoxStokes, LinearSolver, SuperLUBackend<TypeTag>);
#endif

//! the Model property
SET_TYPE_PROP(BoxStokes, Model, StokesModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokes, VolumeVariables, StokesVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokes, FluxVariables, StokesFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(BoxStokes, MassUpwindWeight, 1.0);

//! The fluid system to use by default
SET_PROP(BoxStokes, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::OneP<Scalar, Fluid> type;
};

SET_PROP(BoxStokes, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_TYPE_PROP(BoxStokes, StokesIndices, StokesCommonIndices<TypeTag>);

//! Choose the type of the employed fluid state
SET_PROP(BoxStokes, FluidState)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> type;

};

SET_INT_PROP(BoxStokes, PhaseIndex, 0);

// \}

}
}

#endif
