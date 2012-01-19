/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \file
 *
 * \brief Defines the properties required for the Stokes BOX model.
 */

#ifndef DUMUX_STOKESPROPERTIES_HH
#define DUMUX_STOKESPROPERTIES_HH

#include "stokesindices.hh"

#include <dumux/boxmodels/common/boxproperties.hh>
#include <dumux/freeflow/stokes/stokeslocaljacobian.hh>

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>



namespace Dumux
{
/*!
 * \addtogroup StokesModel
 */
// \{

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
 * \addtogroup StokesModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the stokes problems
NEW_TYPE_TAG(BoxStokes, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//NEW_PROP_TAG(Soil); //!< The type of the soil properties object
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(MassUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(StokesIndices); //!< Enumerations for the Stokes models
NEW_PROP_TAG(Fluid);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(StabilizationAlpha);
NEW_PROP_TAG(StabilizationBeta);

NEW_PROP_TAG(CalculateVariablesAtCenterOfSCV);
NEW_PROP_TAG(PhaseIndex);

NEW_PROP_TAG(SpatialParameters);

NEW_PROP_TAG(Scaling); //!Defines Scaling of the model
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

//// enable jacobian matrix recycling by default
//SET_BOOL_PROP(BoxStokes, EnableJacobianRecycling, false);
//// enable partial reassembling by default
//SET_BOOL_PROP(BoxStokes, EnablePartialReassemble, false);
//// set some Newton properties deviating from the default ones
//
//SET_SCALAR_PROP(BoxStokes, NewtonRelTolerance, 1e-4);
//SET_INT_PROP(BoxStokes, NewtonTargetSteps, 10);

}
}

#endif
