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
 * \ingroup BoxStokesModel
 *
 * \file
 *
 * \brief Defines default properties for the Stokes box model.
 *
 * These can be overwritten at a different place or
 * may be replaced by values of the input file.
 */

#ifndef DUMUX_STOKES_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKES_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include "indices.hh"
//#include "localjacobian.hh"
#include "localresidual.hh"
#include "model.hh"
//#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "newtoncontroller.hh"


//added to resemble geomechanics
#include <dumux/discretization/fem/ipdata.hh>  //include ipData from elastic for FEM
#include "secondaryvariables.hh"


#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/fluidsystems/1p.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! The local jacobian operator for the stokes box scheme
//SET_TYPE_PROP(BoxStokes, LocalJacobian, StokesLocalJacobian<TypeTag>);

SET_PROP(BoxStokes, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
 public:
    static constexpr int value = 1 + dim;
};

SET_SCALAR_PROP(BoxStokes, Scaling, 1); //!< set scaling to 1 by default

//! Use the Stokes local residual function for the Stokes model
SET_TYPE_PROP(BoxStokes, LocalResidual, StokesLocalResidual<TypeTag>);

//! Use the Stokes specific newton controller for the Stokes model
SET_TYPE_PROP(BoxStokes, NewtonController, StokesNewtonController<TypeTag>);

#if HAVE_SUPERLU
SET_TYPE_PROP(BoxStokes, LinearSolver, SuperLUBackend<TypeTag>);
#elif HAVE_UMFPACK
SET_TYPE_PROP(BoxStokes, LinearSolver, UMFPackBackend<TypeTag>);
#endif

//! the Model property
SET_TYPE_PROP(BoxStokes, Model, StokesFemModel<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(BoxStokes, SecondaryVariables, StokesSecondaryVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokes, FluxVariables, StokesFluxVariables<TypeTag>);

//! the upwind factor.
SET_SCALAR_PROP(BoxStokes, ImplicitMassUpwindWeight, 1.0);

//! The fluid system to use by default
SET_PROP(BoxStokes, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::OneP<Scalar, Fluid> type;
};

//! The fluid that is used in the single-phase fluidsystem
SET_PROP(BoxStokes, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

//! Set the indices used by the Stokes model
SET_TYPE_PROP(BoxStokes, Indices, StokesCommonIndices<TypeTag>);

//! Choose the type of the employed fluid state.
SET_PROP(BoxStokes, FluidState)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;

};

//! Enable evaluation of shape function gradients at the sub-control volume center by default
// Used for the computation of the pressure gradients
SET_BOOL_PROP(BoxStokes, EvalGradientsAtSCVCenter, true);

//! Set the phaseIndex per default to zero (for two-phase or two-component system imported
//  from fluid systems, see property defaults in stokesnc model).
SET_INT_PROP(BoxStokes, PhaseIdx, 0);

//! Use symmetrizedVelocityGradient by default
SET_BOOL_PROP(BoxStokes, EnableUnsymmetrizedVelocityGradient, false);

//! Disbale additional wall friction term by default
SET_BOOL_PROP(BoxStokes, EnablePseudoThreeDWallFriction, false);

//! Set calculation to Stokes, not Navier-Stokes
SET_BOOL_PROP(BoxStokes, EnableNavierStokes, false);

//! The mass density is used in the continuity equation
SET_BOOL_PROP(BoxStokes, UseMoles, false);

//! A stabilization factor. Set negative for stabilization and to zero for no stabilization
SET_SCALAR_PROP(BoxStokes, StokesStabilizationAlpha, 0.0);

//! Stabilization factor for the boundaries
SET_SCALAR_PROP(BoxStokes, StokesStabilizationBeta, 0.0);

//! Set gravity by default
SET_BOOL_PROP(BoxStokes, ProblemEnableGravity, true);

//! Set the value for the maximum relative shift
SET_SCALAR_PROP(BoxStokes, NewtonMaxRelativeShift, 1e-6);

//! Set the number of iterations at which the Newton method should aim at.
SET_INT_PROP(BoxStokes, NewtonTargetSteps, 10);

//! Set the number of maximum iterations for the Newton method.
SET_INT_PROP(BoxStokes, NewtonMaxSteps, 18);

}
}

#endif
