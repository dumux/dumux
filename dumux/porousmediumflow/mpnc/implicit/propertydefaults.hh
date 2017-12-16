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
#ifndef DUMUX_MPNC_PROPERTY_DEFAULTS_HH
#define DUMUX_MPNC_PROPERTY_DEFAULTS_HH

#include "indices.hh"
#include "model.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "vtkwritermodule.hh"
#include "vtkwritercommon.hh"
#include "mass/vtkwriter.hh"
#include "energy/vtkwriter.hh"

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/material/spatialparams/implicit.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>


/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup BoxMpNcModel
 * \file
 * \brief  Default properties for the MpNc fully implicit model.
 */
namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// default property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system.
 */
SET_PROP(MPNC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const unsigned int value = FluidSystem::numComponents;
};

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_PROP(MPNC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const unsigned int value = FluidSystem::numPhases;
};

/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_PROP(MPNC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    static const unsigned int value = Indices::numPrimaryVars;
};

/*!
 * \brief Set the thermodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(MPNC, ConstraintSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef CompositionFromFugacities<Scalar, FluidSystem> type;
};


//! Use the MpNc local jacobian operator for the MpNc model
SET_TYPE_PROP(MPNC,
              LocalResidual,
              MPNCLocalResidual<TypeTag>);

//! Use the MpNc specific newton controller for the MpNc model
SET_PROP(MPNC, NewtonController)
{public:
    typedef MPNCNewtonController<TypeTag> type;
};

//! the Model property
SET_TYPE_PROP(MPNC, Model, MPNCModel<TypeTag>);

//! use an isothermal model by default
SET_BOOL_PROP(MPNC, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(MPNC, EnableDiffusion, false);

//! disable kinetic mass transfer by default
SET_BOOL_PROP(MPNC, EnableKinetic, false);

//! disable kinetic energy transfer by default
SET_INT_PROP(MPNC, NumEnergyEquations, 0);

//! disable Maxwell Diffusion by default: use Fick
SET_BOOL_PROP(MPNC, UseMaxwellDiffusion, false);

//! enable smooth upwinding by default
SET_BOOL_PROP(MPNC, ImplicitEnableSmoothUpwinding, false);

//! the VolumeVariables property
SET_TYPE_PROP(MPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(MPNC, FluxVariables, MPNCFluxVariables<TypeTag>);

//! the Base of the Fluxvariables property (Darcy / Forchheimer)
SET_TYPE_PROP(MPNC, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! truncate the newton update for the first few Newton iterations of a time step
SET_BOOL_PROP(MPNC, NewtonEnableChop, true);

//! The indices required by the mpnc model
SET_TYPE_PROP(MPNC, Indices, MPNCIndices<TypeTag, 0>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(MPNC, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(MPNC, ImplicitMobilityUpwindWeight, 1.0);

//! The spatial parameters to be employed.
//! Use FVSpatialParams by default.
SET_TYPE_PROP(MPNC, SpatialParams, FVSpatialParams<TypeTag>);

//! The VTK writer module for common quantities
SET_PROP(MPNC, MPNCVtkCommonModule)
{
    typedef MPNCVtkWriterCommon<TypeTag> type;
};

//! The VTK writer module for quantities which are specific for each
//! mass module
SET_PROP(MPNC, MPNCVtkMassModule)
{
private: enum { enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic) };
public: typedef MPNCVtkWriterMass<TypeTag, enableKinetic> type;
};

//! The VTK writer module for quantities which are specific for each
//! energy module
SET_PROP(MPNC, MPNCVtkEnergyModule)
{
private:
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations) };
public:
    typedef MPNCVtkWriterEnergy<TypeTag, enableEnergy, numEnergyEquations> type;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(MPNC, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivitySomerton<Scalar> type;
};

//! The VTK writer for user specified data (does nothing by default)
SET_PROP(MPNC, MPNCVtkCustomModule)
{ typedef MPNCVtkWriterModule<TypeTag> type; };

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(MPNC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Set the default pressure formulation to the pressure of the (most) wetting phase
SET_INT_PROP(MPNC,
             PressureFormulation,
             MpNcPressureFormulation::mostWettingFirst);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(MPNC, SpatialParamsForchCoeff, 0.55);


//!< Should the averaging of velocities be done in the Model? (By default in the output)
SET_BOOL_PROP(MPNC, VelocityAveragingInModel, false);

//! Specify what to add to the VTK output by default
SET_BOOL_PROP(MPNC, VtkAddPorosity, true);
SET_BOOL_PROP(MPNC, VtkAddPermeability, false);
SET_BOOL_PROP(MPNC, VtkAddBoundaryTypes, false);
SET_BOOL_PROP(MPNC, VtkAddSaturations, true);
SET_BOOL_PROP(MPNC, VtkAddPressures, true);
SET_BOOL_PROP(MPNC, VtkAddVarPressures, false);
SET_BOOL_PROP(MPNC, VtkAddVelocities, false);
SET_BOOL_PROP(MPNC, VtkAddDensities, true);
SET_BOOL_PROP(MPNC, VtkAddMobilities, true);
SET_BOOL_PROP(MPNC, VtkAddAverageMolarMass, false);
SET_BOOL_PROP(MPNC, VtkAddMassFractions, false);
SET_BOOL_PROP(MPNC, VtkAddMoleFractions, true);
SET_BOOL_PROP(MPNC, VtkAddMolarities, false);
SET_BOOL_PROP(MPNC, VtkAddFugacities, false);
SET_BOOL_PROP(MPNC, VtkAddFugacityCoefficients, false);
SET_BOOL_PROP(MPNC, VtkAddTemperatures, false);
SET_BOOL_PROP(MPNC, VtkAddEnthalpies, true);
SET_BOOL_PROP(MPNC, VtkAddInternalEnergies, false);
SET_BOOL_PROP(MPNC, VtkAddReynolds, false);
SET_BOOL_PROP(MPNC, VtkAddPrandtl, false);
SET_BOOL_PROP(MPNC, VtkAddNusselt, false);
SET_BOOL_PROP(MPNC, VtkAddInterfacialArea, false);
SET_BOOL_PROP(MPNC, VtkAddxEquil, false);

// enable gravity by default
SET_BOOL_PROP(MPNC, ProblemEnableGravity, true);

}

}

#endif
