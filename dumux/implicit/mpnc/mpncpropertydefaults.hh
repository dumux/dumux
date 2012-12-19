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

#include "mpncindices.hh"

#include "mpncmodel.hh"
#include "mpncindices.hh"
#include "mpnclocalresidual.hh"
#include "mpncfluxvariables.hh"
#include "mpncvolumevariables.hh"
#include "mpncproperties.hh"
#include "mpncnewtoncontroller.hh"
#include "mpncvtkwritermodule.hh"
#include "mpncvtkwritercommon.hh"
#include "mass/mpncvtkwritermass.hh"
#include "energy/mpncvtkwriterenergy.hh"

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/material/spatialparams/boxspatialparams.hh>


/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxMpNcModel
 * \file
 * \brief  Default properties for the Mp-Nc box model.
 */
namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// default property values
//////////////////////////////////////////////////////////////////


//! Set the default pressure formulation to the pressure of the (most) wetting phase
SET_INT_PROP(BoxMPNC,
             PressureFormulation,
             MpNcPressureFormulation::mostWettingFirst);

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system.
 */
SET_PROP(BoxMPNC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const unsigned int value = FluidSystem::numComponents;
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxMPNC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const unsigned int value = FluidSystem::numPhases;
};

/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_PROP(BoxMPNC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    static const unsigned int value = Indices::NumPrimaryVars;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxMPNC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

/*!
 * \brief Set the thermodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 *
 *        \deprecated version. Use "ConstraintSolver" instead of "CompositionFromFugacitiesSolver"
 */
SET_PROP(BoxMPNC, CompositionFromFugacitiesSolver) // DEPRECATED version. Use "ConstraintSolver" instead of "CompositionFromFugacitiesSolver"
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> type; // DEPRECATED version. Use "ConstraintSolver" instead of "CompositionFromFugacitiesSolver"
};

/*!
 * \brief Set the thermodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(BoxMPNC, ConstraintSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef typename GET_PROP_TYPE(TypeTag, CompositionFromFugacitiesSolver)  type;
};


//! Use the MpNc local jacobian operator for the MpNc model
SET_TYPE_PROP(BoxMPNC,
              LocalResidual,
              MPNCLocalResidual<TypeTag>);

//! Use the MpNc specific newton controller for the MpNc model
SET_PROP(BoxMPNC, NewtonController)
{public:
    typedef MPNCNewtonController<TypeTag> type;
};

//! the Model property
SET_TYPE_PROP(BoxMPNC, Model, MPNCModel<TypeTag>);

//! use an isothermal model by default
SET_BOOL_PROP(BoxMPNC, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(BoxMPNC, EnableDiffusion, false);

//! disable kinetic mass transfer by default
SET_BOOL_PROP(BoxMPNC, EnableKinetic, false);

//! disable kinetic energy transfer by default
SET_BOOL_PROP(BoxMPNC, EnableKineticEnergy, false);

//! enable smooth upwinding by default
SET_BOOL_PROP(BoxMPNC, ImplicitEnableSmoothUpwinding, false);

//! the VolumeVariables property
SET_TYPE_PROP(BoxMPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxMPNC, FluxVariables, MPNCFluxVariables<TypeTag>);

//! the Base of the Fluxvariables property (Darcy / Forchheimer)
SET_TYPE_PROP(BoxMPNC, BaseFluxVariables, BoxDarcyFluxVariables<TypeTag>);

//! truncate the newton update for the first few Newton iterations of a time step
SET_BOOL_PROP(BoxMPNC, NewtonEnableChop, true);

//! The indices required by the mpnc model
SET_TYPE_PROP(BoxMPNC, Indices, MPNCIndices<TypeTag, 0>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(BoxMPNC, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(BoxMPNC, ImplicitMobilityUpwindWeight, 1.0);

//! The spatial parameters to be employed. 
//! Use BoxSpatialParams by default.
SET_TYPE_PROP(BoxMPNC, SpatialParams, BoxSpatialParams<TypeTag>);

//! The VTK writer module for common quantities
SET_PROP(BoxMPNC, MPNCVtkCommonModule)
{
    typedef MPNCVtkWriterCommon<TypeTag> type;
};

//! The VTK writer module for quantities which are specific for each
//! mass module
SET_PROP(BoxMPNC, MPNCVtkMassModule)
{
private: enum { enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic) };
public: typedef MPNCVtkWriterMass<TypeTag, enableKinetic> type;
};

//! The VTK writer module for quantities which are specific for each
//! energy module
SET_PROP(BoxMPNC, MPNCVtkEnergyModule)
{
private:
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableKineticEnergy = GET_PROP_VALUE(TypeTag, EnableKineticEnergy) };
public:
    typedef MPNCVtkWriterEnergy<TypeTag, enableEnergy, enableKineticEnergy> type;
};

//! The VTK writer for user specified data (does nothing by default)
SET_PROP(BoxMPNC, MPNCVtkCustomModule)
{ typedef MPNCVtkWriterModule<TypeTag> type; };

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(BoxMPNC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);


//!< Should the averaging of velocities be done in the Model? (By default in the output)
SET_BOOL_PROP(BoxMPNC, VelocityAveragingInModel, false);

//! Specify what to add to the VTK output by default
SET_BOOL_PROP(BoxMPNC, VtkAddPorosity, true);
SET_BOOL_PROP(BoxMPNC, VtkAddPermeability, false);
SET_BOOL_PROP(BoxMPNC, VtkAddBoundaryTypes, false);
SET_BOOL_PROP(BoxMPNC, VtkAddSaturations, true);
SET_BOOL_PROP(BoxMPNC, VtkAddPressures, true);
SET_BOOL_PROP(BoxMPNC, VtkAddVarPressures, false);
SET_BOOL_PROP(BoxMPNC, VtkAddVelocities, false);
SET_BOOL_PROP(BoxMPNC, VtkAddDensities, true);
SET_BOOL_PROP(BoxMPNC, VtkAddMobilities, true);
SET_BOOL_PROP(BoxMPNC, VtkAddAverageMolarMass, false);
SET_BOOL_PROP(BoxMPNC, VtkAddMassFractions, false);
SET_BOOL_PROP(BoxMPNC, VtkAddMoleFractions, true);
SET_BOOL_PROP(BoxMPNC, VtkAddMolarities, false);
SET_BOOL_PROP(BoxMPNC, VtkAddFugacities, false);
SET_BOOL_PROP(BoxMPNC, VtkAddFugacityCoefficients, false);
SET_BOOL_PROP(BoxMPNC, VtkAddTemperatures, false);
SET_BOOL_PROP(BoxMPNC, VtkAddEnthalpies, true);
SET_BOOL_PROP(BoxMPNC, VtkAddInternalEnergies, false);
SET_BOOL_PROP(BoxMPNC, VtkAddReynolds, false);
SET_BOOL_PROP(BoxMPNC, VtkAddPrandtl, false);
SET_BOOL_PROP(BoxMPNC, VtkAddNusselt, false);
SET_BOOL_PROP(BoxMPNC, VtkAddInterfacialArea, false);
SET_BOOL_PROP(BoxMPNC, VtkAddxEquil, false);

// enable gravity by default
SET_BOOL_PROP(BoxMPNC, ProblemEnableGravity, true);

}

}

#endif
