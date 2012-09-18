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
#include "mpncproblem.hh"
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
 * \brief Set the themodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(BoxMPNC, CompositionFromFugacitiesSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> type;
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
SET_BOOL_PROP(BoxMPNC, ImplicitEnableSmoothUpwinding, GET_PROP_VALUE(TypeTag, EnableSmoothUpwinding));
SET_BOOL_PROP(BoxMPNC, EnableSmoothUpwinding, false);//DEPRECATED

//! the VolumeVariables property
SET_TYPE_PROP(BoxMPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxMPNC, FluxVariables, MPNCFluxVariables<TypeTag>);

//! the Base of the Fluxvariables property (Darcy / Forchheimer)
SET_TYPE_PROP(BoxMPNC, BaseFluxVariables, BoxDarcyFluxVariables<TypeTag>);

//! truncate the newton update for the first few Newton iterations of a time step
SET_BOOL_PROP(BoxMPNC, NewtonEnableChop, true);

//! The indices required by the mpnc model
SET_PROP(BoxMPNC, Indices)
{
    typedef MPNCIndices<TypeTag, 0> type;
};

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(BoxMPNC, ImplicitMassUpwindWeight, GET_PROP_VALUE(TypeTag, MassUpwindWeight));
SET_SCALAR_PROP(BoxMPNC, MassUpwindWeight, 1.0);//DEPRECATED

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(BoxMPNC, ImplicitMobilityUpwindWeight, GET_PROP_VALUE(TypeTag, MobilityUpwindWeight));
SET_SCALAR_PROP(BoxMPNC, MobilityUpwindWeight, 1.0);//DEPRECATED

//! DEPRECATED MPNCIndices property
SET_TYPE_PROP(BoxMPNC, MPNCIndices, typename GET_PROP_TYPE(TypeTag, Indices));

//! DEPRECATED SpatialParameters property
SET_TYPE_PROP(BoxMPNC, SpatialParameters, typename GET_PROP_TYPE(TypeTag, SpatialParams));

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


//!< Should the averaging of velocities be done in the Model? (By default in the output)
SET_BOOL_PROP(BoxMPNC, VelocityAveragingInModel, false);

//! Specify what to add to the VTK output by default
SET_BOOL_PROP(BoxMPNC, VtkAddPorosity, GET_PROP_VALUE(TypeTag, MPNCVtkAddPorosity));
SET_BOOL_PROP(BoxMPNC, VtkAddBoundaryTypes, GET_PROP_VALUE(TypeTag, MPNCVtkAddBoundaryTypes));
SET_BOOL_PROP(BoxMPNC, VtkAddSaturations, GET_PROP_VALUE(TypeTag, MPNCVtkAddSaturations));
SET_BOOL_PROP(BoxMPNC, VtkAddPressures, GET_PROP_VALUE(TypeTag, MPNCVtkAddPressures));
SET_BOOL_PROP(BoxMPNC, VtkAddVarPressures, GET_PROP_VALUE(TypeTag, MPNCVtkAddVarPressures));
SET_BOOL_PROP(BoxMPNC, VtkAddVelocities, GET_PROP_VALUE(TypeTag, MPNCVtkAddVelocities));
SET_BOOL_PROP(BoxMPNC, VtkAddDensities, GET_PROP_VALUE(TypeTag, MPNCVtkAddDensities));
SET_BOOL_PROP(BoxMPNC, VtkAddMobilities, GET_PROP_VALUE(TypeTag, MPNCVtkAddMobilities));
SET_BOOL_PROP(BoxMPNC, VtkAddAverageMolarMass, GET_PROP_VALUE(TypeTag, MPNCVtkAddAverageMolarMass));
SET_BOOL_PROP(BoxMPNC, VtkAddMassFractions, GET_PROP_VALUE(TypeTag, MPNCVtkAddMassFractions));
SET_BOOL_PROP(BoxMPNC, VtkAddMoleFractions, GET_PROP_VALUE(TypeTag, MPNCVtkAddMoleFractions));
SET_BOOL_PROP(BoxMPNC, VtkAddMolarities, GET_PROP_VALUE(TypeTag, MPNCVtkAddMolarities));
SET_BOOL_PROP(BoxMPNC, VtkAddFugacities, GET_PROP_VALUE(TypeTag, MPNCVtkAddFugacities));
SET_BOOL_PROP(BoxMPNC, VtkAddFugacityCoefficients, GET_PROP_VALUE(TypeTag, MPNCVtkAddFugacityCoefficients));
SET_BOOL_PROP(BoxMPNC, VtkAddTemperatures, GET_PROP_VALUE(TypeTag, MPNCVtkAddTemperatures));
SET_BOOL_PROP(BoxMPNC, VtkAddEnthalpies, GET_PROP_VALUE(TypeTag, MPNCVtkAddEnthalpies));
SET_BOOL_PROP(BoxMPNC, VtkAddInternalEnergies, GET_PROP_VALUE(TypeTag, MPNCVtkAddInternalEnergies));
SET_BOOL_PROP(BoxMPNC, VtkAddReynolds, GET_PROP_VALUE(TypeTag, MPNCVtkAddReynolds));
SET_BOOL_PROP(BoxMPNC, VtkAddPrandtl, GET_PROP_VALUE(TypeTag, MPNCVtkAddPrandtl));
SET_BOOL_PROP(BoxMPNC, VtkAddNusselt, GET_PROP_VALUE(TypeTag, MPNCVtkAddNusselt));
SET_BOOL_PROP(BoxMPNC, VtkAddInterfacialArea, GET_PROP_VALUE(TypeTag, MPNCVtkAddInterfacialArea));
SET_BOOL_PROP(BoxMPNC, VtkAddxEquil, GET_PROP_VALUE(TypeTag, MPNCVtkAddxEquil));

SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPorosity, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddBoundaryTypes, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddSaturations, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPressures, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddVarPressures, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddVelocities, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddDensities, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMobilities, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddAverageMolarMass, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMassFractions, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMoleFractions, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMolarities, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddFugacities, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddFugacityCoefficients, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddTemperatures, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddEnthalpies, true);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddInternalEnergies, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddReynolds, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPrandtl, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddNusselt, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddInterfacialArea, false);//DEPRECATED
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddxEquil, false);//DEPRECATED

//Has to be removed if DEPRECATED EnableGravity is removed!
SET_BOOL_PROP(BoxMPNC, ProblemEnableGravity, GET_PROP_VALUE(TypeTag, EnableGravity));

}

}

#endif
