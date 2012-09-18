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
#ifndef DUMUX_MPNC_PROPERTIES_HH
#define DUMUX_MPNC_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>

/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxMpNcModel
 * \file
 * \brief  Defines the properties required for the Mp-Nc box model.
 */
namespace Dumux
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

/*!
 * \brief Define the type tag for the compositional twophase box model.
 */
NEW_TYPE_TAG(BoxMPNC, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(MPNCIndices); //!< DEPRECATED Enumerations for the 2pNc model
NEW_PROP_TAG(MPNCEnergyIndices); //!< DEPRECATED Enumerations for the 2pNc model
NEW_PROP_TAG(Indices); //!< Enumerations for the model

NEW_PROP_TAG(BaseFluxVariables); //!< The type of velocity calculation that is to be used

NEW_PROP_TAG(MPNCVtkCommonModule); //!< Vtk writer module for writing the common quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkMassModule); //!< Vtk writer module for writing the mass related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkEnergyModule); //!< Vtk writer module for writing the energy related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkCustomModule); //!< Vtk writer module for writing the user-specified quantities into the VTK output file

NEW_PROP_TAG(VelocityAveragingInModel);//!< Should the averaging of velocities be done in the model?

//! specify which quantities are written to the vtk output files
NEW_PROP_TAG(VtkAddPorosity);
NEW_PROP_TAG(VtkAddBoundaryTypes);
NEW_PROP_TAG(VtkAddSaturations);
NEW_PROP_TAG(VtkAddPressures);
NEW_PROP_TAG(VtkAddVarPressures);
NEW_PROP_TAG(VtkAddVelocities);
NEW_PROP_TAG(VtkAddDensities);
NEW_PROP_TAG(VtkAddMobilities);
NEW_PROP_TAG(VtkAddAverageMolarMass);
NEW_PROP_TAG(VtkAddMassFractions);
NEW_PROP_TAG(VtkAddMoleFractions);
NEW_PROP_TAG(VtkAddMolarities);
NEW_PROP_TAG(VtkAddFugacities);
NEW_PROP_TAG(VtkAddFugacityCoefficients);
NEW_PROP_TAG(VtkAddTemperatures);
NEW_PROP_TAG(VtkAddEnthalpies);
NEW_PROP_TAG(VtkAddInternalEnergies);

NEW_PROP_TAG(VtkAddxEquil);

NEW_PROP_TAG(VtkAddReynolds);
NEW_PROP_TAG(VtkAddPrandtl);
NEW_PROP_TAG(VtkAddNusselt);
NEW_PROP_TAG(VtkAddInterfacialArea);

NEW_PROP_TAG(MPNCVtkAddPorosity);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddBoundaryTypes);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddSaturations);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddPressures);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddVarPressures);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddVelocities);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddDensities);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddMobilities);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddAverageMolarMass);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddMassFractions);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddMoleFractions);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddMolarities);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddFugacities);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddFugacityCoefficients);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddTemperatures);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddEnthalpies);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddInternalEnergies);//DEPRECATED

NEW_PROP_TAG(MPNCVtkAddxEquil);//DEPRECATED

NEW_PROP_TAG(MPNCVtkAddReynolds);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddPrandtl);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddNusselt);//DEPRECATED
NEW_PROP_TAG(MPNCVtkAddInterfacialArea);//DEPRECATED

NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(SpatialParameters); //!< DEPRECATED The type of the spatial parameters

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the soil)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the soil)

//! The compositional twophase system of fluids which is considered
NEW_PROP_TAG(FluidSystem);

//! The themodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
NEW_PROP_TAG(CompositionFromFugacitiesSolver);

//! Enable the energy equation?
NEW_PROP_TAG(EnableEnergy);

//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! Enable kinetic resolution of mass transfer processes?
NEW_PROP_TAG(EnableKinetic);

//! Enable kinetic resolution of energy transfer processes?
NEW_PROP_TAG(EnableKineticEnergy);

//! Enable gravity?
NEW_PROP_TAG(ProblemEnableGravity);
NEW_PROP_TAG(EnableGravity);//DEPRECATED

//! Use the smooth upwinding method?
NEW_PROP_TAG(ImplicitEnableSmoothUpwinding);
NEW_PROP_TAG(EnableSmoothUpwinding);//DEPRECATED

NEW_PROP_TAG(MassUpwindWeight); //!< DEPRECATED The value of the weight of the upwind direction in the mass conservation equations
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the weight of the upwind direction in the mass conservation equations

NEW_PROP_TAG(MobilityUpwindWeight); //!< DEPRECATED Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation

//! Chop the Newton update at the beginning of the non-linear solver?
NEW_PROP_TAG(NewtonEnableChop);

//! Which type of fluidstate should be used?
NEW_PROP_TAG(FluidState);
}
}

#endif
