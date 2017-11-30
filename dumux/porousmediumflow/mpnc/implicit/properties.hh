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

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup BoxMpNcModel
 * \file
 * \brief  Defines the properties required for the MpNc fully implicit model.
 */
namespace Dumux
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit m-phase n-component problems
NEW_TYPE_TAG(MPNC);
NEW_TYPE_TAG(BoxMPNC, INHERITS_FROM(BoxModel, MPNC));
NEW_TYPE_TAG(CCMPNC, INHERITS_FROM(CCModel, MPNC));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model

NEW_PROP_TAG(BaseFluxVariables); //!< The type of velocity calculation that is to be used

NEW_PROP_TAG(PressureFormulation);   //!< The formulation of the model

NEW_PROP_TAG(MPNCVtkCommonModule); //!< Vtk writer module for writing the common quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkMassModule); //!< Vtk writer module for writing the mass related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkEnergyModule); //!< Vtk writer module for writing the energy related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkCustomModule); //!< Vtk writer module for writing the user-specified quantities into the VTK output file

NEW_PROP_TAG(VelocityAveragingInModel);//!< Should the averaging of velocities be done in the model?

//! specify which quantities are written to the vtk output files
NEW_PROP_TAG(VtkAddPorosity);
NEW_PROP_TAG(VtkAddPermeability);
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

NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatialParams)

//! The compositional twophase system of fluids which is considered
NEW_PROP_TAG(FluidSystem);

//! The thermodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
NEW_PROP_TAG(ConstraintSolver);

//! Enable the energy equation?
NEW_PROP_TAG(EnableEnergy);

//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! Enable kinetic resolution of mass transfer processes?
NEW_PROP_TAG(EnableKinetic);

//! Property for the definition of the number of energy equations (0,1,2,3)
NEW_PROP_TAG(NumEnergyEquations);

//! Enable Maxwell Diffusion? (If false: use Fickian Diffusion) Maxwell incorporated the mutual
//! influences of multiple diffusing components. However, Fick seems to be more robust.
NEW_PROP_TAG(UseMaxwellDiffusion);

//! The model for the effective thermal conductivity
NEW_PROP_TAG(ThermalConductivityModel);

//! Enable gravity?
NEW_PROP_TAG(ProblemEnableGravity);

//! Use the smooth upwinding method?
NEW_PROP_TAG(ImplicitEnableSmoothUpwinding);

NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the weight of the upwind direction in the mass conservation equations

NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation

//! Chop the Newton update at the beginning of the non-linear solver?
NEW_PROP_TAG(NewtonEnableChop);

//! Which type of fluidstate should be used?
NEW_PROP_TAG(FluidState);

//! Property for the forchheimer coefficient
NEW_PROP_TAG(SpatialParamsForchCoeff);
}
}

#endif
