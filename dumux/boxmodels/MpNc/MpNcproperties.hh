// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
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
NEW_PROP_TAG(MPNCIndices); //!< Enumerations for the 2pNc model
NEW_PROP_TAG(MPNCEnergyIndices); //!< Enumerations for the 2pNc model

NEW_PROP_TAG(MPNCVtkCommonModule); //!< Vtk writer module for writing the common quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkMassModule); //!< Vtk writer module for writing the mass related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkEnergyModule); //!< Vtk writer module for writing the energy related quantities into the VTK output file
NEW_PROP_TAG(MPNCVtkCustomModule); //!< Vtk writer module for writing the user-specified quantities into the VTK output file

NEW_PROP_TAG(VelocityAveragingInModel);//!< Should the averaging of velocities be done in the model?

//! specify which quantities are written to the vtk output files
NEW_PROP_TAG(MPNCVtkAddPorosity);
NEW_PROP_TAG(MPNCVtkAddBoundaryTypes);
NEW_PROP_TAG(MPNCVtkAddSaturations);
NEW_PROP_TAG(MPNCVtkAddPressures);
NEW_PROP_TAG(MPNCVtkAddVarPressures);
NEW_PROP_TAG(MPNCVtkAddVelocities);
NEW_PROP_TAG(MPNCVtkAddDensities);
NEW_PROP_TAG(MPNCVtkAddMobilities);
NEW_PROP_TAG(MPNCVtkAddAverageMolarMass);
NEW_PROP_TAG(MPNCVtkAddMassFractions);
NEW_PROP_TAG(MPNCVtkAddMoleFractions);
NEW_PROP_TAG(MPNCVtkAddMolarities);
NEW_PROP_TAG(MPNCVtkAddFugacities);
NEW_PROP_TAG(MPNCVtkAddFugacityCoefficients);
NEW_PROP_TAG(MPNCVtkAddTemperatures);
NEW_PROP_TAG(MPNCVtkAddEnthalpies);
NEW_PROP_TAG(MPNCVtkAddInternalEnergies);

NEW_PROP_TAG(MPNCVtkAddReynolds);
NEW_PROP_TAG(MPNCVtkAddPrandtl);
NEW_PROP_TAG(MPNCVtkAddNusselt);
NEW_PROP_TAG(MPNCVtkAddInterfacialArea);

NEW_PROP_TAG(SpatialParameters); //!< The type of the soil properties object

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
NEW_PROP_TAG(EnableGravity);

//! Use the smooth upwinding method?
NEW_PROP_TAG(EnableSmoothUpwinding);

//! Chop the Newton update at the beginning of the non-linear solver?
NEW_PROP_TAG(NewtonEnableChop);
}
}

#endif
