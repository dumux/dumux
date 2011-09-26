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
#ifndef DUMUX_MPNC_PROPERTY_DEFAULTS_HH
#define DUMUX_MPNC_PROPERTY_DEFAULTS_HH

#include "MpNcindices.hh"

#include "MpNcmodel.hh"

#include "MpNcproblem.hh"
#include "MpNcindices.hh"
#include "MpNclocalresidual.hh"
#include "MpNcfluxvariables.hh"
#include "MpNcvolumevariables.hh"
#include "MpNcproperties.hh"
#include "MpNcnewtoncontroller.hh"

#include "MpNcvtkwritermodule.hh"
#include "MpNcvtkwritercommon.hh"
#include "mass/MpNcvtkwritermass.hh"
#include "energy/MpNcvtkwriterenergy.hh"

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
};

/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_PROP(BoxMPNC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

public:
    static const int value = Indices::NumPrimaryVars;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxMPNC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

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
SET_BOOL_PROP(BoxMPNC, EnableSmoothUpwinding, true);

//! the VolumeVariables property
SET_TYPE_PROP(BoxMPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxMPNC, FluxVariables, MPNCFluxVariables<TypeTag>);

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxMPNC, EnableJacobianRecycling, false);
// enable partial reassembling by default
SET_BOOL_PROP(BoxMPNC, EnablePartialReassemble, true);
// truncate the newton update in the beginning
SET_BOOL_PROP(BoxMPNC, NewtonEnableChop, true);


//! The indices required by the compositional twophase model
SET_PROP(BoxMPNC, MPNCIndices)
{
    typedef MPNCIndices<TypeTag, 0> type;
};

//! The VTK writer module for common quantities
SET_PROP(BoxMPNC, MPNCVtkCommonModule)
{
    typedef MPNCVtkWriterCommon<TypeTag> type;
};

//! The VTK writer module for quantities which are specific for each
//! mass module
SET_PROP(BoxMPNC, MPNCVtkMassModule)
{
private: enum { enableKinetic = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)) };
public: typedef MPNCVtkWriterMass<TypeTag, enableKinetic> type;
};

//! The VTK writer module for quantities which are specific for each
//! energy module
SET_PROP(BoxMPNC, MPNCVtkEnergyModule)
{
private:
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)) };
    enum { enableKineticEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)) };
public:
    typedef MPNCVtkWriterEnergy<TypeTag, enableEnergy, enableKineticEnergy> type;
};

//! The VTK writer for user specified data (does nothing by default)
SET_PROP(BoxMPNC, MPNCVtkCustomModule)
{ typedef MPNCVtkWriterModule<TypeTag> type; };


//!< Should the averaging of velocities be done in the problem? (By default in the output)
SET_BOOL_PROP(BoxMPNC, VelocityAveragingInProblem, false);

//! Specify what to add to the VTK output by default
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPorosity, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddBoundaryTypes, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddSaturations, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPressures, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddVarPressures, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddVelocities, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddDensities, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMobilities, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMeanMolarMass, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMassFractions, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMoleFractions, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddMolarities, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddFugacities, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddFugacityCoefficients, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddTemperatures, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddEnthalpies, true);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddInternalEnergies, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddReynolds, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddPrandtl, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddNusselt, false);
SET_BOOL_PROP(BoxMPNC, MPNCVtkAddInterfacialArea, false);
}

}

#endif
