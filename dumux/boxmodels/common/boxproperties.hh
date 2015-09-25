// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#ifndef DUMUX_BOX_PROPERTIES_HH
#define DUMUX_BOX_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>

#include <dumux/common/basicproperties.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxModel.
 */
namespace Dumux
{

namespace Properties
{
/*!
 * \addtogroup BoxModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the box-scheme
NEW_TYPE_TAG(BoxModel, INHERITS_FROM(ImplicitModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Grid);     //!< The type of the DUNE grid
NEW_PROP_TAG(GridView); //!< The type of the grid view

NEW_PROP_TAG(FVElementGeometry); //! The type of the finite-volume geometry in the box scheme

NEW_PROP_TAG(Problem);       //!< The type of the problem
NEW_PROP_TAG(Model);         //!< The type of the discretization
NEW_PROP_TAG(NumEq);         //!< Number of equations in the system of PDEs
NEW_PROP_TAG(LocalResidual); //!< The type of the local residual function
NEW_PROP_TAG(LocalJacobian); //!< The type of the local jacobian operator

NEW_PROP_TAG(JacobianAssembler); //!< Assembles the global jacobian matrix
NEW_PROP_TAG(JacobianMatrix); //!< Type of the global jacobian matrix
NEW_PROP_TAG(BoundaryTypes); //!< Stores the boundary types of a single degree of freedom
NEW_PROP_TAG(ElementBoundaryTypes); //!< Stores the boundary types on an element

NEW_PROP_TAG(PrimaryVariables); //!< A vector of primary variables within a sub-control volume
NEW_PROP_TAG(SolutionVector); //!< Vector containing all primary variable vector of the grid
NEW_PROP_TAG(ElementSolutionVector); //!< A vector of primary variables within a sub-control volume

NEW_PROP_TAG(VolumeVariables);  //!< The secondary variables within a sub-control volume
NEW_PROP_TAG(ElementVolumeVariables); //!< The secondary variables of all sub-control volumes in an element
NEW_PROP_TAG(FluxVariables); //!< Data required to calculate a flux over a face
NEW_PROP_TAG(BoundaryVariables); //!< Data required to calculate fluxes over boundary faces (outflow)

// high level simulation control
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time
NEW_PROP_TAG(NewtonMethod);     //!< The type of the newton method
NEW_PROP_TAG(NewtonController); //!< The type of the newton controller

#if HAVE_DUNE_PDELAB
// properties for the PDELab wrapper
NEW_PROP_TAG(LocalFEMSpace); //!< The local finite element space used for the finite element interpolation
NEW_PROP_TAG(ScalarGridFunctionSpace); //!< The used grid function space for a single finite element function
NEW_PROP_TAG(GridFunctionSpace); //!< The used grid function space
NEW_PROP_TAG(Constraints); //!< The constraints on the grid function space
NEW_PROP_TAG(ConstraintsTrafo); //!< The type of PDELab's constraints transformation
NEW_PROP_TAG(LocalOperator); //!< The type of the local operator used by PDELab
NEW_PROP_TAG(GridOperatorSpace); //!< The used grid operator space
#endif

//! Specify whether the jacobian matrix of the last iteration of a
//! time step should be re-used as the jacobian of the first iteration
//! of the next time step.
NEW_PROP_TAG(EnableJacobianRecycling);

//! Specify whether the jacobian matrix should be only reassembled for
//! elements where at least one vertex is above the specified
//! tolerance
NEW_PROP_TAG(EnablePartialReassemble);

//! Specify whether the time step should be increased in between
//! newton iterations to achive larger time step sizes
NEW_PROP_TAG(EnableTimeStepRampUp);

/*!
 * \brief Specify which kind of method should be used to numerically
 * calculate the partial derivatives of the residual.
 *
 * -1 means backward differences, 0 means central differences, 1 means
 * forward differences. By default we use central differences.
 */
NEW_PROP_TAG(NumericDifferenceMethod);

// mappers from local to global indices

//! maper for vertices
NEW_PROP_TAG(VertexMapper);
//! maper for elements
NEW_PROP_TAG(ElementMapper);
//! maper for degrees of freedom
NEW_PROP_TAG(DofMapper);
}
}

// \}

#endif
