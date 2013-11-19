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
#ifndef DUMUX_COUPLED_PROPERTIES_HH
#define DUMUX_COUPLED_PROPERTIES_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dumux/nonlinear/newtonmethod.hh>

#include <dumux/implicit/common/implicitproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

#include <dumux/common/timemanager.hh>

/*!
 * \file
 * \brief Specify properties required for the coupled model
 */
namespace Dumux
{
namespace Properties
{
/*!
 * \addtogroup ModelCoupling
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for problems which utilize the coupling approach
NEW_TYPE_TAG(MultiDomain, INHERITS_FROM(LinearSolverTypeTag, ImplicitBase));

//! The type tag from which sub-problems of coupling models inherit
NEW_TYPE_TAG(CoupledSubProblem);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
// FROM COUPLEDPROPERTIES.HH
//////////////////////////////////////////////////////////////////

//! Specifies the type tag of the first sub-problem
NEW_PROP_TAG(SubProblem1TypeTag);

//! Specifies the type tag of the second sub-problem
NEW_PROP_TAG(SubProblem2TypeTag);

//! Specifies the type tag of the other sub-problem
NEW_PROP_TAG(OtherSubProblemTypeTag);

//! Specifies the type tag of coupled problem
NEW_PROP_TAG(CoupledProblemTypeTag);

//! Specifies the local jacobian of a meta element
NEW_PROP_TAG(CoupledLocalJacobian);

//! Specifies the jacobian assembler
NEW_PROP_TAG(JacobianAssembler);

//! Specifies the type of the jacobian matrix as used for the linear
//! solver
NEW_PROP_TAG(JacobianMatrix);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every newton iteration
NEW_PROP_TAG(NewtonWriteConvergence);

//! the maximum allowed number of timestep divisions for the
//! Newton solver
NEW_PROP_TAG(NewtonMaxTimeStepDivisions);

//! Specifies the model
NEW_PROP_TAG(Model);

//! Specifies the time manager
NEW_PROP_TAG(TimeManager);

//! Specifies the number of equations in the coupled model
NEW_PROP_TAG(NumEq);

//! Specifies the number of equations in submodel 1
NEW_PROP_TAG(NumEq1);

//! Specifies the number of equations in submodel 2
NEW_PROP_TAG(NumEq2);

//! Specifies the fluidsystem that is used in the subdomains
NEW_PROP_TAG(FluidSystem);

//! Specifies whether the enriched(mortar) coupling is used (set to false by default)
NEW_PROP_TAG(DoEnrichedCoupling);


//////////////////////////////////////////////////////////////////
// FROM MULTIDOMAINPROPERTIES.HH
//////////////////////////////////////////////////////////////////

//! Specifies the host grid
NEW_PROP_TAG(Grid);

//! Specifies the multidomain grid
NEW_PROP_TAG(MDGrid);

//! Specifies the multidomain grid function space
NEW_PROP_TAG(MDGridFunctionSpace);

//! Specifies the equality conditions
NEW_PROP_TAG(MDCondition);

//! Specifies the multidomain type based subproblem for subdomain 1
NEW_PROP_TAG(MDSubProblem1);

//! Specifies the multidomain type based subproblem for subdomain 2
NEW_PROP_TAG(MDSubProblem2);

//! the local coupling operator for use with dune-multidomain
NEW_PROP_TAG(MDCouplingLocalOperator);

//! Specifies the multidomain coupling
NEW_PROP_TAG(MDCoupling);

//! Property tag for the multidomain constraints transformation
NEW_PROP_TAG(MDConstraintsTrafo);

//! Specifies the multidomain grid operator
NEW_PROP_TAG(MDGridOperator);

NEW_PROP_TAG(ConstraintsTrafo);

}
}
#endif
