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
 * \file
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup MultidomainModel
 * \brief Specify properties required for the coupled model
 */
#ifndef DUMUX_MULTIDOMAIN_PROPERTIES_HH
#define DUMUX_MULTIDOMAIN_PROPERTIES_HH

#include <dumux/implicit/common/implicitproperties.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/linear/linearsolverproperties.hh>

namespace Dumux
{

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for problems which utilize the coupling approach
NEW_TYPE_TAG(MultiDomain, INHERITS_FROM(NewtonMethod, LinearSolverTypeTag, ImplicitBase));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Specifies the model
NEW_PROP_TAG(Model);

//! Specifies the type tag of the first sub-problem
NEW_PROP_TAG(SubDomain1TypeTag);

//! Specifies the type tag of the second sub-problem
NEW_PROP_TAG(SubDomain2TypeTag);

//! Specifies the type tag of the other sub-problem
NEW_PROP_TAG(OtherSubDomainTypeTag);

//! Specifies the type tag of coupled problem
NEW_PROP_TAG(MultiDomainTypeTag);

//! Specifies the host grid
NEW_PROP_TAG(Grid);

//! Specifies the multidomain grid
NEW_PROP_TAG(MultiDomainGrid);

//! Specifies the number of equations in submodel 1
NEW_PROP_TAG(NumEq1);

//! Specifies the number of equations in submodel 2
NEW_PROP_TAG(NumEq2);

//! Specifies the fluidsystem that is used in the subdomains
NEW_PROP_TAG(FluidSystem);

//! the maximum allowed number of timestep divisions for the
//! Newton solver
NEW_PROP_TAG(NewtonMaxTimeStepDivisions);

//! Specifies the multidomain grid function space
NEW_PROP_TAG(MultiDomainGridFunctionSpace);

//! Specifies the multidomain grid operator
NEW_PROP_TAG(MultiDomainGridOperator);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(NewtonWriteConvergence);

//! Specifies the equality conditions
NEW_PROP_TAG(MultiDomainCondition);

//! Specifies the multidomain type based subproblem for subdomain 1
NEW_PROP_TAG(MultiDomainSubProblem1);

//! Specifies the multidomain type based subproblem for subdomain 2
NEW_PROP_TAG(MultiDomainSubProblem2);

//! the local coupling operator for use with dune-multidomain
NEW_PROP_TAG(MultiDomainCouplingLocalOperator);

//! Specifies the multidomain coupling
NEW_PROP_TAG(MultiDomainCoupling);

//! Property tag for the multidomain constraints transformation
NEW_PROP_TAG(MultiDomainConstraintsTrafo);
NEW_PROP_TAG(ConstraintsTrafo);

//! the routines that are used to split and merge solution vectors
NEW_PROP_TAG(SplitAndMerge);

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_PROPERTIES_HH
