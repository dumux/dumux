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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
#ifndef DUMUX_MULTIDOMAIN_NAVIERSTOKESTWOCNI_DARCYTWOCNI_PROPERTIES_HH
#define DUMUX_MULTIDOMAIN_NAVIERSTOKESTWOCNI_DARCYTWOCNI_PROPERTIES_HH

#include <dumux/porousmediumflow/2p2c/implicit/propertydefaults.hh>

//#include <appl/staggeredgrid/freeflow/navierstokes/navierstokes2cni/navierstokes2cnipropertydefaults.hh> // TODO ok?

//#include <appl/staggeredgrid/multidomain/common/multidomainpropertydefaults.hh> // TODO set in subproblems
//#include <appl/staggeredgrid/multidomain/common/subdomainpropertydefaults.hh> // TODO - not necessary?

namespace Dumux
{
namespace Properties
{
// Define how the Beavers-Joseph Condition is applied
NEW_PROP_TAG(CouplingBeaversJosephAsSolDependentDirichlet);

// Minimum accepted time step size to proceed simulation
NEW_PROP_TAG(NewtonAbortTimeStepSize);

// Factor to adapt the time step, if Newton failed
NEW_PROP_TAG(NewtonTimeStepReduction);

// Output frequency of Vtk files
NEW_PROP_TAG(VtkFreqOutput);

// Reference pressure used for pn_rel in Vtk files
NEW_PROP_TAG(VtkReferencePressure);

// Criteria for the interface Newton
NEW_PROP_TAG(CouplingSolverMaxSteps);
NEW_PROP_TAG(CouplingSolverResidualReduction);
NEW_PROP_TAG(CouplingSolverMaxRelativeShift);
NEW_PROP_TAG(CouplingSlopeLimitingFactor);

// Limit range of k and epsilon to positive values
NEW_PROP_TAG(NewtonKEpsilonLimiter);
} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_MULTIDOMAIN_NAVIERSTOKESTWOCNI_DARCYTWOCNI_PROPERTIES_HH
