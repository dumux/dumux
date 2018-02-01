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
 *
 * \brief A test problem for the coupled problem (1p)
 */
#ifndef DUMUX_STOKES_DARCY_1P_HORIZONTAL_TEST_PROBLEM_HH
#define DUMUX_STOKES_DARCY_1P_HORIZONTAL_TEST_PROBLEM_HH

#include "darcytestproblem.hh"
#include "stokestestproblem.hh"

#include <dumux/multidomain/problem.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/map.hh>

#include <dumux/multidomain/staggered-ccfv/properties.hh>
#include <dumux/multidomain/staggered-ccfv/model.hh>
#include <dumux/multidomain/staggered-ccfv/assembler.hh>
#include <dumux/multidomain/staggered-ccfv/newtoncontroller.hh>
#include <dumux/multidomain/staggered-ccfv/subproblemlocaljacobian.hh>


namespace Dumux
{
template <class TypeTag>
class HorizontalFlowCoupledStokesDarcyProblem;

namespace Properties
{
NEW_TYPE_TAG(HorizontalFlowCoupledStokesDarcyProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, Problem, Dumux::HorizontalFlowCoupledStokesDarcyProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, DarcyProblemTypeTag, TTAG(DarcyTestProblem));
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, StokesProblemTypeTag, TTAG(StokesTestProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(DarcyTestProblem, GlobalProblemTypeTag, TTAG(HorizontalFlowCoupledStokesDarcyProblem));
SET_TYPE_PROP(StokesTestProblem, GlobalProblemTypeTag, TTAG(HorizontalFlowCoupledStokesDarcyProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(DarcyTestProblem, ParameterTree) : GET_PROP(TTAG(HorizontalFlowCoupledStokesDarcyProblem), ParameterTree) {};
SET_PROP(StokesTestProblem, ParameterTree) : GET_PROP(TTAG(HorizontalFlowCoupledStokesDarcyProblem), ParameterTree) {};

NEW_PROP_TAG(DarcyToStokesMapValue);
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);

#if HAVE_UMFPACK
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif

NEW_PROP_TAG(StokesData);
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, StokesData, StokesData<TypeTag>);

NEW_PROP_TAG(DarcyData);
SET_TYPE_PROP(HorizontalFlowCoupledStokesDarcyProblem, DarcyData, DarcyData<TypeTag>);

}//end namespace properties

template <class TypeTag>
class HorizontalFlowCoupledStokesDarcyProblem : public MultiDomainProblem<TypeTag>
{
    using ParentType = MultiDomainProblem<TypeTag>;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    // obtain types from the sub problem type tags
    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

public:
    HorizontalFlowCoupledStokesDarcyProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
        verbose_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Verbose);
    }


    bool shouldWriteOutput() const //define output
    {
        return true;
    }

private:
    bool verbose_;
};

} //end namespace

#endif
