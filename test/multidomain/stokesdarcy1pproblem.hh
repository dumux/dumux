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
 * \brief A test problem for the coupled problem
 */
#ifndef DUMUX_STOKES_DARCY_1P_TEST_PROBLEM_HH
#define DUMUX_STOKES_DARCY_1P_TEST_PROBLEM_HH

#include "darcytestproblem.hh"
#include "stokestestproblem.hh"

#include <dumux/multidomain/boundingboxtree/problem.hh>
#include <dumux/multidomain/boundingboxtree/couplingmanager.hh>
#include <dumux/multidomain/boundingboxtree/map.hh>

#include <dumux/multidomain/staggeredgrid/properties.hh>
#include <dumux/multidomain/staggeredgrid/model.hh>
#include <dumux/multidomain/staggeredgrid/assembler.hh>
#include <dumux/multidomain/staggeredgrid/newtoncontroller.hh>
#include <dumux/multidomain/staggeredgrid/subproblemlocaljacobian.hh>


namespace Dumux
{
template <class TypeTag>
class TestCoupledStokesDarcyProblem;

namespace Properties
{
NEW_TYPE_TAG(TestCoupledStokesDarcyProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, Problem, Dumux::TestCoupledStokesDarcyProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, DarcyProblemTypeTag, TTAG(DarcyTestProblem));
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, StokesProblemTypeTag, TTAG(StokesTestProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(DarcyTestProblem, GlobalProblemTypeTag, TTAG(TestCoupledStokesDarcyProblem));
SET_TYPE_PROP(StokesTestProblem, GlobalProblemTypeTag, TTAG(TestCoupledStokesDarcyProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(DarcyTestProblem, ParameterTree) : GET_PROP(TTAG(TestCoupledStokesDarcyProblem), ParameterTree) {};
SET_PROP(StokesTestProblem, ParameterTree) : GET_PROP(TTAG(TestCoupledStokesDarcyProblem), ParameterTree) {};

// SET_BOOL_PROP(TestCoupledStokesDarcyProblem, MultiDimensionUseIterativeSolver, true);

NEW_PROP_TAG(DarcyToStokesMapValue); // TODO: make specialized map value class
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);

#if HAVE_UMFPACK
SET_TYPE_PROP(TestCoupledStokesDarcyProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif

}//end namespace properties

template <class TypeTag>
class TestCoupledStokesDarcyProblem : public MultiDomainProblem<TypeTag>
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
    TestCoupledStokesDarcyProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
        verbose_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Verbose);
    }


    bool shouldWriteOutput() const //define output
    {
        return true;
    }

    void postTimeStep()
    {

    }
private:
    bool verbose_;
};

} //end namespace

#endif
