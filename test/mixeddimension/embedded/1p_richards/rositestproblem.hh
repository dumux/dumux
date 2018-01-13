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
 * \brief A test problem the rootsystem coupled with richards in the bulk
 */
#ifndef DUMUX_ROSI_TEST_PROBLEM_HH
#define DUMUX_ROSI_TEST_PROBLEM_HH

#include "rootsystemtestproblem.hh"
#include "richardstestproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanager.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanagersimple.hh>
#include <dumux/mixeddimension/integrationpointsource.hh>

namespace Dumux
{
template <class TypeTag>
class RosiTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RosiTestProblem, INHERITS_FROM(MixedDimension));

// Set the problem property
SET_TYPE_PROP(RosiTestProblem, Problem, Dumux::RosiTestProblem<TypeTag>);

// Set the coupling manager
//SET_TYPE_PROP(RosiTestProblem, CouplingManager, Dumux::CCBBoxTreeEmbeddedCouplingManager<TypeTag>);
SET_TYPE_PROP(RosiTestProblem, CouplingManager, Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(RosiTestProblem, LowDimProblemTypeTag, TTAG(RootsystemTestCCProblem));
SET_TYPE_PROP(RosiTestProblem, BulkProblemTypeTag, TTAG(RichardsTestCCProblem));

// publish this problem in the sub problems
SET_TYPE_PROP(RootsystemTestCCProblem, GlobalProblemTypeTag, TTAG(RosiTestProblem));
SET_TYPE_PROP(RichardsTestCCProblem, GlobalProblemTypeTag, TTAG(RosiTestProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(RootsystemTestCCProblem, ParameterTree) : GET_PROP(TTAG(RosiTestProblem), ParameterTree) {};
SET_PROP(RichardsTestCCProblem, ParameterTree) : GET_PROP(TTAG(RosiTestProblem), ParameterTree) {};

// Set the point source type of the subproblems to an integration point source
SET_TYPE_PROP(RootsystemTestCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(RootsystemTestCCProblem), unsigned int>);
SET_TYPE_PROP(RootsystemTestCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(RootsystemTestCCProblem)>);
SET_TYPE_PROP(RichardsTestCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(RichardsTestCCProblem), unsigned int>);
SET_TYPE_PROP(RichardsTestCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(RichardsTestCCProblem)>);

SET_TYPE_PROP(RosiTestProblem, LinearSolver, ILU0BiCGSTABBackend);

}//end namespace properties

template <class TypeTag>
class RosiTestProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    // obtain types from the sub problem type tags
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

public:
    RosiTestProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimgridView)
    : ParentType(timeManager, bulkGridView, lowDimgridView)
    {}
};

} // end namespace Dumux

#endif
