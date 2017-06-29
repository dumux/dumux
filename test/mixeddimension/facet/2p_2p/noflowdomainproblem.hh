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
#ifndef DUMUX_NOFLOWDOMAIN_TEST_PROBLEM_HH
#define DUMUX_NOFLOWDOMAIN_TEST_PROBLEM_HH

#include "noflowfractureproblem.hh"
#include "noflowmatrixproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/facet/properties.hh>
#include <dumux/mixeddimension/facet/gmshdualfacetgridcreator.hh>
#include <dumux/mixeddimension/facet/mpfa/couplingmanager.hh>

namespace Dumux
{
template <class TypeTag>
class NoFlowDomainFacetCouplingProblem;

namespace Properties
{
// Set the type tag and properties
NEW_TYPE_TAG(NoFlowDomain2pFacetCoupling, INHERITS_FROM(MixedDimensionFacetCoupling));

// the problem property
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, Problem, Dumux::NoFlowDomainFacetCouplingProblem<TypeTag>);

// set the facet coupling grid creator
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, BulkProblemTypeTag, TTAG(NoFlowDomainMatrixProblem));
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, LowDimProblemTypeTag, TTAG(NoFlowDomainFractureProblem));

// The coupling manager
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);

// The linear solver to be used
SET_TYPE_PROP(NoFlowDomain2pFacetCoupling, LinearSolver, ILU0BiCGSTABBackend<TypeTag>);

// The sub-problems need to know the global problem's type tag
SET_TYPE_PROP(NoFlowDomainMatrixProblem, GlobalProblemTypeTag, TTAG(NoFlowDomain2pFacetCoupling));
SET_TYPE_PROP(NoFlowDomainFractureProblem, GlobalProblemTypeTag, TTAG(NoFlowDomain2pFacetCoupling));

// The subproblems inherit the parameter tree from this problem
SET_PROP(NoFlowDomainMatrixProblem, ParameterTree) : GET_PROP(TTAG(NoFlowDomain2pFacetCoupling), ParameterTree) {};
SET_PROP(NoFlowDomainFractureProblem, ParameterTree) : GET_PROP(TTAG(NoFlowDomain2pFacetCoupling), ParameterTree) {};

// Set the grids for the two sub-problems
SET_TYPE_PROP(NoFlowDomainMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
SET_TYPE_PROP(NoFlowDomainFractureProblem, Grid, Dune::FoamGrid<1,2>);
}

template <class TypeTag>
class NoFlowDomainFacetCouplingProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

public:

    NoFlowDomainFacetCouplingProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimGridView)
    : ParentType(timeManager, bulkGridView, lowDimGridView)
    {}

    bool shouldWriteOutput() const
    { return true; }
};

} //end namespace

#endif
