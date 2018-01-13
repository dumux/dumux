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
 * \brief A test problem for the 1d2d coupled problem:
 *        1d fractures living on the element facets of a 2d matrix
 */
#ifndef DUMUX_1D2D_FACET_TEST_PROBLEM_HH
#define DUMUX_1D2D_FACET_TEST_PROBLEM_HH

#include "fractureproblem.hh"
#include "matrixproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/facet/gmshdualfacetgridcreator.hh>
#include <dumux/mixeddimension/facet/mpfa/couplingmanager.hh>

namespace Dumux
{
template <class TypeTag>
class OnePTwoCFacetCouplingProblem;

namespace Properties
{
// Type tag of the isothermal and non-isotherman global Problem
NEW_TYPE_TAG(OnePTwoCFacetCoupling, INHERITS_FROM(MixedDimension));
NEW_TYPE_TAG(OnePTwoCIFacetCoupling, INHERITS_FROM(OnePTwoCFacetCoupling));
NEW_TYPE_TAG(OnePTwoCNIFacetCoupling, INHERITS_FROM(OnePTwoCFacetCoupling));

// Set the problem property
SET_TYPE_PROP(OnePTwoCFacetCoupling, Problem, Dumux::OnePTwoCFacetCouplingProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(OnePTwoCFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(OnePTwoCIFacetCoupling, BulkProblemTypeTag, TTAG(OnePTwoCICCMpfaMatrixProblem));
SET_TYPE_PROP(OnePTwoCIFacetCoupling, LowDimProblemTypeTag, TTAG(OnePTwoCICCFractureProblem));
SET_TYPE_PROP(OnePTwoCNIFacetCoupling, BulkProblemTypeTag, TTAG(OnePTwoCNICCMpfaMatrixProblem));
SET_TYPE_PROP(OnePTwoCNIFacetCoupling, LowDimProblemTypeTag, TTAG(OnePTwoCNICCFractureProblem));

// The coupling manager
SET_TYPE_PROP(OnePTwoCFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);

// The linear solver to be used
SET_TYPE_PROP(OnePTwoCFacetCoupling, LinearSolver, ILU0BiCGSTABBackend);

// The sub-problems need to know the global problem's type tag
SET_TYPE_PROP(OnePTwoCICCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(OnePTwoCIFacetCoupling));
SET_TYPE_PROP(OnePTwoCICCFractureProblem, GlobalProblemTypeTag, TTAG(OnePTwoCIFacetCoupling));
SET_TYPE_PROP(OnePTwoCNICCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(OnePTwoCNIFacetCoupling));
SET_TYPE_PROP(OnePTwoCNICCFractureProblem, GlobalProblemTypeTag, TTAG(OnePTwoCNIFacetCoupling));

// The subproblems inherit the parameter tree from this problem
SET_PROP(OnePTwoCICCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(OnePTwoCIFacetCoupling), ParameterTree) {};
SET_PROP(OnePTwoCICCFractureProblem, ParameterTree) : GET_PROP(TTAG(OnePTwoCIFacetCoupling), ParameterTree) {};
SET_PROP(OnePTwoCNICCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(OnePTwoCNIFacetCoupling), ParameterTree) {};
SET_PROP(OnePTwoCNICCFractureProblem, ParameterTree) : GET_PROP(TTAG(OnePTwoCNIFacetCoupling), ParameterTree) {};

// Set the grids for the two sub-problems
SET_TYPE_PROP(OnePTwoCICCMpfaMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(OnePTwoCICCFractureProblem, Grid, Dune::FoamGrid<1,2>);
SET_TYPE_PROP(OnePTwoCNICCMpfaMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(OnePTwoCNICCFractureProblem, Grid, Dune::FoamGrid<1,2>);

}

template <class TypeTag>
class OnePTwoCFacetCouplingProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

public:

    OnePTwoCFacetCouplingProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimGridView)
    : ParentType(timeManager, bulkGridView, lowDimGridView)
    {}
};

} //end namespace

#endif
