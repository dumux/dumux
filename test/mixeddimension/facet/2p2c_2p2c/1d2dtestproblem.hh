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
class TwoPTwoCFacetCouplingProblem;

namespace Properties
{
// Type tag of the isothermal and non-isotherman global Problem
NEW_TYPE_TAG(TwoPTwoCFacetCoupling, INHERITS_FROM(MixedDimension));
NEW_TYPE_TAG(TwoPTwoCIFacetCoupling, INHERITS_FROM(TwoPTwoCFacetCoupling));
NEW_TYPE_TAG(TwoPTwoCNIFacetCoupling, INHERITS_FROM(TwoPTwoCFacetCoupling));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCFacetCoupling, Problem, Dumux::TwoPTwoCFacetCouplingProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(TwoPTwoCFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoPTwoCIFacetCoupling, BulkProblemTypeTag, TTAG(TwoPTwoCICCMpfaMatrixProblem));
SET_TYPE_PROP(TwoPTwoCIFacetCoupling, LowDimProblemTypeTag, TTAG(TwoPTwoCICCFractureProblem));
SET_TYPE_PROP(TwoPTwoCNIFacetCoupling, BulkProblemTypeTag, TTAG(TwoPTwoCNICCMpfaMatrixProblem));
SET_TYPE_PROP(TwoPTwoCNIFacetCoupling, LowDimProblemTypeTag, TTAG(TwoPTwoCNICCFractureProblem));

// The coupling manager
SET_TYPE_PROP(TwoPTwoCFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);

// The linear solver to be used
SET_TYPE_PROP(TwoPTwoCFacetCoupling, LinearSolver, SuperLUBackend<TypeTag>);

// The sub-problems need to know the global problem's type tag
SET_TYPE_PROP(TwoPTwoCICCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(TwoPTwoCIFacetCoupling));
SET_TYPE_PROP(TwoPTwoCICCFractureProblem, GlobalProblemTypeTag, TTAG(TwoPTwoCIFacetCoupling));
SET_TYPE_PROP(TwoPTwoCNICCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(TwoPTwoCNIFacetCoupling));
SET_TYPE_PROP(TwoPTwoCNICCFractureProblem, GlobalProblemTypeTag, TTAG(TwoPTwoCNIFacetCoupling));

// The subproblems inherit the parameter tree from this problem
SET_PROP(TwoPTwoCICCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(TwoPTwoCIFacetCoupling), ParameterTree) {};
SET_PROP(TwoPTwoCICCFractureProblem, ParameterTree) : GET_PROP(TTAG(TwoPTwoCIFacetCoupling), ParameterTree) {};
SET_PROP(TwoPTwoCNICCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(TwoPTwoCNIFacetCoupling), ParameterTree) {};
SET_PROP(TwoPTwoCNICCFractureProblem, ParameterTree) : GET_PROP(TTAG(TwoPTwoCNIFacetCoupling), ParameterTree) {};

// Set the grids for the two sub-problems
SET_TYPE_PROP(TwoPTwoCICCMpfaMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(TwoPTwoCICCFractureProblem, Grid, Dune::FoamGrid<1,2>);
SET_TYPE_PROP(TwoPTwoCNICCMpfaMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(TwoPTwoCNICCFractureProblem, Grid, Dune::FoamGrid<1,2>);

}

template <class TypeTag>
class TwoPTwoCFacetCouplingProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

public:

    TwoPTwoCFacetCouplingProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimGridView)
    : ParentType(timeManager, bulkGridView, lowDimGridView)
    {}

    bool shouldWriteOutput() const
    { return true; }
};

} //end namespace

#endif
