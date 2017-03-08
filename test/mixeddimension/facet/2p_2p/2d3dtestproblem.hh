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
 * \brief A test problem for the 2d3d coupled problem:
 *        2d fractures living on the element facets of a 3d matrix
 */
#ifndef DUMUX_2D3D_2P_FACET_TEST_PROBLEM_HH
#define DUMUX_2D3D_2P_FACET_TEST_PROBLEM_HH

#include "fractureproblem.hh"
#include "matrixproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/facet/gmshdualfacetgridcreator.hh>
#include <dumux/mixeddimension/facet/mpfa/couplingmanager.hh>

namespace Dumux
{
template <class TypeTag>
class TwoPFacetCouplingProblem;

namespace Properties
{
// Type tag of the isothermal and non-isotherman global Problem
NEW_TYPE_TAG(TwoPFacetCoupling, INHERITS_FROM(MixedDimension));
NEW_TYPE_TAG(TwoPNIFacetCoupling, INHERITS_FROM(MixedDimension));

// Set the problem property
SET_TYPE_PROP(TwoPFacetCoupling, Problem, Dumux::TwoPFacetCouplingProblem<TypeTag>);
SET_TYPE_PROP(TwoPNIFacetCoupling, Problem, Dumux::TwoPFacetCouplingProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(TwoPFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);
SET_TYPE_PROP(TwoPNIFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoPFacetCoupling, BulkProblemTypeTag, TTAG(TwoPCCMpfaMatrixProblem));
SET_TYPE_PROP(TwoPFacetCoupling, LowDimProblemTypeTag, TTAG(TwoPCCMpfaFractureProblem));
SET_TYPE_PROP(TwoPNIFacetCoupling, BulkProblemTypeTag, TTAG(TwoPNICCMpfaMatrixProblem));
SET_TYPE_PROP(TwoPNIFacetCoupling, LowDimProblemTypeTag, TTAG(TwoPNICCMpfaFractureProblem));

// The coupling manager
SET_TYPE_PROP(TwoPFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);
SET_TYPE_PROP(TwoPNIFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);

// The linear solver to be used
SET_TYPE_PROP(TwoPFacetCoupling, LinearSolver, ILU0BiCGSTABBackend<TypeTag>);
SET_TYPE_PROP(TwoPNIFacetCoupling, LinearSolver, ILU0BiCGSTABBackend<TypeTag>);

// The sub-problems need to know the global problem's type tag
SET_TYPE_PROP(TwoPCCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(TwoPFacetCoupling));
SET_TYPE_PROP(TwoPCCMpfaFractureProblem, GlobalProblemTypeTag, TTAG(TwoPFacetCoupling));
SET_TYPE_PROP(TwoPNICCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(TwoPNIFacetCoupling));
SET_TYPE_PROP(TwoPNICCMpfaFractureProblem, GlobalProblemTypeTag, TTAG(TwoPNIFacetCoupling));

// The subproblems inherit the parameter tree from this problem
SET_PROP(TwoPCCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(TwoPFacetCoupling), ParameterTree) {};
SET_PROP(TwoPCCMpfaFractureProblem, ParameterTree) : GET_PROP(TTAG(TwoPFacetCoupling), ParameterTree) {};
SET_PROP(TwoPNICCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(TwoPNIFacetCoupling), ParameterTree) {};
SET_PROP(TwoPNICCMpfaFractureProblem, ParameterTree) : GET_PROP(TTAG(TwoPNIFacetCoupling), ParameterTree) {};

// Set the grids for the two sub-problems
SET_TYPE_PROP(TwoPCCMpfaMatrixProblem, Grid, Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(TwoPCCMpfaFractureProblem, Grid, Dune::FoamGrid<2, 3>);
SET_TYPE_PROP(TwoPNICCMpfaMatrixProblem, Grid, Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>);
SET_TYPE_PROP(TwoPNICCMpfaFractureProblem, Grid, Dune::FoamGrid<2, 3>);

}

template <class TypeTag>
class TwoPFacetCouplingProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

public:

    TwoPFacetCouplingProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimGridView)
    : ParentType(timeManager, bulkGridView, lowDimGridView)
    {}

    bool shouldWriteOutput() const
    { return true; }
};

} //end namespace

#endif
