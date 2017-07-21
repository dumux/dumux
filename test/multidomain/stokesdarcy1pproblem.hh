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
#include <dumux/multidomain/boundingboxtree/couplingdata.hh>
#include <dumux/multidomain/boundingboxtree/map.hh>
#include <dumux/multidomain/staggeredgrid/properties.hh>

// #include<dumux/porenetworkflow/common/functions.hh>

#include <dumux/multidomain/staggeredgrid/model.hh>
#include <dumux/multidomain/staggeredgrid/assembler.hh>
#include <dumux/multidomain/staggeredgrid/newtoncontroller.hh>
#include <dumux/multidomain/staggeredgrid/subproblemlocaljacobian.hh>


namespace Dumux
{
template <class TypeTag>
class TestCoupledStokesDarcyBBTProblem; // boundingboxtree

namespace Properties
{
NEW_TYPE_TAG(TestCoupledStokesDarcyBBTProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, Problem, Dumux::TestCoupledStokesDarcyBBTProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, DarcyProblemTypeTag, TTAG(DarcyTestProblem));
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, StokesProblemTypeTag, TTAG(StokesTestProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(DarcyTestProblem, GlobalProblemTypeTag, TTAG(TestCoupledStokesDarcyBBTProblem));
SET_TYPE_PROP(StokesTestProblem, GlobalProblemTypeTag, TTAG(TestCoupledStokesDarcyBBTProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(DarcyTestProblem, ParameterTree) : GET_PROP(TTAG(TestCoupledStokesDarcyBBTProblem), ParameterTree) {};
SET_PROP(StokesTestProblem, ParameterTree) : GET_PROP(TTAG(TestCoupledStokesDarcyBBTProblem), ParameterTree) {};

// SET_BOOL_PROP(TestCoupledStokesDarcyBBTProblem, MultiDimensionUseIterativeSolver, true);

NEW_PROP_TAG(DarcyToStokesMapValue); // TODO: make specialized map value class
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);

NEW_PROP_TAG(StokesData);
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, StokesData, StokesData<TypeTag>);

NEW_PROP_TAG(DarcyData);
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, DarcyData, DarcyData<TypeTag>);

//! Set the BaseModel to MultiDomainModel
// SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, Model, MultiDomainModelForStaggered<TypeTag>);
// SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, JacobianAssembler,MultiDomainAssemblerForStaggered<TypeTag>);
// SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, NewtonController, MultiDomainNewtonControllerForStaggered<TypeTag>);

#if HAVE_UMFPACK
SET_TYPE_PROP(TestCoupledStokesDarcyBBTProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif

}//end namespace properties

template <class TypeTag>
class TestCoupledStokesDarcyBBTProblem : public MultiDomainProblem<TypeTag>
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

    using StokesPrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, PrimaryVariables);
    using DarcyPrimaryVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, PrimaryVariables);

public:
    TestCoupledStokesDarcyBBTProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
        verbose_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Verbose);
    }


    bool shouldWriteOutput() const //define output
    {
        return true;
//         return (this->timeManager().willBeFinished());
    }

    void postTimeStep()
    {
        // calculate the mass flux entering the stokes domain
//         StokesPrimaryVariables inflow;
//         const auto left = this->couplingManager().stokesProblem().bBoxMin()[0] + 1e-3;
//         this->couplingManager().stokesProblem().calculateFluxAcrossLayer(inflow, 0, left);
//
        //calculate the mass flux leaving  the darcy domain
//         const DarcyPrimaryVariables outflow = Functions<DarcyProblemTypeTag>::boundaryFlux(this->couplingManager().darcyProblem(), "max", 1);
//
//         if(verbose_)
//         {
//             std::cout << "IN Darcy: " << inflow << " kg/s" << std::endl;
//             std::cout << "OUT PNM : " << outflow <<  " kg/s"<<  std::endl;
//         }
//         // at the end of the simulation, check if everything is mass conservative
//         if(this->timeManager().willBeFinished())
//         {
//             if(Dune::FloatCmp::ne<Scalar>(inflow, outflow, 1e-6))
//             {
//                 std::cout << "Test failed: Inflow differs from outflow: " << inflow - outflow << std::endl;
//                 std::exit(EXIT_FAILURE);
//             }
//         }
    }
private:
    bool verbose_;
};

} //end namespace

#endif
