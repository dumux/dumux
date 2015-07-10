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
 * \brief Sets default values for the MultiDomain properties
 */

#ifndef DUMUX_MULTIDOMAIN_PROPERTY_DEFAULTS_HH
#define DUMUX_MULTIDOMAIN_PROPERTY_DEFAULTS_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/multidomain/subdomainset.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>

#include "subdomainpropertydefaults.hh"
#include "multidomainmodel.hh"
#include "multidomainproperties.hh"
#include "multidomainnewtoncontroller.hh"
#include "splitandmerge.hh"

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/timemanager.hh>

namespace Dumux
{
template <class TypeTag> class MultiDomainModel;
template <class TypeTag> class MultiDomainJacobianAssembler;
template <class TypeTag> class MultiDomainNewtonController;

namespace Properties
{

SET_PROP(MultiDomain, MultiDomainGrid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) HostGrid;
    typedef typename Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,4> MDGridTraits;
public:
    typedef typename Dune::MultiDomainGrid<HostGrid, MDGridTraits> type;
};

SET_PROP(MultiDomain, MultiDomainGridFunctionSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename Dune::PDELab::LexicographicOrderingTag OrderingTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;
    typedef typename GET_PROP_TYPE(SubTypeTag1, GridFunctionSpace) GridFunctionSpace1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, GridFunctionSpace) GridFunctionSpace2;
public:
    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<MDGrid,
                                                                    Dune::PDELab::ISTLVectorBackend<>,
                                                                    OrderingTag,
                                                                    GridFunctionSpace1,
                                                                    GridFunctionSpace2> type;
};

// set the subdomain equality condition by default
SET_PROP(MultiDomain, MultiDomainCondition)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
public:
    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<MDGrid> type;
};

SET_PROP(MultiDomain, MultiDomainSubProblem1)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridFunctionSpace) MDGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(SubTypeTag1, Constraints) Constraints1;
    typedef typename GET_PROP_TYPE(SubTypeTag1, LocalOperator) LocalOperator1;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCondition) MDCondition;
    typedef typename GET_PROP_TYPE(SubTypeTag1, GridFunctionSpace) GridFunctionSpace1;
public:
    typedef Dune::PDELab::MultiDomain::SubProblem<MDGridFunctionSpace,
                                                  MDGridFunctionSpace,
                                                  LocalOperator1, MDCondition,
                                                  0> type;
};

SET_PROP(MultiDomain, MultiDomainSubProblem2)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridFunctionSpace) MDGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Constraints) Constraints2;
    typedef typename GET_PROP_TYPE(SubTypeTag2, LocalOperator) LocalOperator2;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCondition) MDCondition;
    typedef typename GET_PROP_TYPE(SubTypeTag2, GridFunctionSpace) GridFunctionSpace2;
public:
    typedef Dune::PDELab::MultiDomain::SubProblem<MDGridFunctionSpace,
                                                  MDGridFunctionSpace,
                                                  LocalOperator2, MDCondition,
                                                  1> type;
};

SET_PROP(MultiDomain, MultiDomainCoupling)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem1) MDSubProblem1;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem2) MDSubProblem2;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCouplingLocalOperator) MDCouplingLocalOperator;
public:
    typedef Dune::PDELab::MultiDomain::Coupling<MDSubProblem1,
                                                MDSubProblem2,
                                                MDCouplingLocalOperator> type;
};

// set trivial constraints transformation by default
// TODO create a proper constraint transformation
SET_PROP(MultiDomain, MultiDomainConstraintsTrafo)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridFunctionSpace) MDGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef typename MDGridFunctionSpace::template ConstraintsContainer<Scalar>::Type type;
};

SET_PROP(MultiDomain, MultiDomainGridOperator)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridFunctionSpace) MDGridFunctionSpace;
    typedef Dune::PDELab::ISTLMatrixBackend MBE;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem1) MDSubProblem1;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem2) MDSubProblem2;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCoupling) MDCoupling;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainConstraintsTrafo) MDConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dune::PDELab::MultiDomain::GridOperator<
                    MDGridFunctionSpace, MDGridFunctionSpace,
                    MBE, Scalar, Scalar, Scalar,
                    MDConstraintsTrafo, MDConstraintsTrafo,
                    MDSubProblem1, MDSubProblem2, MDCoupling> type;
};

SET_PROP(MultiDomain, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridOperator) MDGridOperator;
public:
    typedef typename MDGridOperator::Traits::Jacobian type;
};

SET_INT_PROP(MultiDomain, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));


// Set property values for the coupled model
//SET_BOOL_PROP(MultiDomain, DoEnrichedCoupling, false);
SET_TYPE_PROP(MultiDomain, Model, MultiDomainModel<TypeTag>);

SET_PROP(MultiDomain, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numEq> > type;
};

// use the plain newton method for the coupled problems by default
SET_TYPE_PROP(MultiDomain, NewtonMethod, NewtonMethod<TypeTag>);

// use the plain newton controller for coupled problems by default
SET_TYPE_PROP(MultiDomain, NewtonController, MultiDomainNewtonController<TypeTag>);

// Set the default type of the time manager for coupled models
SET_TYPE_PROP(MultiDomain, TimeManager, TimeManager<TypeTag>);

// needed to define size of ImplicitBase's PrimaryVariables
SET_PROP(MultiDomain, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) TypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TypeTag2;

    enum {
        numEq1 = GET_PROP_VALUE(TypeTag1, NumEq),
        numEq2 = GET_PROP_VALUE(TypeTag2, NumEq)
    };
public:
    static const int value = numEq1;
};

SET_PROP(MultiDomain, NumEq1)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) TypeTag1;
    enum {numEq = GET_PROP_VALUE(TypeTag1, NumEq)};
public:
    static const int value = numEq;
};

SET_PROP(MultiDomain, NumEq2)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TypeTag2;
    enum {numEq = GET_PROP_VALUE(TypeTag2, NumEq)};
public:
    static const int value = numEq;
};

// set the type of the linear solver
SET_TYPE_PROP(MultiDomain, LinearSolver, BoxBiCGStabILU0Solver<TypeTag>);

// set the minimum residual reduction of the linear solver
SET_SCALAR_PROP(MultiDomain, LinearSolverResidualReduction, 1e-6);

// set the default number of maximum iterations for the linear solver
SET_INT_PROP(MultiDomain, LinearSolverMaxIterations, 250);

// set the maximum time step divisions
SET_INT_PROP(MultiDomain, NewtonMaxTimeStepDivisions, 10);

// set the routines for splitting and merging solution vectors
SET_TYPE_PROP(MultiDomain, SplitAndMerge, SplitAndMerge<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_PROPERTY_DEFAULTS_HH

