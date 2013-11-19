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
 * \brief Sets default values for the MultiDomain properties
 */
#ifndef DUMUX_MULTIDOMAIN_PROPERTY_DEFAULTS_HH
#define DUMUX_MULTIDOMAIN_PROPERTY_DEFAULTS_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/timemanager.hh>

#include "splitandmerge.hh"
#include "multidomainmodel.hh"
#include "multidomainnewtoncontroller.hh"
#include "multidomainproperties.hh"

namespace Dumux
{
template <class TypeTag> class MultiDomain;
template <class TypeTag> class MultiDomainJacobianAssembler;
template <class TypeTag> class MultiDomainNewtonController;

namespace Properties
{

// Set property values for the coupled model
//SET_BOOL_PROP(MultiDomain, DoEnrichedCoupling, false);
SET_TYPE_PROP(MultiDomain, Model, MultiDomainModel<TypeTag>);

SET_TYPE_PROP(MultiDomain, JacobianAssembler, CoupledJacobianAssembler<TypeTag>);

SET_PROP(MultiDomain, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numEq> > type;
};

//! use the plain newton method for the coupled problems by default
SET_TYPE_PROP(MultiDomain, NewtonMethod, NewtonMethod<TypeTag>);

//! use the plain newton controller for coupled problems by default
SET_TYPE_PROP(MultiDomain, NewtonController, CoupledNewtonController<TypeTag>);

//! Set the default type of the time manager for coupled models
SET_TYPE_PROP(MultiDomain, TimeManager, TimeManager<TypeTag>);

SET_PROP(MultiDomain, JacobianMatrix)
{ private:
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numEq, numEq> > type;
};

SET_PROP(MultiDomain, NumEq)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) TypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) TypeTag2;

    enum {
        numEq1 = GET_PROP_VALUE(TypeTag1, NumEq),
        numEq2 = GET_PROP_VALUE(TypeTag2, NumEq)
    };
public:
    static const int value = numEq1; //TODO: why??
};

SET_PROP(MultiDomain, NumEq1)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) TypeTag1;
    enum {numEq = GET_PROP_VALUE(TypeTag1, NumEq)};
public:
    static const int value = numEq;
};

SET_PROP(MultiDomain, NumEq2)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) TypeTag2;
    enum {numEq = GET_PROP_VALUE(TypeTag2, NumEq)};
public:
    static const int value = numEq;
};

// Linear solver settings
SET_TYPE_PROP(MultiDomain, LinearSolver, BoxBiCGStabILU0Solver<TypeTag>};


SET_SCALAR_PROP(MultiDomain, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(MultiDomain, LinearSolverMaxIterations, 250);

SET_INT_PROP(MultiDomain, NewtonMaxTimeStepDivisions, 10);


// TODO: move to subdomainpropertydefaults.hh
// use the time manager for the coupled problem in the sub problems
SET_PROP(SubDomainProblem, TimeManager)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomain) CoupledTypeTag;
public:
    typedef typename GET_PROP_TYPE(CoupledTypeTag, TimeManager) type;
};

// use the time manager for the coupled problem in the sub problems
SET_PROP(SubDomainProblem, ParameterTree)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomain) MultiDomain;
    typedef typename GET_PROP(MultiDomain, ParameterTree) ParameterTree;
public:
    typedef typename ParameterTree::type type;

    static type &tree()
    { return ParameterTree::tree(); }

    static type &compileTimeParams()
    { return ParameterTree::compileTimeParams(); }

    static type &runTimeParams()
    { return ParameterTree::runTimeParams(); }

    static type &deprecatedRunTimeParams()
    { return ParameterTree::deprecatedRunTimeParams(); }

    static type &unusedNewRunTimeParams()
    { return ParameterTree::unusedNewRunTimeParams(); }
};

// \}
}
}
#endif

