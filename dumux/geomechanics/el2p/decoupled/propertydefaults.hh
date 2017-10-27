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
 * \brief Defines the properties required for the two phase linear-elastic model.
 *
 * This class inherits from the properties of the two-phase model and
 * from the properties of the simple linear-elastic model
 */

#ifndef DUMUX_DECOUPLED_ELASTIC_PROPERTY_DEFAULTS_HH
#define DUMUX_DECOUPLED_ELASTIC_PROPERTY_DEFAULTS_HH

#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>


#include "model.hh"
#include "fluxvariables.hh"
#include "elementvolumevariables.hh"
#include "volumevariables.hh"
#include "localoperator.hh"
#include "assembler.hh"
#include "newtoncontroller.hh"
#include "dumux/geomechanics/el2p/basemodel.hh"
#include "dumux/geomechanics/el2p/indices.hh"
#include "dumux/geomechanics/el2p/localjacobian.hh"
#include "dumux/geomechanics/el2p/localresidual.hh"
#include "dumux/geomechanics/el2p/properties.hh"
#include <dumux/implicit/box/propertydefaults.hh>
#include <dumux/porousmediumflow/2p/implicit/propertydefaults.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/amgbackend.hh>

#include <dumux/material/fluidmatrixinteractions/permeabilityrutqvisttsang.hh>

namespace Dumux
{

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////

namespace Properties
{
SET_PROP(BoxElasticTwoP, NumEq) //!< set the number of equations to dim + 2
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
public:
    static const int value = dim;
};

SET_INT_PROP(BoxElasticTwoP, NumPhases, 2); //!< The number of fluid phases in the elastic 2p model is 2

//! Use the elastic local jacobian operator for the two-phase linear-elastic model
SET_TYPE_PROP(BoxElasticTwoP,
              LocalResidual,
              ElTwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxElasticTwoP, Model, DecoupledElasticModel<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(BoxElasticTwoP, ElementVolumeVariables, DecoupledElasticElementVolumeVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxElasticTwoP, VolumeVariables, DecoupledVolumeVariables<TypeTag>);

//! Set the default formulation to pWsN
SET_INT_PROP(BoxElasticTwoP,
             Formulation,
             0);

//! The indices required by the two-phase linear-elastic model
// Set the ElasticIndices to 0 (instead of 2), as we got rid of the two primary variables of
// the two-phase model
SET_PROP(BoxElasticTwoP, Indices)
{
    typedef ElTwoPIndices<TypeTag, 0, 0> type;
};

//! The FluxVariables required by the two-phase linear-elastic model
SET_TYPE_PROP(BoxElasticTwoP, FluxVariables, DecoupledElasticFluxVariables<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxElasticTwoP, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(BoxElasticTwoP, ImplicitMobilityUpwindWeight, 1.0);

//! enable gravity by default
SET_BOOL_PROP(BoxElasticTwoP, ProblemEnableGravity, true);


//! Enable evaluation of shape function gradients at the sub-control volume center by default
// Used for the computation of the pressure gradients
SET_BOOL_PROP(BoxElasticTwoP, EvalGradientsAtSCVCenter, true);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxElasticTwoP, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

SET_PROP(BoxElasticTwoP, EffectivePermeabilityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef PermeabilityRutqvistTsang<Scalar> type;
};

// SET_TYPE_PROP(BoxElasticTwoP, EffectivePermeabilityModel, PermeabilityRutqvistTsang<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Gridview)::dimension>);

// use the SuperLU linear solver by default
#if HAVE_SUPERLU
SET_TYPE_PROP(BoxElasticTwoP, LinearSolver, SuperLUBackend<TypeTag> );
#else
#warning no SuperLU detected, defaulting to ILU0BiCGSTAB. For many problems, the el2p model requires a direct solver.
SET_TYPE_PROP(BoxElasticTwoP, LinearSolver, ILU0BiCGSTABBackend<TypeTag> );
#endif

// set the grid operator
SET_PROP(BoxElasticTwoP, GridOperator)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsTrafo) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, LocalOperator) LocalOperator;
    typedef typename Dune::PDELab::ISTLMatrixBackend MatrixBackend;

    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:

    typedef Dune::PDELab::GridOperator<GridFunctionSpace,
            GridFunctionSpace,
            LocalOperator,
            MatrixBackend,
            Scalar, Scalar, Scalar,
            ConstraintsTrafo,
            ConstraintsTrafo,
            true
            > type;
};

SET_PROP(BoxElasticTwoP, JacobianMatrix)
{
private:
    //typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GFSU;
    typedef typename GFSU::template ConstraintsContainer<Scalar>::Type CU;
      //! The global assembler type
      typedef Dune::PDELab::DefaultAssembler<GFSU,GFSU,CU,CU,true> Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,Scalar>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,Scalar>::Type Range;
      //! The type of the jacobian.
    typedef typename Dune::PDELab::ISTLMatrixBackend MB;
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,Scalar>::Type Jacobian;

      //! The local assembler type
    typedef typename GET_PROP_TYPE(TypeTag, LocalOperator) LOP;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
      typedef Dune::PDELab::DefaultLocalAssembler<GridOperator,LOP,true>
      LocalAssembler;
      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSU,MB,Scalar,Scalar,Scalar,CU,CU,Assembler,LocalAssembler> Traits;
public:
    typedef typename Traits::Jacobian type;
};

SET_PROP(BoxElasticTwoP, SolutionVector)
{
private:
    //typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GFSU;
    typedef typename GFSU::template ConstraintsContainer<Scalar>::Type CU;
      //! The global assembler type
      typedef Dune::PDELab::DefaultAssembler<GFSU,GFSU,CU,CU,true> Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,Scalar>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,Scalar>::Type Range;
      //! The type of the jacobian.
    typedef typename Dune::PDELab::ISTLMatrixBackend MB;
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,Scalar>::Type Jacobian;

      //! The local assembler type
    typedef typename GET_PROP_TYPE(TypeTag, LocalOperator) LOP;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
      typedef Dune::PDELab::DefaultLocalAssembler<GridOperator,LOP,true>
      LocalAssembler;
      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSU,MB,Scalar,Scalar,Scalar,CU,CU,Assembler,LocalAssembler> Traits;
public:
    typedef typename Traits::Domain type;
};

SET_PROP(BoxElasticTwoP, GridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, DisplacementFEM) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename Dune::PDELab::EntityBlockedOrderingTag OrderingTag;
    typedef typename Dune::PDELab::ISTLVectorBackend<> VBE;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::NoConstraints Constraints;

    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints, VBE>
        ScalarGridFunctionSpace;

    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, dim, VBE, OrderingTag>
        type;

    typedef typename type::template ConstraintsContainer<Scalar>::Type
        ConstraintsTrafo;
};

SET_PROP(BoxElasticTwoP, ConstraintsTrafo)
{private:
  typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type type;
};

// set the grid function space for the sub-models
SET_TYPE_PROP(BoxElasticTwoP, Constraints, Dune::PDELab::NoConstraints);

SET_TYPE_PROP(BoxElasticTwoP, JacobianAssembler, PDELab::DecoupledElasticAssembler<TypeTag>);

SET_PROP(BoxElasticTwoP, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

SET_PROP(BoxElasticTwoP, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

SET_PROP(BoxElasticTwoP, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(BoxElasticTwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxElasticTwoP, ImplicitEnableJacobianRecycling, false);
// enable partial reassembling by default
SET_BOOL_PROP(BoxElasticTwoP, ImplicitEnablePartialReassemble, false);

SET_TYPE_PROP(BoxElasticTwoP, NewtonController, DecoupledElasticNewtonController<TypeTag>);

SET_PROP(BoxElasticTwoP, LocalOperator)
{
    typedef PDELab::DecoupledElasticLocalOperator<TypeTag> type;
};

//! use the local FEM space associated with cubes by default
SET_PROP(BoxElasticTwoP, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

/*!
 * \brief A vector of primary variables.
 */
SET_PROP(BoxElasticTwoP, PrimaryVariables)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
public:
    typedef Dune::FieldVector<Scalar, numEq> type;
};

template <class TypeTag, class MType, class VType, bool isParallel>
class ElasticTwoPSolverTraits
: public NonoverlappingSolverTraits<MType, VType, isParallel>
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
};

template <class TypeTag, class MType, class VType>
class ElasticTwoPSolverTraits<TypeTag, MType, VType, true>
: public NonoverlappingSolverTraits<MType, VType, true>
{
public:
    typedef MType JacobianMatrix;
};

//! define the traits for the AMGBackend
SET_PROP(BoxElasticTwoP, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { dofCodim = Grid::dimension,
           isNonOverlapping = true };
    enum { isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v };

    static const int numEq = isParallel ? GET_PROP_VALUE(TypeTag, NumEq)
            : GET_PROP_TYPE(TypeTag, JacobianMatrix)::block_type::rows;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef ElasticTwoPSolverTraits<TypeTag, MType, VType, isParallel> SolverTraits;
    typedef typename SolverTraits::Comm Comm;
    typedef typename SolverTraits::LinearOperator LinearOperator;
    typedef typename SolverTraits::ScalarProduct ScalarProduct;
    typedef typename SolverTraits::Smoother Smoother;
    typedef typename SolverTraits::JacobianMatrix JacobianMatrix;
};

//! The local jacobian operator
SET_TYPE_PROP(BoxElasticTwoP, LocalJacobian, ElTwoPLocalJacobian<TypeTag>);

SET_TYPE_PROP(BoxElasticTwoP, BaseModel, ElTwoPBaseModel<TypeTag>);

//! set number of equations of the mathematical model as default
SET_INT_PROP(BoxElasticTwoP, LinearSolverBlockSize, 1);

// write the stress and displacement output according to rock mechanics sign convention (compressive stresses > 0)
SET_BOOL_PROP(BoxElasticTwoP, VtkRockMechanicsSignConvention, true);

// \}
}
}

#endif
