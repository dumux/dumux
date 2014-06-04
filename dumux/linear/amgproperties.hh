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
#ifndef DUMUXAMGPROPERTIES_HH
#define DUMUXAMGPROPERTIES_HH


#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dumux/decoupled/common/pressureproperties.hh>
#include "linearsolverproperties.hh"

namespace Dumux
{

// forward declaration for the property definitions
template <class TypeTag> class AMGBackend;

namespace Properties
{
//! the PDELab finite element map used for the gridfunctionspace
NEW_PROP_TAG(AMGLocalFemMap);

//! box: use the (multi-)linear local FEM space associated with cubes by default
SET_PROP(BoxModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
 public:
    enum{
        //! \brief The codimension that the degrees of freedom are attached to.
        dofCodim = GridView::dimension
    };
};

//! cell-centered: use the element-wise constant local FEM space by default
SET_PROP(CCModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
 public:
    enum{
        //! \brief The codimension that the degrees of freedom are attached to.
        dofCodim = 0
    };
};

//! decoupled models: use the element-wise constant local FEM space by default
SET_PROP(DecoupledModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
 public:
    enum{
        //! \brief The codimension that the degrees of freedom are attached to.
        dofCodim = 0
    };
};

//! set a property JacobianMatrix also for the decoupled models
SET_PROP(DecoupledModel, JacobianMatrix)
{
public:
typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) type;
};

 //! the type of the employed PDELab backend
NEW_PROP_TAG(AMGPDELabBackend);

//! box: use the non-overlapping AMG backend
SET_PROP(BoxModel, AMGPDELabBackend)
{
 public:
    //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GridOperator> type;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        isNonOverlapping = true
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int> Comm;
    typedef Dune::NonoverlappingSchwarzOperator<MType,VType, VType,Comm> LinearOperator;
    typedef Dune::NonoverlappingSchwarzScalarProduct<VType,Comm> ScalarProduct;
    typedef Dune::NonoverlappingBlockPreconditioner<Comm,Dune::SeqSSOR<MType,VType, VType> > Smoother;
#else
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType,VType> LinearOperator;
    typedef Dune::SeqScalarProduct<VType> ScalarProduct;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
#endif
};

//! cell-centered: use the overlapping PDELab AMG backend
SET_PROP(CCModel, AMGPDELabBackend)
{
 public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        isNonOverlapping = true
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int> Comm;
    typedef Dune::OverlappingSchwarzOperator<MType,VType, VType,Comm> LinearOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<VType,Comm> ScalarProduct;
    typedef Dune::BlockPreconditioner<VType,VType,Comm,Dune::SeqSSOR<MType,VType, VType> > Smoother;
#else
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType,VType> LinearOperator;
    typedef Dune::SeqScalarProduct<VType> ScalarProduct;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
#endif
};

//! decoupled model: use the overlapping PDELab AMG backend
SET_PROP(DecoupledModel, AMGPDELabBackend)
{
    //typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
 public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        isNonOverlapping = false
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int> Comm;
    typedef Dune::OverlappingSchwarzOperator<MType,VType, VType,Comm> LinearOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<VType,Comm> ScalarProduct;
    typedef Dune::BlockPreconditioner<VType,VType,Comm,Dune::SeqSSOR<MType,VType, VType> > Smoother;
#else
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType,VType> LinearOperator;
    typedef Dune::SeqScalarProduct<VType> ScalarProduct;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
#endif
};

//! box: reset the type of solution vector to be PDELab conforming
SET_PROP(BoxModel, SolutionVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > type;
};

//! cell-centered: reset the type of solution vector to be PDELab conforming
SET_PROP(CCModel, SolutionVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > type;
};


//! decoupled model: reset the type of solution vector to be PDELab conforming
SET_PROP(DecoupledModel, PressureSolutionVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > type;
};

//! decoupled model: reset the type of solution vector to be PDELab conforming
SET_PROP(DecoupledModel, PressureRHSVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > type;
};

} // namespace Properties
} // namespace Dumux
#endif
