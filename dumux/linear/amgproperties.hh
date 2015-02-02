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
 * \ingroup Linear
 *
 * \brief Defines some fundamental properties for the AMG backend.
 */
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
//! the type traits required for using the AMG backend
NEW_PROP_TAG(AmgTraits);

//! box: use the non-overlapping AMG
SET_PROP(BoxModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = GridView::dimension,
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

//! cell-centered: use the overlapping AMG
SET_PROP(CCModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
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

//! decoupled model: use the overlapping AMG
SET_PROP(DecoupledModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) JacobianMatrix;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
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

} // namespace Properties
} // namespace Dumux
#endif
