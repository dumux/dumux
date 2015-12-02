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
#include <dune/grid/common/capabilities.hh>
#include <dune/common/version.hh>

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dumux/decoupled/common/pressureproperties.hh>
#include "linearsolverproperties.hh"

namespace Dumux
{

// Forward declaration for the property definitions
template <class TypeTag> class AMGBackend;

namespace Properties
{
//! The type traits required for using the AMG backend
NEW_PROP_TAG(AmgTraits);

template <class MType, class VType, bool isParallel>
class NonoverlappingSolverTraits
{
public:
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType,VType> LinearOperator;
    typedef Dune::SeqScalarProduct<VType> ScalarProduct;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
};

template <class MType, class VType>
class NonoverlappingSolverTraits<MType, VType, true>
{
public:
    typedef Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int> Comm;
    typedef Dune::NonoverlappingSchwarzOperator<MType,VType, VType,Comm> LinearOperator;
    typedef Dune::NonoverlappingSchwarzScalarProduct<VType,Comm> ScalarProduct;
    typedef Dune::NonoverlappingBlockPreconditioner<Comm,Dune::SeqSSOR<MType,VType, VType> > Smoother;
};

//! Box: use the non-overlapping AMG
SET_PROP(BoxModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = Grid::dimension,
        isNonOverlapping = true,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef NonoverlappingSolverTraits<MType, VType, isParallel> SolverTraits;
    typedef typename SolverTraits::Comm Comm;
    typedef typename SolverTraits::LinearOperator LinearOperator;
    typedef typename SolverTraits::ScalarProduct ScalarProduct;
    typedef typename SolverTraits::Smoother Smoother;
};

template <class MType, class VType, bool isParallel>
class OverlappingSolverTraits
{
public:
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType,VType> LinearOperator;
    typedef Dune::SeqScalarProduct<VType> ScalarProduct;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
};

template <class MType, class VType>
class OverlappingSolverTraits<MType, VType, true>
{
public:
    typedef Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int> Comm;
    typedef Dune::OverlappingSchwarzOperator<MType,VType, VType,Comm> LinearOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<VType,Comm> ScalarProduct;
    typedef Dune::BlockPreconditioner<VType,VType,Comm,Dune::SeqSSOR<MType,VType, VType> > Smoother;
};

//! Cell-centered: use the overlapping AMG
SET_PROP(CCModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
        isNonOverlapping = false,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef OverlappingSolverTraits<MType, VType, isParallel> SolverTraits;
    typedef typename SolverTraits::Comm Comm;
    typedef typename SolverTraits::LinearOperator LinearOperator;
    typedef typename SolverTraits::ScalarProduct ScalarProduct;
    typedef typename SolverTraits::Smoother Smoother;
};

//! Decoupled model: use the overlapping AMG
SET_PROP(DecoupledModel, AmgTraits)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
        isNonOverlapping = false,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef OverlappingSolverTraits<MType, VType, isParallel> SolverTraits;
    typedef typename SolverTraits::Comm Comm;
    typedef typename SolverTraits::LinearOperator LinearOperator;
    typedef typename SolverTraits::ScalarProduct ScalarProduct;
    typedef typename SolverTraits::Smoother Smoother;
};

} // namespace Properties
} // namespace Dumux
#endif
