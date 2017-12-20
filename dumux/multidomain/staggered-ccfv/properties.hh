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
 * \ingroup MultiDomainModel
 * \brief Base properties for equal dimensional models
 *        including a staggered grid discretization for the stokes model
 */
#ifndef DUMUX_STAGGERED_MULTIDOMAIN_PROPERTIES_HH
#define DUMUX_STAGGERED_MULTIDOMAIN_PROPERTIES_HH

#include <dumux/implicit/staggered/properties.hh>

#include <dumux/multidomain/subproblemproperties.hh>
#include <dumux/multidomain/staggered-ccfv/model.hh>
#include <dumux/multidomain/staggered-ccfv/assembler.hh>
#include <dumux/multidomain/staggered-ccfv/newtoncontroller.hh>
#include <dumux/multidomain/staggered-ccfv/subproblemlocaljacobian.hh>

namespace Dumux
{
namespace Properties
{
NEW_TYPE_TAG(CouplingStokesStaggeredModel, INHERITS_FROM(StaggeredModel));

SET_TYPE_PROP(CouplingStokesStaggeredModel, StokesLocalJacobian, StokesLocalJacobianForStaggered<TypeTag>);
SET_TYPE_PROP(CouplingStokesStaggeredModel, DarcyLocalJacobian, DarcyLocalJacobianForStaggered<TypeTag>);

SET_PROP(CouplingStokesStaggeredModel, SolutionVector)
{
private:
    using StokesSolutionVector = typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag), SolutionVector);
    using DarcySolutionVector = typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag), SolutionVector);
public:
    using type = typename Dune::MultiTypeBlockVector<StokesSolutionVector, DarcySolutionVector>;
};


//! Set the type of a global jacobian matrix from the solution types
SET_PROP(CouplingStokesStaggeredModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);
    enum {
        numEqStokesCC = GET_PROP_VALUE(StokesProblemTypeTag, NumEqCellCenter),
        numEqStokesFace = GET_PROP_VALUE(StokesProblemTypeTag, NumEqFace),
        numEqDarcy = GET_PROP_VALUE(DarcyProblemTypeTag, NumEq)
    };

public:

    // the problems' actual jacobian matrices
    using StokesJacobianMatrix = typename GET_PROP_TYPE(StokesProblemTypeTag, JacobianMatrix);
    using DarcyJacobianMatrix = typename GET_PROP_TYPE(DarcyProblemTypeTag, JacobianMatrix);

    // the coupling matrices
    using MatrixLittleBlockStokesCCToDarcy = typename  Dune::FieldMatrix<Scalar, numEqStokesCC, numEqDarcy>;
    using MatrixBlockStokesCCToDarcy = typename Dune::BCRSMatrix<MatrixLittleBlockStokesCCToDarcy>;

    using MatrixLittleBlockStokesFaceToDarcy = typename  Dune::FieldMatrix<Scalar, numEqStokesFace, numEqDarcy>;
    using MatrixBlockStokesFaceToDarcy = typename Dune::BCRSMatrix<MatrixLittleBlockStokesFaceToDarcy>;

    using MatrixLittleBlockDarcyToStokesCC = typename  Dune::FieldMatrix<Scalar, numEqDarcy, numEqStokesCC>;
    using MatrixBlockDarcyToStokesCC = typename Dune::BCRSMatrix<MatrixLittleBlockDarcyToStokesCC>;

    using MatrixLittleBlockDarcyToStokesFace = typename  Dune::FieldMatrix<Scalar, numEqDarcy, numEqStokesFace>;
    using MatrixBlockDarcyToStokesFace = typename Dune::BCRSMatrix<MatrixLittleBlockDarcyToStokesFace>;

    // the final four matrix blocks for the global Jacobian
    using StokesBlock = StokesJacobianMatrix;

    using DarcyBlock = typename Dune::MultiTypeBlockMatrix<typename Dune::MultiTypeBlockVector<DarcyJacobianMatrix>>;

    using StokesCellCenterToDarcyRow = typename Dune::MultiTypeBlockVector<MatrixBlockStokesCCToDarcy>;
    using StokesFaceToDarcyRow = typename Dune::MultiTypeBlockVector<MatrixBlockStokesFaceToDarcy>;
    using StokesToDarcyCouplingBlock = typename Dune::MultiTypeBlockMatrix<StokesCellCenterToDarcyRow, StokesFaceToDarcyRow>;

    using DarcyToStokesRow = typename Dune::MultiTypeBlockVector<MatrixBlockDarcyToStokesCC, MatrixBlockDarcyToStokesFace>;
    using DarcyToStokesCouplingBlock = typename Dune::MultiTypeBlockMatrix<DarcyToStokesRow>;

    // the two rows for the global Jacobian
    using StokesRow = typename Dune::MultiTypeBlockVector<StokesBlock, StokesToDarcyCouplingBlock>;
    using DarcyRow = typename Dune::MultiTypeBlockVector<DarcyToStokesCouplingBlock, DarcyBlock>;

    // the actual global Jacobian
    using type = Dune::MultiTypeBlockMatrix<StokesRow, DarcyRow>;

    // local index to acces the correct blocks within the darcy matrix
    using localDarcyIdx = Dune::index_constant<0>;
};

//! Set the BaseModel to MultiDomainModel
SET_TYPE_PROP(CouplingStokesStaggeredModel, Model, MultiDomainModelForStaggered<TypeTag>);
SET_TYPE_PROP(CouplingStokesStaggeredModel, JacobianAssembler, MultiDomainAssemblerForStaggered<TypeTag>);
SET_TYPE_PROP(CouplingStokesStaggeredModel, NewtonController, MultiDomainNewtonControllerForStaggered<TypeTag>);

}//end namespace Properties

}//end namespace Dumux

#endif
