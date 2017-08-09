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
 * \ingroup MultiDomain
 * \brief An assembler for the global Jacobian matrix for models of equal dimension.
 *        We assume two domains on on separate grids.
 */
#ifndef DUMUX_MULTIDOMAIN_ASSEMBLER_FOR_STAGGERED_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLER_FOR_STAGGERED_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/multidomain/properties.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief An assembler for the global Jacobian matrix for models of equal dimension.
 *        We assume two domains on separate grids.
 */
template<class TypeTag>
class MultiDomainAssemblerForStaggered
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);

    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SubProblemBlockIndices = typename GET_PROP(TypeTag, SubProblemBlockIndices);

    using CCToCCMatrixBlock = typename GET_PROP(StokesProblemTypeTag, JacobianMatrix)::MatrixBlockCCToCC;
    using CCToFaceMatrixBlock = typename GET_PROP(StokesProblemTypeTag, JacobianMatrix)::MatrixBlockCCToFace;

    using FaceToFaceMatrixBlock = typename GET_PROP(StokesProblemTypeTag, JacobianMatrix)::MatrixBlockFaceToFace;
    using FaceToCCMatrixBlock = typename GET_PROP(StokesProblemTypeTag, JacobianMatrix)::MatrixBlockFaceToCC;

    using DarcyMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::DarcyJacobianMatrix;

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;


    typename SubProblemBlockIndices::StokesIdx stokesIdx;
    typename SubProblemBlockIndices::DarcyIdx darcyIdx;

    typename GET_PROP(TypeTag, JacobianMatrix)::localDarcyIdx localDarcyIdx;

    using MatrixBlockStokesCCToDarcy = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockStokesCCToDarcy;
    using MatrixBlockStokesFaceToDarcy = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockStokesFaceToDarcy;
    using MatrixBlockDarcyToStokesCC = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockDarcyToStokesCC;
    using MatrixBlockDarcyToStokesFace = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockDarcyToStokesFace;

    enum {
        dim = StokesGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

public:
    // copying the jacobian assembler is not a good idea
    MultiDomainAssemblerForStaggered(const MultiDomainAssemblerForStaggered &) = delete;

    MultiDomainAssemblerForStaggered() : problemPtr_(nullptr) {}

    /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated :. We can not assume that they are
     * fully initialized, though.
     *
     * \param problem The problem object
     */
    void init(Problem& problem)
    {
        // save problem pointer
        problemPtr_ = &problem;

        std::cout << "in init(Problem& problem) of multidomain/staggeredgrid/assembler" << std::endl;

        // initialize the multitype matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix with zeros
        *matrix_ = 0.0;

        // allocate the residual vector (right-hand-side)
        residual_[stokesIdx][cellCenterIdx].resize(problem.stokesProblem().model().numCellCenterDofs());
        residual_[stokesIdx][faceIdx].resize(problem.stokesProblem().model().numFaceDofs());
        residual_[darcyIdx].resize(problem.model().darcyNumDofs());



//         printmatrix(std::cout, (*matrix_)[stokesIdx][stokesIdx][cellCenterIdx][cellCenterIdx], "A11", "");
//         printmatrix(std::cout, (*matrix_)[stokesIdx][stokesIdx][cellCenterIdx][faceIdx], "A12", "");
//         printmatrix(std::cout, (*matrix_)[stokesIdx][stokesIdx][faceIdx][cellCenterIdx], "A21", "");
//         printmatrix(std::cout, (*matrix_)[stokesIdx][stokesIdx][faceIdx][faceIdx], "A22", "");
//         printmatrix(std::cout,(*matrix_)[stokesIdx][darcyIdx][cellCenterIdx][localDarcyIdx], "A13", "");
//         printmatrix(std::cout,(*matrix_)[stokesIdx][darcyIdx][faceIdx][localDarcyIdx], "A23", "");
//         printmatrix(std::cout,(*matrix_)[darcyIdx][darcyIdx][localDarcyIdx][localDarcyIdx], "A33", "");
//         printmatrix(std::cout,(*matrix_)[darcyIdx][stokesIdx][localDarcyIdx][cellCenterIdx], "A31", "");
//         printmatrix(std::cout,(*matrix_)[darcyIdx][stokesIdx][localDarcyIdx][faceIdx], "A32", "");
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool succeeded;
        try
        {
            asImp_().assemble_();
            succeeded = true;
        }
        catch (Dumux::NumericalProblem &e)
        {
            std::cout << " caught an exception while assembling:" << e.what() << std::endl;
            succeeded = false;
        }

        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const JacobianMatrix& matrix() const
    { return *matrix_; }

    JacobianMatrix& matrix()
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return residual_; }

    SolutionVector& residual()
    { return residual_; }

protected:
    // reset the global linear system of equations.
    void resetSystem_()
    {
        // reset the right hand side.
        residual_ = 0.0;

        // reset the matrix
        *matrix_ = 0.0;
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // assemble the elements of the Stokes problem
        for (const auto& element : elements(problem_().stokesGridView()))
            problem_().model().stokesLocalJacobian().assemble(element,
                                                            (*matrix_)[stokesIdx][stokesIdx],
                                                            (*matrix_)[stokesIdx][darcyIdx],
                                                            residual_[stokesIdx]);

        // assemble the elements of the Darcy problem
        for (const auto& element : elements(problem_().darcyGridView()))
            problem_().model().darcyLocalJacobian().assemble(element,
                                                              (*matrix_)[darcyIdx][darcyIdx][localDarcyIdx][localDarcyIdx],
                                                              (*matrix_)[darcyIdx][stokesIdx],
                                                              residual_[darcyIdx]);
    }

    // Construct the multitype matrix for the global jacobian
    void createMatrix_()
    {
        // create the multitype matrix
        matrix_ = std::make_shared<JacobianMatrix>();

        // get sub matrix sizes for staggered model
        const auto cellCenterSize = this->problem_().stokesProblem().model().numCellCenterDofs();
        const auto faceSize = this->problem_().stokesProblem().model().numFaceDofs();
        const auto darcySize = problem_().model().darcyNumDofs();

        // allocate the sub matrices (BCRS matrices)
        auto A11 = CCToCCMatrixBlock(cellCenterSize, cellCenterSize, CCToCCMatrixBlock::random);
        auto A12 = CCToFaceMatrixBlock(cellCenterSize, faceSize, CCToFaceMatrixBlock::random);

        auto A21 = FaceToCCMatrixBlock(faceSize, cellCenterSize, FaceToCCMatrixBlock::random);
        auto A22 = FaceToFaceMatrixBlock(faceSize, faceSize, FaceToFaceMatrixBlock::random);

        auto A13 = MatrixBlockStokesCCToDarcy(cellCenterSize, darcySize, MatrixBlockStokesCCToDarcy::random);
        auto A23 = MatrixBlockStokesFaceToDarcy(faceSize, darcySize, MatrixBlockStokesFaceToDarcy::random);
        setupStokesMatrices_(A11, A12, A21, A22, A13, A23);

        (*matrix_)[stokesIdx][stokesIdx][cellCenterIdx][cellCenterIdx] = A11;
        (*matrix_)[stokesIdx][stokesIdx][cellCenterIdx][faceIdx] = A12;
        (*matrix_)[stokesIdx][stokesIdx][faceIdx][cellCenterIdx] = A21;
        (*matrix_)[stokesIdx][stokesIdx][faceIdx][faceIdx] = A22;

        (*matrix_)[stokesIdx][darcyIdx][cellCenterIdx][localDarcyIdx] = A13;
        (*matrix_)[stokesIdx][darcyIdx][faceIdx][localDarcyIdx] = A23;

        // now treat the Darcy problem
        auto A31 = MatrixBlockDarcyToStokesCC(darcySize, cellCenterSize, MatrixBlockDarcyToStokesCC::random);
        auto A32 = MatrixBlockDarcyToStokesFace(darcySize, faceSize, MatrixBlockDarcyToStokesFace::random);
        auto A33 = DarcyMatrixBlock(darcySize, darcySize, DarcyMatrixBlock::random);

        setupDarcyMatrices_(A33, A31, A32);
        (*matrix_)[darcyIdx][darcyIdx][localDarcyIdx][localDarcyIdx] = A33;
        (*matrix_)[darcyIdx][stokesIdx][localDarcyIdx][cellCenterIdx] = A31;
        (*matrix_)[darcyIdx][stokesIdx][localDarcyIdx][faceIdx] = A32;
    }

    void setupStokesMatrices_(CCToCCMatrixBlock& A11,
                            CCToFaceMatrixBlock& A12,
                            FaceToCCMatrixBlock& A21,
                            FaceToFaceMatrixBlock& A22,
                            MatrixBlockStokesCCToDarcy& A13,
                            MatrixBlockStokesFaceToDarcy& A23)
    {
        // get sub matrix sizes
        const auto numDofsCC = problem_().stokesProblem().model().numCellCenterDofs();
        const auto numDofsFace = problem_().stokesProblem().model().numFaceDofs();

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPatternA11;
        Dune::MatrixIndexSet occupationPatternA12;
        Dune::MatrixIndexSet occupationPatternA21;
        Dune::MatrixIndexSet occupationPatternA22;
        occupationPatternA11.resize(numDofsCC, numDofsCC);
        occupationPatternA12.resize(numDofsCC, numDofsFace);
        occupationPatternA21.resize(numDofsFace, numDofsCC);
        occupationPatternA22.resize(numDofsFace, numDofsFace);

        Dune::MatrixIndexSet occupationPatternA13, occupationPatternA23;
        occupationPatternA13.resize(A13.N(), A13.M());
        occupationPatternA23.resize(A23.N(), A23.M());

        const auto& assemblyMap = problem_().stokesProblem().model().localJacobian().assemblyMap();

        for (const auto& element : elements(problem_().stokesProblem().gridView()))
        {
            // the global index of the element at hand
            const auto ccGlobalI = problem_().stokesProblem().elementMapper().index(element);

            const auto& couplingCCStencil = problem_().couplingManager().couplingStencil(element);

            for (auto&& ccGlobalJ : assemblyMap(cellCenterIdx, cellCenterIdx, ccGlobalI))
                occupationPatternA11.add(ccGlobalI, ccGlobalJ);
            for (auto&& faceGlobalJ : assemblyMap(cellCenterIdx, faceIdx, ccGlobalI))
                occupationPatternA12.add(ccGlobalI, faceGlobalJ);

            for (auto&& globalJ : couplingCCStencil)
                occupationPatternA13.add(ccGlobalI, globalJ);

            auto fvGeometry = localView(problem_().stokesProblem().model().globalFvGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceGlobalI = scvf.dofIndex();
                for (auto&& ccGlobalJ : assemblyMap(faceIdx, cellCenterIdx, scvf.index()))
                    occupationPatternA21.add(faceGlobalI, ccGlobalJ);
                for (auto&& faceGlobalJ : assemblyMap(faceIdx, faceIdx, scvf.index()))
                    occupationPatternA22.add(faceGlobalI, faceGlobalJ);

                const auto& couplingFaceStencil = problem_().couplingManager().couplingStencil(scvf);
                for (auto&& globalJ : couplingFaceStencil)
                    occupationPatternA23.add(faceGlobalI, globalJ);
            }
        }
        // export patterns to matrices
        occupationPatternA11.exportIdx(A11);
        occupationPatternA12.exportIdx(A12);
        occupationPatternA21.exportIdx(A21);
        occupationPatternA22.exportIdx(A22);

        occupationPatternA13.exportIdx(A13);
        occupationPatternA23.exportIdx(A23);
    }

    void setupDarcyMatrices_(DarcyMatrixBlock& m,
                              MatrixBlockDarcyToStokesCC& darcyToCC,
                              MatrixBlockDarcyToStokesFace& darcyToFace)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet darcyPattern, darcyToCCPattern, darcyToFacePattern;
        darcyPattern.resize(m.N(), m.M());
        darcyToCCPattern.resize(darcyToCC.N(), darcyToCC.M());
        darcyToFacePattern.resize(darcyToFace.N(), darcyToFace.M());

        for (const auto& element : elements(problem_().darcyGridView()))
        {
            const auto& couplingCCStencil = problem_().couplingManager().couplingStencil(element, cellCenterIdx);
            const auto& couplingFaceStencil = problem_().couplingManager().couplingStencil(element, faceIdx);

            for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
            {
                const auto globalI = problem_().model().darcyVertexMapper().subIndex(element, vIdx, dim);
                for (unsigned int vIdx2 = vIdx; vIdx2 < element.subEntities(dim); ++vIdx2)
                {
                    const auto globalJ = problem_().model().darcyVertexMapper().subIndex(element, vIdx2, dim);
                    darcyPattern.add(globalI, globalJ);
                    darcyPattern.add(globalJ, globalI);
                }

                for (auto&& globalJ : couplingCCStencil)
                    darcyToCCPattern.add(globalI, globalJ);
                for (auto&& globalJ : couplingFaceStencil)
                    darcyToFacePattern.add(globalI, globalJ);
            }
        }

        // export occupation patterns to matrices
        darcyPattern.exportIdx(m);
        darcyToCCPattern.exportIdx(darcyToCC);
        darcyToFacePattern.exportIdx(darcyToFace);
    }



    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }

    // the multidimension problem
    Problem *problemPtr_;

    // the jacobian matrix
    std::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    SolutionVector residual_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
