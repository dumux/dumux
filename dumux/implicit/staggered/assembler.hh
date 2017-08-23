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
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
#ifndef DUMUX_STAGGERED_ASSEMBLER_HH
#define DUMUX_STAGGERED_ASSEMBLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/assembler.hh>
#include <dune/istl/matrixindexset.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class StaggeredAssembler : public ImplicitAssembler<TypeTag>
{
    typedef ImplicitAssembler<TypeTag> ParentType;
    friend class ImplicitAssembler<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet::IndexType IndexType;

    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToCC CCToCCMatrixBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToFace CCToFaceMatrixBlock;

    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToFace FaceToFaceMatrixBlock;
    typedef typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToCC FaceToCCMatrixBlock;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

     /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated. We can not assume that they are
     * fully initialized, though.
     *
     * \param problem The problem object
     */
    void init(Problem& problem)
    {
        std::cout << "init(Problem& problem)" << std::endl;
        this->problemPtr_ = &problem;

        // initialize the BCRS matrix
        createMatrix_();

        // initialize the jacobian matrix with zeros
        *this->matrix_ = 0;

        // allocate the residual vector
        this->residual_[cellCenterIdx].resize(problem.model().numCellCenterDofs());
        this->residual_[faceIdx].resize(problem.model().numFaceDofs());
//         printmatrix(std::cout, matrix(), "", "");
    }


private:

    // Construct the multitype matrix for the global jacobian
    void createMatrix_()
    {
        // create the multitype matrix
        this->matrix_ = std::make_shared<JacobianMatrix>();

        // get sub matrix sizes
        const auto cellCenterSize = this->problem_().model().numCellCenterDofs();
        const auto faceSize = this->problem_().model().numFaceDofs();

        // allocate the sub matrices (BCRS matrices)
        auto A11 = CCToCCMatrixBlock(cellCenterSize, cellCenterSize, CCToCCMatrixBlock::random);
        auto A12 = CCToFaceMatrixBlock(cellCenterSize, faceSize, CCToFaceMatrixBlock::random);

        auto A21 = FaceToCCMatrixBlock(faceSize, cellCenterSize, FaceToCCMatrixBlock::random);
        auto A22 = FaceToFaceMatrixBlock(faceSize, faceSize, FaceToFaceMatrixBlock::random);

        setupMatrices_(A11, A12, A21, A22);

        (*this->matrix_)[cellCenterIdx][cellCenterIdx] = A11;
        (*this->matrix_)[cellCenterIdx][faceIdx] = A12;
        (*this->matrix_)[faceIdx][cellCenterIdx] = A21;
        (*this->matrix_)[faceIdx][faceIdx] = A22;

//         printmatrix(std::cout, A11, "A11", "");
//         printmatrix(std::cout, A12, "A12", "");
//         printmatrix(std::cout, A21, "A21", "");
//         printmatrix(std::cout, A22, "A22", "");
    }

    void setupMatrices_(CCToCCMatrixBlock &A11, CCToFaceMatrixBlock &A12,
                     FaceToCCMatrixBlock &A21, FaceToFaceMatrixBlock &A22)
    {
        // get sub matrix sizes
        const auto numDofsCC = this->problem_().model().numCellCenterDofs();
        const auto numDofsFace = this->problem_().model().numFaceDofs();

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPatternA11;
        Dune::MatrixIndexSet occupationPatternA12;
        Dune::MatrixIndexSet occupationPatternA21;
        Dune::MatrixIndexSet occupationPatternA22;
        occupationPatternA11.resize(numDofsCC, numDofsCC);
        occupationPatternA12.resize(numDofsCC, numDofsFace);
        occupationPatternA21.resize(numDofsFace, numDofsCC);
        occupationPatternA22.resize(numDofsFace, numDofsFace);

        const auto& assemblyMap = this->problem_().model().localJacobian().assemblyMap();

        for (const auto& element : elements(this->gridView_()))
        {
            // the global index of the element at hand
            const auto ccGlobalI = this->elementMapper_().index(element);

            for (auto&& ccGlobalJ : assemblyMap(cellCenterIdx, cellCenterIdx, ccGlobalI))
                occupationPatternA11.add(ccGlobalI, ccGlobalJ);
            for (auto&& faceGlobalJ : assemblyMap(cellCenterIdx, faceIdx, ccGlobalI))
                occupationPatternA12.add(ccGlobalI, faceGlobalJ);

            auto fvGeometry = localView(this->problem_().model().fvGridGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceGlobalI = scvf.dofIndex();
                for (auto&& ccGlobalJ : assemblyMap(faceIdx, cellCenterIdx, scvf.index()))
                    occupationPatternA21.add(faceGlobalI, ccGlobalJ);
                for (auto&& faceGlobalJ : assemblyMap(faceIdx, faceIdx, scvf.index()))
                    occupationPatternA22.add(faceGlobalI, faceGlobalJ);
            }
        }
        // export patterns to matrices
        occupationPatternA11.exportIdx(A11);
        occupationPatternA12.exportIdx(A12);
        occupationPatternA21.exportIdx(A21);
        occupationPatternA22.exportIdx(A22);
    }


};

} // namespace Dumux

#endif
