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
 * \ingroup MixedDimension
 * \brief An assembler for the global Jacobian matrix for models of mixed dimension.
 *        We assume a bulk domain and a lower dimensional domain on separate grids.
 */
#ifndef DUMUX_MIXEDDIMENSION_ASSEMBLER_HH
#define DUMUX_MIXEDDIMENSION_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/mixeddimension/properties.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \brief An assembler for the global Jacobian matrix for models of mixed dimension.
 *        We assume a bulk domain and a lower dimensional domain on separate grids.
 */
template<class TypeTag>
class MixedDimensionAssembler
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);

    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SubProblemBlockIndices = typename GET_PROP(TypeTag, SubProblemBlockIndices);
    using BulkMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulk;
    using BulkCouplingMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulkCoupling;
    using LowDimMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDim;
    using LowDimCouplingMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDimCoupling;

    typename SubProblemBlockIndices::BulkIdx bulkIdx;
    typename SubProblemBlockIndices::LowDimIdx lowDimIdx;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

public:
    // copying the jacobian assembler is not a good idea
    MixedDimensionAssembler(const MixedDimensionAssembler &) = delete;

    MixedDimensionAssembler() : problemPtr_(nullptr) {}

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

        // initialize the multitype matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix with zeros
        *matrix_ = 0.0;

        // allocate the residual vector (right-hand-side)
        residual_[bulkIdx].resize(problem.model().bulkNumDofs());
        residual_[lowDimIdx].resize(problem.model().lowDimNumDofs());
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

        // assemble the elements of the bulk problem
        for (const auto& element : elements(problem_().bulkGridView()))
            problem_().model().bulkLocalJacobian().assemble(element,
                                                            (*matrix_)[bulkIdx][bulkIdx],
                                                            (*matrix_)[bulkIdx][lowDimIdx],
                                                            residual_[bulkIdx]);

        // assemble the elements of the lowdim problem
        for (const auto& element : elements(problem_().lowDimGridView()))
            problem_().model().lowDimLocalJacobian().assemble(element,
                                                              (*matrix_)[lowDimIdx][lowDimIdx],
                                                              (*matrix_)[lowDimIdx][bulkIdx],
                                                              residual_[lowDimIdx]);
    }

    template <class T = BulkProblemTypeTag, typename std::enable_if<((!GET_PROP_VALUE(T, ImplicitIsBox))), bool>::type = 0>
    void buildBulkMatrixBlocks_(BulkMatrixBlock& m, BulkCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet bulkPattern, bulkCouplingPattern;
        bulkPattern.resize(m.N(), m.M());
        bulkCouplingPattern.resize(cm.N(), cm.M());

        const auto& assemblyMap = this->problem_().bulkProblem().model().localJacobian().assemblyMap();
        for (const auto& element : elements(problem_().bulkGridView()))
        {
            const auto globalI = problem_().model().bulkElementMapper().index(element);

            bulkPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap[globalI])
                bulkPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            const auto& additionalDofDependencies = problem_().bulkProblem().getAdditionalDofDependencies(globalI);
            for (auto globalJ : additionalDofDependencies)
                bulkPattern.add(globalI, globalJ);

            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            for (auto globalJ : couplingStencil)
                bulkCouplingPattern.add(globalI, globalJ);
        }

        // export occupation patterns to matrices
        bulkPattern.exportIdx(m);
        bulkCouplingPattern.exportIdx(cm);
    }

    template <class T = LowDimProblemTypeTag, typename std::enable_if<((!GET_PROP_VALUE(T, ImplicitIsBox))), bool>::type = 0>
    void buildLowDimMatrixBlocks_(LowDimMatrixBlock& m, LowDimCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet lowDimPattern, lowDimCouplingPattern;
        lowDimPattern.resize(m.N(), m.M());
        lowDimCouplingPattern.resize(cm.N(), cm.M());

        const auto& assemblyMap = this->problem_().lowDimProblem().model().localJacobian().assemblyMap();
        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            const auto globalI = problem_().model().lowDimElementMapper().index(element);

            lowDimPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap[globalI])
                lowDimPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            const auto& additionalDofDependencies = problem_().lowDimProblem().getAdditionalDofDependencies(globalI);
            for (auto globalJ : additionalDofDependencies)
                lowDimPattern.add(globalI, globalJ);

            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            for (auto globalJ : couplingStencil)
                lowDimCouplingPattern.add(globalI, globalJ);
        }

        // export occupation patterns to matrices
        lowDimPattern.exportIdx(m);
        lowDimCouplingPattern.exportIdx(cm);
    }

    template <class T = BulkProblemTypeTag, typename std::enable_if<((GET_PROP_VALUE(T, ImplicitIsBox))), bool>::type = 0>
    void buildBulkMatrixBlocks_(BulkMatrixBlock& m, BulkCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet bulkPattern, bulkCouplingPattern;
        bulkPattern.resize(m.N(), m.M());
        bulkCouplingPattern.resize(cm.N(), cm.M());

        for (const auto& element : elements(problem_().bulkGridView()))
        {
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int vIdx = 0; vIdx < element.subEntities(bulkDim); ++vIdx)
            {
                const auto globalI = problem_().model().bulkVertexMapper().subIndex(element, vIdx, bulkDim);
                for (unsigned int vIdx2 = vIdx; vIdx2 < element.subEntities(bulkDim); ++vIdx2)
                {
                    const auto globalJ = problem_().model().bulkVertexMapper().subIndex(element, vIdx2, bulkDim);
                    bulkPattern.add(globalI, globalJ);
                    bulkPattern.add(globalJ, globalI);
                }

                for (auto&& globalJ : couplingStencil)
                    bulkCouplingPattern.add(globalI, globalJ);

                //TODO: additionalDofDependencies
            }
        }

        // export occupation patterns to matrices
        bulkPattern.exportIdx(m);
        bulkCouplingPattern.exportIdx(cm);
    }

    template <class T = LowDimProblemTypeTag, typename std::enable_if<((GET_PROP_VALUE(T, ImplicitIsBox))), bool>::type = 0>
    void buildLowDimMatrixBlocks_(LowDimMatrixBlock& m, LowDimCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet lowDimPattern, lowDimCouplingPattern;
        lowDimPattern.resize(m.N(), m.M());
        lowDimCouplingPattern.resize(cm.N(), cm.M());

        for (const auto& element : elements(problem_().lowDimGridView()))
        {
            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

            for (unsigned int vIdx = 0; vIdx < element.subEntities(lowDimDim); ++vIdx)
            {
                const auto globalI = problem_().model().lowDimVertexMapper().subIndex(element, vIdx, lowDimDim);
                for (unsigned int vIdx2 = vIdx; vIdx2 < element.subEntities(lowDimDim); ++vIdx2)
                {
                    const auto globalJ = problem_().model().lowDimVertexMapper().subIndex(element, vIdx2, lowDimDim);
                    lowDimPattern.add(globalI, globalJ);
                    lowDimPattern.add(globalJ, globalI);
                }

                for (auto&& globalJ : couplingStencil)
                    lowDimCouplingPattern.add(globalI, globalJ);

                //TODO: additionalDofDependencies
            }
        }

        // export occupation patterns to matrices
        lowDimPattern.exportIdx(m);
        lowDimCouplingPattern.exportIdx(cm);
    }

    // Construct the multitype matrix for the global jacobian
    void createMatrix_()
    {
        // create the multitype matrix
        matrix_ = std::make_shared<JacobianMatrix>();

        // get sub matrix sizes
        const auto bulkSize = problem_().model().bulkNumDofs();
        const auto lowDimSize = problem_().model().lowDimNumDofs();

        // allocate the sub matrices
        auto A11 = BulkMatrixBlock(bulkSize, bulkSize, BulkMatrixBlock::random);
        auto A12 = BulkCouplingMatrixBlock(bulkSize, lowDimSize, BulkCouplingMatrixBlock::random);
        auto A22 = LowDimMatrixBlock(lowDimSize, lowDimSize, LowDimMatrixBlock::random);
        auto A21 = LowDimCouplingMatrixBlock(lowDimSize, bulkSize, LowDimCouplingMatrixBlock::random);

        buildBulkMatrixBlocks_(A11, A12);
        buildLowDimMatrixBlocks_(A22, A21);

        (*matrix_)[bulkIdx][bulkIdx] = A11;
        (*matrix_)[bulkIdx][lowDimIdx] = A12;
        (*matrix_)[lowDimIdx][bulkIdx] = A21;
        (*matrix_)[lowDimIdx][lowDimIdx] = A22;
    }

    // assemble a bulk element
    // void assembleElement_(const BulkElement &element)
    // {
    //     problem_().model().bulkLocalJacobian().assemble(element);

    //     if (!bulkIsBox)
    //     {
    //         int globalI = problem_().model().bulkElementMapper().index(element);

    //         // update the right hand side
    //         residual_[bulkIdx][globalI] =  problem_().model().bulkLocalJacobian().residual(0);
    //         for (int j = 0; j < residual_[bulkIdx][globalI].dimension; ++j)
    //             assert(std::isfinite(residual_[bulkIdx][globalI][j]));

    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         int j = 0;
    //         for (auto&& globalJ : stencil)
    //             (*matrix_)[bulkIdx][bulkIdx][globalI][globalJ] = problem_().model().bulkLocalJacobian().mat(bulkIdx, 0, j++);

    //         j = 0;
    //         for (auto&& globalJ : couplingStencil)
    //             (*matrix_)[bulkIdx][lowDimIdx][globalI][globalJ] = problem_().model().bulkLocalJacobian().mat(lowDimIdx, 0, j++);
    //     }
    //     else
    //     {
    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         for (unsigned int scvIdx = 0; scvIdx < element.subEntities(bulkDim); ++scvIdx)
    //         {
    //             auto globalI = problem_().model().bulkVertexMapper().subIndex(element, scvIdx, bulkDim);

    //             // update the right hand side
    //             residual_[bulkIdx][globalI] += problem_().model().bulkLocalJacobian().residual(scvIdx);
    //             for (int j = 0; j < residual_[bulkIdx][globalI].dimension; ++j)
    //                 assert(std::isfinite(residual_[bulkIdx][globalI][j]));

    //             int j = 0;
    //             for (auto&& globalJ : stencil)
    //                 (*matrix_)[bulkIdx][bulkIdx][globalI][globalJ] += problem_().model().bulkLocalJacobian().mat(bulkIdx, scvIdx, j++);

    //             j = 0;
    //             for (auto&& globalJ : couplingStencil)
    //                 (*matrix_)[bulkIdx][lowDimIdx][globalI][globalJ] += problem_().model().bulkLocalJacobian().mat(lowDimIdx, scvIdx, j++);
    //         }
    //     }
    // }

    // // assemble a bulk element
    // void assembleElement_(const LowDimElement &element)
    // {
    //     problem_().model().lowDimLocalJacobian().assemble(element);

    //     if (!lowDimIsBox)
    //     {
    //         int globalI = problem_().model().lowDimElementMapper().index(element);

    //         // update the right hand side
    //         residual_[lowDimIdx][globalI] =  problem_().model().lowDimLocalJacobian().residual(0);
    //         for (int j = 0; j < residual_[lowDimIdx][globalI].dimension; ++j)
    //             assert(std::isfinite(residual_[lowDimIdx][globalI][j]));

    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         int j = 0;
    //         for (auto&& globalJ : stencil)
    //             (*matrix_)[lowDimIdx][lowDimIdx][globalI][globalJ] = problem_().model().lowDimLocalJacobian().mat(lowDimIdx, 0, j++);

    //         j = 0;
    //         for (auto&& globalJ : couplingStencil)
    //             (*matrix_)[lowDimIdx][bulkIdx][globalI][globalJ] = problem_().model().lowDimLocalJacobian().mat(bulkIdx, 0, j++);
    //     }
    //     else
    //     {
    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         for (unsigned int scvIdx = 0; scvIdx < element.subEntities(lowDimDim); ++scvIdx)
    //         {
    //             auto globalI = problem_().model().lowDimVertexMapper().subIndex(element, scvIdx, lowDimDim);

    //             // update the right hand side
    //             residual_[lowDimIdx][globalI] +=  problem_().model().lowDimLocalJacobian().residual(scvIdx);
    //             for (int j = 0; j < residual_[lowDimIdx][globalI].dimension; ++j)
    //                 assert(std::isfinite(residual_[lowDimIdx][globalI][j]));

    //             int j = 0;
    //             for (auto&& globalJ : stencil)
    //                 (*matrix_)[lowDimIdx][lowDimIdx][globalI][globalJ] += problem_().model().lowDimLocalJacobian().mat(lowDimIdx, scvIdx, j++);

    //             j = 0;
    //             for (auto&& globalJ : couplingStencil)
    //                 (*matrix_)[lowDimIdx][bulkIdx][globalI][globalJ] += problem_().model().lowDimLocalJacobian().mat(bulkIdx, scvIdx, j++);
    //         }
    //     }
    // }

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
