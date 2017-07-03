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
 * \brief An assembler for the global Jacobian matrix for coupled models.
 *        We assume a Stokes domain and a Darcy domain on separate grids.
 */
#ifndef DUMUX_MIXEDDIMENSION_ASSEMBLER_HH
#define DUMUX_MIXEDDIMENSION_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/multidomain/properties.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief An assembler for the global Jacobian matrix for coupled models.
 *        We assume a Stokes domain and a Darcy domain on separate grids.
 */
template<class TypeTag>
class MultiDomainAssembler
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
    using StokesMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockStokes;
    using StokesCouplingMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockStokesCoupling;
    using DarcyMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockDarcy;
    using DarcyCouplingMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockDarcyCoupling;

    typename SubProblemBlockIndices::StokesIdx stokesIdx;
    typename SubProblemBlockIndices::DarcyIdx darcyIdx;

    enum {
        dim = StokesGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

public:
    // copying the jacobian assembler is not a good idea
    MultiDomainAssembler(const MultiDomainAssembler &) = delete;

    MultiDomainAssembler() : problemPtr_(nullptr) {}

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
        residual_[stokesIdx].resize(problem.model().stokesNumDofs());
        residual_[darcyIdx].resize(problem.model().darcyNumDofs());
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

        // assemble the elements of the stokes problem
        for (const auto& element : elements(problem_().stokesGridView()))
            problem_().model().stokesLocalJacobian().assemble(element,
                                                            (*matrix_)[stokesIdx][stokesIdx],
                                                            (*matrix_)[stokesIdx][darcyIdx],
                                                            residual_[stokesIdx]);

        // assemble the elements of the lowdim problem
        for (const auto& element : elements(problem_().darcyGridView()))
            problem_().model().darcyLocalJacobian().assemble(element,
                                                              (*matrix_)[darcyIdx][darcyIdx],
                                                              (*matrix_)[darcyIdx][stokesIdx],
                                                              residual_[darcyIdx]);
    }

    template <class T = StokesProblemTypeTag>
    void buildStokesMatrixBlocks_(StokesMatrixBlock& m, StokesCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet stokesPattern, stokesCouplingPattern;
        stokesPattern.resize(m.N(), m.M());
        stokesCouplingPattern.resize(cm.N(), cm.M());

        const auto& assemblyMap = this->problem_().stokesProblem().model().localJacobian().assemblyMap();
        for (const auto& element : elements(problem_().stokesGridView()))
        {
            const auto globalI = problem_().model().stokesElementMapper().index(element);

            stokesPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap[globalI])
                stokesPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            const auto& additionalDofDependencies = problem_().stokesProblem().getAdditionalDofDependencies(globalI);
            for (auto globalJ : additionalDofDependencies)
                stokesPattern.add(globalI, globalJ);

            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            for (auto globalJ : couplingStencil)
                stokesCouplingPattern.add(globalI, globalJ);
        }

        // export occupation patterns to matrices
        stokesPattern.exportIdx(m);
        stokesCouplingPattern.exportIdx(cm);
    }

    template <class T = DarcyProblemTypeTag>
    void buildDarcyMatrixBlocks_(DarcyMatrixBlock& m, DarcyCouplingMatrixBlock& cm)
    {
        // get occupation pattern of the matrix
        Dune::MatrixIndexSet darcyPattern, darcyCouplingPattern;
        darcyPattern.resize(m.N(), m.M());
        darcyCouplingPattern.resize(cm.N(), cm.M());

        const auto& assemblyMap = this->problem_().darcyProblem().model().localJacobian().assemblyMap();
        for (const auto& element : elements(problem_().darcyGridView()))
        {
            const auto globalI = problem_().model().darcyElementMapper().index(element);

            darcyPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap[globalI])
                darcyPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            const auto& additionalDofDependencies = problem_().darcyProblem().getAdditionalDofDependencies(globalI);
            for (auto globalJ : additionalDofDependencies)
                darcyPattern.add(globalI, globalJ);

            const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);
            for (auto globalJ : couplingStencil)
                darcyCouplingPattern.add(globalI, globalJ);
        }

        // export occupation patterns to matrices
        darcyPattern.exportIdx(m);
        darcyCouplingPattern.exportIdx(cm);
    }

    // Construct the multitype matrix for the global jacobian
    void createMatrix_()
    {
        // create the multitype matrix
        matrix_ = std::make_shared<JacobianMatrix>();

        // get sub matrix sizes
        const auto stokesSize = problem_().model().stokesNumDofs();
        const auto darcySize = problem_().model().darcyNumDofs();

        // allocate the sub matrices
        auto A11 = StokesMatrixBlock(stokesSize, stokesSize, StokesMatrixBlock::random);
        auto A12 = StokesCouplingMatrixBlock(stokesSize, darcySize, StokesCouplingMatrixBlock::random);
        auto A22 = DarcyMatrixBlock(darcySize, darcySize, DarcyMatrixBlock::random);
        auto A21 = DarcyCouplingMatrixBlock(darcySize, stokesSize, DarcyCouplingMatrixBlock::random);

        buildStokesMatrixBlocks_(A11, A12);
        buildDarcyMatrixBlocks_(A22, A21);

        (*matrix_)[stokesIdx][stokesIdx] = A11;
        (*matrix_)[stokesIdx][darcyIdx] = A12;
        (*matrix_)[darcyIdx][stokesIdx] = A21;
        (*matrix_)[darcyIdx][darcyIdx] = A22;
    }

    // assemble a stokes element
    // void assembleElement_(const StokesElement &element)
    // {
    //     problem_().model().stokesLocalJacobian().assemble(element);

    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         for (unsigned int scvIdx = 0; scvIdx < element.subEntities(dim); ++scvIdx)
    //         {
    //             auto globalI = problem_().model().stokesVertexMapper().subIndex(element, scvIdx, dim);

    //             // update the right hand side
    //             residual_[stokesIdx][globalI] += problem_().model().stokesLocalJacobian().residual(scvIdx);
    //             for (int j = 0; j < residual_[stokesIdx][globalI].dimension; ++j)
    //                 assert(std::isfinite(residual_[stokesIdx][globalI][j]));

    //             int j = 0;
    //             for (auto&& globalJ : stencil)
    //                 (*matrix_)[stokesIdx][stokesIdx][globalI][globalJ] += problem_().model().stokesLocalJacobian().mat(stokesIdx, scvIdx, j++);

    //             j = 0;
    //             for (auto&& globalJ : couplingStencil)
    //                 (*matrix_)[stokesIdx][darcyIdx][globalI][globalJ] += problem_().model().stokesLocalJacobian().mat(darcyIdx, scvIdx, j++);
    //         }
    // }

    // // assemble a stokes element
    // void assembleElement_(const DarcyElement &element)
    // {
    //     problem_().model().darcyLocalJacobian().assemble(element);

    //         const auto& stencil = problem_().couplingManager().stencil(element);
    //         const auto& couplingStencil = problem_().couplingManager().couplingStencil(element);

    //         for (unsigned int scvIdx = 0; scvIdx < element.subEntities(dim); ++scvIdx)
    //         {
    //             auto globalI = problem_().model().darcyVertexMapper().subIndex(element, scvIdx, dim);

    //             // update the right hand side
    //             residual_[darcyIdx][globalI] +=  problem_().model().darcyLocalJacobian().residual(scvIdx);
    //             for (int j = 0; j < residual_[darcyIdx][globalI].dimension; ++j)
    //                 assert(std::isfinite(residual_[darcyIdx][globalI][j]));

    //             int j = 0;
    //             for (auto&& globalJ : stencil)
    //                 (*matrix_)[darcyIdx][darcyIdx][globalI][globalJ] += problem_().model().darcyLocalJacobian().mat(darcyIdx, scvIdx, j++);

    //             j = 0;
    //             for (auto&& globalJ : couplingStencil)
    //                 (*matrix_)[darcyIdx][stokesIdx][globalI][globalJ] += problem_().model().darcyLocalJacobian().mat(stokesIdx, scvIdx, j++);
    //         }
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
