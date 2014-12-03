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
#ifndef DUMUX_IMPLICIT_ASSEMBLER_HH
#define DUMUX_IMPLICIT_ASSEMBLER_HH

#include "implicitproperties.hh"

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class ImplicitAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum{ dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    // copying the jacobian assembler is not a good idea
    ImplicitAssembler(const ImplicitAssembler &);

public:
    /*!
     * \brief The colors of elements and vertices required for partial
     *        Jacobian reassembly.
     */
    enum EntityColor {
        /*!
         * Vertex/element that needs to be reassembled because some
         * relative error is above the tolerance
         */
        Red = 0,

        /*!
         * Vertex/element that needs to be reassembled because a
         * neighboring element/vertex is red
         */
        Yellow = 1,

        /*!
         * Yellow vertex has only non-green neighbor elements.
         *
         * This means that its relative error is below the tolerance,
         * but its defect can be linearized without any additional
         * cost. This is just an "internal" color which is not used
         * ouside of the jacobian assembler.
         */
        Orange = 2,

        /*!
         * Vertex/element that does not need to be reassembled
         */
        Green = 3
    };

    ImplicitAssembler()
    : problemPtr_(0)
    {
        // set reassemble accuracy to 0, so that if partial reassembly
        // of the jacobian matrix is disabled, the reassemble accuracy
        // is always smaller than the current relative tolerance
        reassembleAccuracy_ = 0.0;
    }

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
        problemPtr_ = &problem;

        // initialize the BCRS matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseMatrix_ = false;

        int numVertices = gridView_().size(dim);
        int numElements = gridView_().size(0);
        int numDofs = problem.model().numDofs();

        residual_.resize(numDofs);

        // initialize the storage part of the Jacobian matrix. Since
        // we only need this if Jacobian matrix recycling is enabled,
        // we do not waste space if it is disabled
        if (enableJacobianRecycling_()) {
            storageJacobian_.resize(numDofs);
            storageTerm_.resize(numDofs);
        }

        if (gridView_().comm().size() > 1)
            totalElems_ = gridView_().comm().sum(numElements);
        else
            totalElems_ = numElements;

        // initialize data needed for partial reassembly
        if (enablePartialReassemble_()) {
            delta_.resize(numDofs);
            elementColor_.resize(numElements);
            if (isBox)
                vertexColor_.resize(numVertices);
        }
        reassembleAll();
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool printReassembleStatistics = enablePartialReassemble_() && !reuseMatrix_;
        int succeeded;
        try {
            asImp_().assemble_();
            succeeded = 1;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }
        catch (Dumux::NumericalProblem &e)
        {
            std::cout << "rank " << problem_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = 0;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            DUNE_THROW(NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }

        if (printReassembleStatistics)
        {
            if (gridView_().comm().size() > 1)
            {
                greenElems_ = gridView_().comm().sum(greenElems_);
                reassembleAccuracy_ = gridView_().comm().max(nextReassembleAccuracy_);
            }
            else
            {
                reassembleAccuracy_ = nextReassembleAccuracy_;
            }

            problem_().newtonController().endIterMsg()
                << ", reassembled "
                << totalElems_ - greenElems_ << "/" << totalElems_
                << " (" << 100*Scalar(totalElems_ - greenElems_)/totalElems_ << "%) elems @accuracy="
                << reassembleAccuracy_;
        }

        // reset all vertex colors to green
        for (unsigned int i = 0; i < vertexColor_.size(); ++i) {
            vertexColor_[i] = Green;
        }
    }

    /*!
     * \brief If Jacobian matrix recycling is enabled, this method
     *        specifies whether the next call to assemble() just
     *        rescales the storage term or does a full reassembly
     *
     * \param yesno If true, only rescale; else do full Jacobian assembly.
     */
    void setMatrixReuseable(const bool yesno = true)
    {
        if (enableJacobianRecycling_())
            reuseMatrix_ = yesno;
    }

    /*!
     * \brief If partial Jacobian matrix reassembly is enabled, this
     *        method causes all elements to be reassembled in the next
     *        assemble() call.
     */
    void reassembleAll()
    {
        // do not reuse the current linearization
        reuseMatrix_ = false;

        // do not use partial reassembly for the next iteration
        nextReassembleAccuracy_ = 0.0;
        if (enablePartialReassemble_()) {
            std::fill(vertexColor_.begin(),
                      vertexColor_.end(),
                      Red);
            std::fill(elementColor_.begin(),
                      elementColor_.end(),
                      Red);
            std::fill(delta_.begin(),
                      delta_.end(),
                      0.0);
        }
    }

    /*!
     * \brief Returns the largest relative error of a "green" vertex
     *        for the most recent call of the assemble() method.
     *
     * This only has an effect if partial Jacobian reassembly is
     * enabled. If it is disabled, then this method always returns 0.
     *
     * This returns the _actual_ relative computed seen by
     * computeColors(), not the tolerance which it was given.
     */
    Scalar reassembleAccuracy() const
    { return reassembleAccuracy_; }

    /*!
     * \brief Update the distance where the non-linear system was
     *        originally insistently linearized and the point where it
     *        will be linerized the next time.
     *
     * This only has an effect if partial reassemble is enabled.
     */
    void updateDiscrepancy(const SolutionVector &u,
                           const SolutionVector &uDelta)
    {
        if (!enablePartialReassemble_())
            return;

        // update the vector with the distances of the current
        // evaluation point used for linearization from the original
        // evaluation point
        for (unsigned int i = 0; i < delta_.size(); ++i) {
            PrimaryVariables currentPriVars(u[i]);
            PrimaryVariables nextPriVars(currentPriVars);
            nextPriVars -= uDelta[i];

            // we need to add the distance the solution was moved for
            // this vertex
            Scalar dist = model_().relativeErrorDof(i,
                                                    currentPriVars,
                                                    nextPriVars);
            delta_[i] += std::abs(dist);
        }

    }

    /*!
     * \brief Force to reassemble a given degree of freedom
     * next time the assemble() method is called.
     *
     * \param dofIdxGlobal The global index of the degree of freedom
     */
    void markDofRed(const int dofIdxGlobal)
    {
        if (!enablePartialReassemble_())
            return;

        if (isBox)
            vertexColor_[dofIdxGlobal] = Red;
        else 
            elementColor_[dofIdxGlobal] = Red;
    }

    /*!
     * \brief Determine the colors of vertices and elements for partial
     *        reassembly given a relative tolerance.
     *
     * \param relTol The relative error below which a vertex won't be
     *               reassembled. Note that this specifies the
     *               worst-case relative error between the last
     *               linearization point and the current solution and
     *               _not_ the delta vector of the Newton iteration!
     */
    void computeColors(const Scalar relTol)
    {
        asImp_().computeColors_(relTol);
    }

    /*!
     * \brief Returns the reassemble color of a vertex
     *
     * \param element An element which contains the vertex
     * \param vIdx The local index of the vertex in the element.
     */
    int vertexColor(const Element &element, const int vIdx) const
    {
        if (!enablePartialReassemble_())
            return Red; // reassemble unconditionally!

        int vIdxGlobal = vertexMapper_().map(element, vIdx, dim);
        return vertexColor_[vIdxGlobal];
    }

    /*!
     * \brief Returns the reassemble color of a vertex
     *
     * \param vIdxGlobal The global index of the vertex.
     */
    int vertexColor(const int vIdxGlobal) const
    {
        if (!enablePartialReassemble_())
            return Red; // reassemble unconditionally!
        return vertexColor_[vIdxGlobal];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \param element The Codim-0 DUNE entity
     */
    int elementColor(const Element &element) const
    {
        if (!enablePartialReassemble_())
            return Red; // reassemble unconditionally!

        int eIdxGlobal = elementMapper_().map(element);
        return elementColor_[eIdxGlobal];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \param globalElementIdx The global index of the element.
     */
    int elementColor(const int globalElementIdx) const
    {
        if (!enablePartialReassemble_())
            return Red; // reassemble unconditionally!
        return elementColor_[globalElementIdx];
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
    static bool enableJacobianRecycling_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableJacobianRecycling); }
    static bool enablePartialReassemble_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnablePartialReassemble); }


    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // do not do anything if we can re-use the current linearization
        if (reuseMatrix_)
            return;

        // reset the right hand side.
        residual_ = 0.0;

        if (!enablePartialReassemble_()) {
            // If partial reassembly of the jacobian is not enabled,
            // we can just reset everything!
            (*matrix_) = 0;

            // reset the parts needed for Jacobian recycling
            if (enableJacobianRecycling_()) {
                for (unsigned int i = 0; i < matrix_->N(); i++) {
                    storageJacobian_[i] = 0;
                    storageTerm_[i] = 0;
                }
            }

            return;
        }

        // reset all entries corrosponding to a red or yellow vertex
        for (unsigned int rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if ((isBox && vertexColor_[rowIdx] != Green)
                || (!isBox && elementColor_[rowIdx] != Green))
            {
                // reset the parts needed for Jacobian recycling
                if (enableJacobianRecycling_()) {
                    storageJacobian_[rowIdx] = 0;
                    storageTerm_[rowIdx] = 0;
                }

                // set all matrix entries in the row to 0
                typedef typename JacobianMatrix::ColIterator ColIterator;
                ColIterator colIt = (*matrix_)[rowIdx].begin();
                const ColIterator &colEndIt = (*matrix_)[rowIdx].end();
                for (; colIt != colEndIt; ++colIt) {
                    (*colIt) = 0.0;
                }
            }
        }
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().timeManager().timeStepSize();
        if (reuseMatrix_) {
            for (unsigned int i = 0; i < matrix_->N(); i++) {
                // rescale the mass term of the jacobian matrix
                MatrixBlock &J_i_i = (*matrix_)[i][i];

                J_i_i -= storageJacobian_[i];
                storageJacobian_[i] *= oldDt_/curDt;
                J_i_i += storageJacobian_[i];

                // use the flux term plus the source term as the new
                // residual (since the delta in the d(storage)/dt is 0
                // for the first iteration and the residual is
                // approximately 0 in the last iteration, the flux
                // term plus the source term must be equal to the
                // negative change of the storage term of the last
                // iteration of the last time step...)
                residual_[i] = storageTerm_[i];
                residual_[i] *= -1.0;
            }

            reuseMatrix_ = false;
            oldDt_ = curDt;
            return;
        }

        oldDt_ = curDt;
        greenElems_ = 0;

        // reassemble the elements...
        ElementIterator eIt = gridView_().template begin<0>();
        ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            const Element &elem = *eIt;
            if (elem.partitionType() == Dune::GhostEntity)
            {
                asImp_().assembleGhostElement_(elem);
            }
            else
            {
                asImp_().assembleElement_(elem);
            }
        }
    }

protected:
    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return problem_().model(); }
    Model &model_()
    { return problem_().model(); }
    const GridView &gridView_() const
    { return problem_().gridView(); }
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); }
    const ElementMapper &elementMapper_() const
    { return problem_().elementMapper(); }

    Problem *problemPtr_;

    // the jacobian matrix
    Dune::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    SolutionVector residual_;

    // attributes required for jacobian matrix recycling
    bool reuseMatrix_;
    // The storage part of the local Jacobian
    std::vector<MatrixBlock> storageJacobian_;
    std::vector<VectorBlock> storageTerm_;
    // time step size of last assembly
    Scalar oldDt_;


    // attributes required for partial jacobian reassembly
    std::vector<EntityColor> vertexColor_;
    std::vector<EntityColor> elementColor_;
    std::vector<Scalar> delta_;

    int totalElems_;
    int greenElems_;

    Scalar nextReassembleAccuracy_;
    Scalar reassembleAccuracy_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
