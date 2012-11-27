// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief An assembler for the global Jacobian matrix for models using the box discretization.
 */
#ifndef DUMUX_CC_ASSEMBLER_HH
#define DUMUX_CC_ASSEMBLER_HH

#include <dune/grid/common/gridenums.hh>

#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/foreignoverlapfrombcrsmatrix.hh>
#include <dumux/parallel/vertexhandles.hh>

namespace Dumux {

/*!
 * \brief An assembler for the global Jacobian matrix for models using the box discretization.
 */
template<class TypeTag>
class CCAssembler
{
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
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    // copying the jacobian assembler is not a good idea
    CCAssembler(const CCAssembler &);

public:
    /*!
     * \brief The colors of elements required for partial
     *        Jacobian reassembly.
     */
    enum EntityColor {
        /*!
         * Element that needs to be reassembled because some
         * relative error is above the tolerance
         */
        Red = 0,

        /*!
         * Element that does not need to be reassembled
         */
        Green = 3
    };

    CCAssembler()
    {
        problemPtr_ = 0;
        matrix_ = 0;

        // set reassemble accuracy to 0, so that if partial reassembly
        // of the jacobian matrix is disabled, the reassemble accuracy
        // is always smaller than the current relative tolerance
        reassembleAccuracy_ = 0.0;
    }

    ~CCAssembler()
    {
        delete matrix_;
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
        createMatrix_();

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseMatrix_ = false;

        int numElems = gridView_().size(0);

        residual_.resize(numElems);

        // initialize the storage part of the Jacobian matrix. Since
        // we only need this if Jacobian matrix recycling is enabled,
        // we do not waste space if it is disabled
        if (enableJacobianRecycling_()) {
            storageJacobian_.resize(numElems);
            storageTerm_.resize(numElems);
        }

        if (gridView_().comm().size() > 1)
            totalElems_ = gridView_().comm().sum(numElems);
        else
            totalElems_ = numElems;

        // initialize data needed for partial reassembly
        if (enablePartialReassemble_()) 
        {
            elementColor_.resize(numElems);
            elementDelta_.resize(numElems);
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
            assemble_();
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
    }

   /*!
    * \brief If Jacobian matrix recycling is enabled, this method
    *        specifies whether the next call to assemble() just
    *        rescales the storage term or does a full reassembly
    *
    * \param yesno If true, only rescale; else do full Jacobian assembly.
    */
   void setMatrixReuseable(bool yesno = true) 
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
            std::fill(elementColor_.begin(),
                      elementColor_.end(),
                      Red);
            std::fill(elementDelta_.begin(),
                      elementDelta_.end(),
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
        for (unsigned int i = 0; i < elementDelta_.size(); ++i) {
            PrimaryVariables currentPriVars(u[i]);
            PrimaryVariables nextPriVars(currentPriVars);
            nextPriVars -= uDelta[i];

            // we need to add the distance the solution was moved for
            // this vertex
            Scalar dist = model_().relativeErrorVertex(i,
                                                       currentPriVars,
                                                       nextPriVars);
            elementDelta_[i] += std::abs(dist);
        }
    }
    
    /*!
     * \brief Determine the colors of elements for partial
     *        reassembly given a relative tolerance.
     *
     * The following approach is used:
     *
     * - Set all elements to 'green'
     * - Mark all elements as 'red' which exhibit an relative error above
     *   the tolerance
     * - Mark all neighbors of 'red' elements also 'red'
     *
     * \param relTol The relative error below which an element won't be
     *               reassembled. Note that this specifies the
     *               worst-case relative error between the last
     *               linearization point and the current solution and
     *               _not_ the delta vector of the Newton iteration!
     */
    void computeColors(Scalar relTol) 
    {
        if (!enablePartialReassemble_())
            return;

        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();

        // mark the red elements and update the tolerance of the
        // linearization which actually will get achieved
        nextReassembleAccuracy_ = 0;
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementDelta_[elemIdx] > relTol)
            {
                // mark element as red if discrepancy is larger than
                // the relative tolerance
                elementColor_[elemIdx] = Red;
            }
            else
            {
                elementColor_[elemIdx] = Green;
                nextReassembleAccuracy_ =
                    std::max(nextReassembleAccuracy_, elementDelta_[elemIdx]);
            }
        }

        // mark the neighbors also red
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Red) 
                continue; // element is red already!

            if (elementDelta_[elemIdx] > relTol)
            {
                // also mark the neighbors 
               IntersectionIterator endIsIt = gridView_().iend(*elemIt);
               for (IntersectionIterator isIt = gridView_().ibegin(*elemIt); isIt != endIsIt; ++isIt)
               {
                   if (isIt->neighbor())
                   {
                       int neighborIdx = this->elementMapper_().map(*isIt->outside());
                       elementColor_[neighborIdx] = Red;
                   }
               }
            }
        }
        
        // set the discrepancy of the red elements to zero
        for (int i = 0; i < elementDelta_.size(); i++)
            if (elementColor_[i] == Red)
                elementDelta_[i] = 0;
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

        int globalIdx = elementMapper_().map(element);
        return elementColor_[globalIdx];
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
    
    // functions needed for interface consistence
    void markVertexRed(const int globalVertIdx) {}

private:
    static bool enableJacobianRecycling_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableJacobianRecycling); }
    static bool enablePartialReassemble_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnablePartialReassemble); }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int nElems = gridView_().size(0);

        // allocate raw matrix
        matrix_ = new JacobianMatrix(nElems, nElems, JacobianMatrix::random);

        // find out the global indices of the neighboring elements of
        // each element
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(nElems);
        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            const Element &elem = *eIt;
            
            int globalI = elementMapper_().map(elem);
            neighbors[globalI].insert(globalI);

            // if the element is ghost,
            // all dofs just contain main-diagonal entries
            //if (elem.partitionType() == Dune::GhostEntity)
            //    continue;

            // loop over all neighbors
            IntersectionIterator isIt = gridView_().ibegin(elem);
            const IntersectionIterator &endIt = gridView_().iend(elem);
            for (; isIt != endIt; ++isIt)
            {
                if (isIt->neighbor())
                {
                    int globalJ = elementMapper_().map(*(isIt->outside()));
                    neighbors[globalI].insert(globalJ);
                }
            }
        }

        // allocate space for the rows of the matrix
        for (int i = 0; i < nElems; ++i) {
            matrix_->setrowsize(i, neighbors[i].size());
        }
        matrix_->endrowsizes();

        // fill the rows with indices. each element talks to all of its
        // neighbors and itself.
        for (int i = 0; i < nElems; ++i) {
            typename NeighborSet::iterator nIt = neighbors[i].begin();
            typename NeighborSet::iterator nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                matrix_->addindex(i, *nIt);
            }
        }
        matrix_->endindices();
    }
    
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
                int numElementsGlobal = matrix_->N();
                for (int i=0; i < numElementsGlobal; ++ i) {
                    storageJacobian_[i] = 0;
                    storageTerm_[i] = 0;
                }
            }

            return;
        }

        // reset all entries corrosponding to a red or yellow vertex
        for (unsigned int rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if (elementColor_[rowIdx] == Green)
                continue; // the equations for this control volume are
                          // already below the treshold

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

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().timeManager().timeStepSize();
        if (reuseMatrix_) {
            int numElementsGlobal = storageJacobian_.size();
            for (int i = 0; i < numElementsGlobal; ++i) {
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
                residual_[i] *= -1;
            }

            reuseMatrix_ = false;
            oldDt_ = curDt;
            return;
        }

        oldDt_ = curDt;
        greenElems_ = 0;
        
        // reassemble the elements...
        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element &elem = *elemIt;
            if (elem.partitionType() == Dune::GhostEntity)
            {
                assembleGhostElement_(elem);
            }
            else
            {
                assembleElement_(elem);
            }
        }
    }

    // assemble a non-ghost element
    void assembleElement_(const Element &elem)
    {
        if (enablePartialReassemble_()) {
            int globalElemIdx = model_().elementMapper().map(elem);
            if (elementColor_[globalElemIdx] == Green) {
                ++greenElems_;

                assembleGreenElement_(elem);
                return;
            }
        }

        model_().localJacobian().assemble(elem);

        int globalI = elementMapper_().map(elem);

        // update the right hand side
        residual_[globalI] = model_().localJacobian().residual(0);
        for (int j = 0; j < residual_[globalI].dimension; ++j)
            assert(std::isfinite(residual_[globalI][j]));
        if (enableJacobianRecycling_()) {
            storageTerm_[globalI] +=
                model_().localJacobian().storageTerm(0);
        }

        if (enableJacobianRecycling_())
            storageJacobian_[globalI] +=
                model_().localJacobian().storageJacobian(0);

        // update the diagonal entry 
        (*matrix_)[globalI][globalI] = model_().localJacobian().mat(0,0);

        IntersectionIterator isIt = gridView_().ibegin(elem);
        const IntersectionIterator &endIt = gridView_().iend(elem);
        for (int j = 0; isIt != endIt; ++isIt)
        {
            if (isIt->neighbor())
            {
                int globalJ = elementMapper_().map(*(isIt->outside()));
                (*matrix_)[globalI][globalJ] = model_().localJacobian().mat(0,++j);
            }
        }
    }

    // "assemble" a green element. green elements only get the
    // residual updated, but the jacobian is left alone...
    void assembleGreenElement_(const Element &elem)
    {
        model_().localResidual().eval(elem);

        int globalI = elementMapper_().map(elem);

        // update the right hand side
        residual_[globalI] += model_().localResidual().residual(0);
        if (enableJacobianRecycling_())
            storageTerm_[globalI] += model_().localResidual().storageTerm(0);
    }

    // "assemble" a ghost element
    void assembleGhostElement_(const Element &elem)
    {
        int globalI = elementMapper_().map(elem);

        // update the right hand side
        residual_[globalI] = 0.0;

        // update the diagonal entry 
        typedef typename JacobianMatrix::block_type BlockType;
        BlockType &J = (*matrix_)[globalI][globalI];
        for (int j = 0; j < BlockType::rows; ++j)
            J[j][j] = 1.0;
    }


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
    JacobianMatrix *matrix_;
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
    std::vector<EntityColor> elementColor_;
    std::vector<Scalar> elementDelta_;

    int totalElems_;
    int greenElems_;

    Scalar nextReassembleAccuracy_;
    Scalar reassembleAccuracy_;
};

} // namespace Dumux

#endif
