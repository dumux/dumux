// $Id$
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief An assembler for the Jacobian matrix based on PDELab.
 */
#ifndef DUMUX_BOX_ASSEMBLER_HH
#define DUMUX_BOX_ASSEMBLER_HH

#if HAVE_DUNE_PDELAB
#include "pdelabboxlocaloperator.hh"
#endif

//#include "overlapmatrix.hh"
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/foreignoverlapfrombcrsmatrix.hh>

#include <dumux/parallel/vertexhandles.hh>

namespace Dumux {

/*!
 * \brief An assembler for the Jacobian matrix based on PDELab.
 */
template<class TypeTag>
class BoxAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

#if HAVE_DUNE_PDELAB
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Constraints)) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ScalarGridFunctionSpace)) ScalarGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalOperator)) LocalOperator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridOperatorSpace)) GridOperatorSpace;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    enum{dim = GridView::dimension};
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef SolutionVector Vector;
    typedef JacobianMatrix Matrix;
    typedef Matrix RepresentationType;

    enum {
        enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)),
        enableJacobianRecycling = GET_PROP_VALUE(TypeTag, PTAG(EnableJacobianRecycling)),

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))
    };

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    // copying the jacobian assembler is not a good idea
    BoxAssembler(const BoxAssembler &);

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

    BoxAssembler()
    {
#if HAVE_DUNE_PDELAB
        fem_ = 0;
        cn_ = 0;
        scalarGridFunctionSpace_ = 0;
        gridFunctionSpace_ = 0;
        constraintsTrafo_ = 0;
        localOperator_ = 0;
        gridOperatorSpace_ = 0;
#endif // HAVE_DUNE_PDELAB

        problemPtr_ = 0;
        matrix_ = 0;

        // set reassemble accuracy to 0, so that if partial reassembly
        // of the jacobian matrix is disabled, the reassemble accuracy
        // is always smaller than the current relative tolerance
        reassembleAccuracy_ = 0.0;
    }

    ~BoxAssembler()
    {
        delete matrix_;

#if HAVE_DUNE_PDELAB
        delete gridOperatorSpace_;
        delete localOperator_;
        delete constraintsTrafo_;
        delete gridFunctionSpace_;
        delete scalarGridFunctionSpace_;
        delete cn_;
        delete fem_;
#endif
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

#if !HAVE_DUNE_PDELAB
        // initialize the BCRS matrix
        createMatrix_();
#else
        fem_ = new FEM();
        //cn_ = new Constraints(*problemPtr_);
        cn_ = new Constraints();
        scalarGridFunctionSpace_ = new ScalarGridFunctionSpace(problemPtr_->gridView(), *fem_, *cn_);
        gridFunctionSpace_ = new GridFunctionSpace(*scalarGridFunctionSpace_);
        //cn_->compute_ghosts(*gridFunctionSpace_);

        //typedef BoundaryIndexHelper<TypeTag> BoundaryFunction;
        //BoundaryFunction *bTypes = new BoundaryFunction();
        constraintsTrafo_ = new ConstraintsTrafo();
        //Dune::PDELab::constraints(*bTypes, *gridFunctionSpace_, *constraintsTrafo_, false);

        // initialize the grid operator spaces
        localOperator_ = new LocalOperator(problemPtr_->model());
        gridOperatorSpace_ =
            new GridOperatorSpace(*gridFunctionSpace_, *constraintsTrafo_,
                                  *gridFunctionSpace_, *constraintsTrafo_, *localOperator_);
        matrix_ = new Matrix(*gridOperatorSpace_);
#endif

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseMatrix_ = false;

        int numVerts = gridView_().size(dim);
        int numElems = gridView_().size(0);

        residual_.resize(numVerts);

        // initialize the storage part of the Jacobian matrix. Since
        // we only need this if Jacobian matrix recycling is enabled,
        // we do not waste space if it is disabled
        if (enableJacobianRecycling) {
            storageJacobian_.resize(numVerts);
            storageTerm_.resize(numVerts);
        }

        totalElems_ = gridView_().comm().sum(numElems);

        // initialize data needed for partial reassembly
        if (enablePartialReassemble) {
            vertexColor_.resize(numVerts);
            vertexDelta_.resize(numVerts);
            elementColor_.resize(numElems);
        }
        reassembleAll();
    }

    /*!
     * \brief Assemble the local jacobian of the problem.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        resetSystem_();
        
        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().timeManager().timeStepSize();
        if (reuseMatrix_) {
            int numVertices = storageJacobian_.size();
            for (int i = 0; i < numVertices; ++i) {
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
            };
            
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
            if (elem.partitionType() != Dune::InteriorEntity  &&
                elem.partitionType() != Dune::BorderEntity)
            {
                assembleGhostElement_(elem);
            }
            else
            {
                assembleElement_(elem);
            }
        };

        // if partial reassembly is enabled, print some statistics at
        // the end of the iteration
        if (enablePartialReassemble) {
            greenElems_ = gridView_().comm().sum(greenElems_);
            reassembleAccuracy_ = gridView_().comm().max(nextReassembleAccuracy_);

            problem_().newtonController().endIterMsg()
                << ", reassembled "
                << totalElems_ - greenElems_ << "/" << totalElems_
                << " (" << 100*Scalar(totalElems_ - greenElems_)/totalElems_ << "%) elems @accuracy="
                << reassembleAccuracy_;
        }

        return;
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
        if (enableJacobianRecycling)
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
        if (enablePartialReassemble) {
            std::fill(vertexColor_.begin(),
                      vertexColor_.end(),
                      Red);
            std::fill(elementColor_.begin(),
                      elementColor_.end(),
                      Red);
            std::fill(vertexDelta_.begin(),
                      vertexDelta_.end(),
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
        if (!enablePartialReassemble)
            return;

        // update the vector with the distances of the current
        // evaluation point used for linearization from the original
        // evaluation point
        for (int i = 0; i < vertexDelta_.size(); ++i) {
            PrimaryVariables uCurrent(u[i]);
            PrimaryVariables uNext(uCurrent);
            uNext -= uDelta[i];

            // we need to add the distance the solution was moved for
            // this vertex
            Scalar dist = model_().relativeErrorVertex(i,
                                                       uCurrent,
                                                       uNext);
            vertexDelta_[i] += std::abs(dist);
        }

    }

    /*!
     * \brief Determine the colors of vertices and elements for partial
     *        reassembly given a relative tolerance.
     *
     * The following approach is used:
     *
     * - Set all vertices and elements to 'green'
     * - Mark all vertices as 'red' which exhibit an relative error above
     *   the tolerance
     * - Mark all elements which feature a 'red' vetex as 'red'
     * - Mark all vertices which are not 'red' and are part of a
     *   'red' element as 'yellow'
     * - Mark all elements which are not 'red' and contain a
     *   'yellow' vertex as 'yellow'
     *
     * \param relTol The relative error below which a vertex won't be
     *               reassembled. Note that this specifies the
     *               worst-case relative error between the last
     *               linearization point and the current solution and
     *               _not_ the delta vector of the Newton iteration!
     */
    void computeColors(Scalar relTol)
    {
        if (!enablePartialReassemble)
            return;

        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();

        // mark the red vertices and update the tolerance of the
        // linearization which actually will get achieved
        nextReassembleAccuracy_ = 0;
        for (int i = 0; i < vertexColor_.size(); ++i) {
            vertexColor_[i] = Green;
            if (vertexDelta_[i] > relTol)
                // mark vertex as red if discrepancy is larger than
                // the relative tolerance
                vertexColor_[i] = Red;
            else
                nextReassembleAccuracy_ = 
                    std::max(nextReassembleAccuracy_, vertexDelta_[i]);
        };

        // Mark all red elements
        for (; elemIt != elemEndIt; ++elemIt) {
            // find out whether the current element features a red
            // vertex
            bool isRed = false;
            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                if (vertexColor_[globalI] == Red) {
                    isRed = true;
                    break;
                }
            };

            // if yes, the element color is also red, else it is not
            // red, i.e. green for the mean time
            int globalElemIdx = elementMapper_().map(*elemIt);
            if (isRed)
                elementColor_[globalElemIdx] = Red;
            else
                elementColor_[globalElemIdx] = Green;
        }

        // Mark yellow vertices (as orange for the mean time)
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Red)
                continue; // non-red elements do not tint vertices
                          // yellow!

            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                // if a vertex is already red, don't recolor it to
                // yellow!
                if (vertexColor_[globalI] != Red) {
                    vertexColor_[globalI] = Orange;
                }
            };
        }

        // at this point we communicate the yellow vertices to the
        // neighboring processes because a neigbor process may not see
        // the red vertex for yellow border vertices
        VertexHandleMin<EntityColor, std::vector<EntityColor>,  VertexMapper>
            minHandle(vertexColor_, vertexMapper_());
        gridView_().communicate(minHandle, 
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);
        
        // Mark yellow elements
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Red) {
                continue; // element is red already!
            }

            // check whether the element features a yellow
            // (resp. orange at this point) vertex
            bool isYellow = false;
            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                if (vertexColor_[globalI] == Orange) {
                    isYellow = true;
                    break;
                }
            };

            if (isYellow)
                elementColor_[elemIdx] = Yellow;
        }

        // Demote orange vertices to yellow ones if it has at least
        // one green element as a neighbor.
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Green)
                continue; // yellow and red elements do not make
                          // orange vertices yellow!

            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                // if a vertex is orange, recolor it to yellow!
                if (vertexColor_[globalI] == Orange)
                    vertexColor_[globalI] = Yellow;
            };
        }

        // demote the border orange vertices
        VertexHandleMax<EntityColor, std::vector<EntityColor>,  VertexMapper>
            maxHandle(vertexColor_,
                      vertexMapper_());
        gridView_().communicate(maxHandle, 
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // promote the remaining orange vertices to red
        for (int i=0; i < vertexColor_.size(); ++i) {
            // if a vertex is green or yellow don't do anything!
            if (vertexColor_[i] == Green || vertexColor_[i] == Yellow)
                continue;

            // make sure the vertex is red (this is a no-op vertices
            // which are already red!)
            vertexColor_[i] = Red;

            // set the error of this vertex to 0 because the system
            // will be consistently linearized at this vertex
            vertexDelta_[i] = 0.0;
        };
    };

    /*!
     * \brief Returns the reassemble color of a vertex
     *
     * \param element An element which contains the vertex
     * \param vertIdx The local index of the vertex in the element.
     */
    int vertexColor(const Element &element, int vertIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!

        int globalIdx = vertexMapper_().map(element, vertIdx, dim);
        return vertexColor_[globalIdx];
    }

    /*!
     * \brief Returns the reassemble color of a vertex
     *
     * \param globalVertIdx The global index of the vertex.
     */
    int vertexColor(int globalVertIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        return vertexColor_[globalVertIdx];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \param element The Codim-0 DUNE entity
     */
    int elementColor(const Element &element) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!

        int globalIdx = elementMapper_().map(element);
        return elementColor_[globalIdx];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \param globalElementIdx The global index of the element.
     */
    int elementColor(int globalElementIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        return elementColor_[globalElementIdx];
    }

#if HAVE_DUNE_PDELAB
    /*!
     * \brief Returns a pointer to the PDELab's grid function space.
     */
    const GridFunctionSpace& gridFunctionSpace() const
    {
        return *gridFunctionSpace_;
    }

    /*!
     * \brief Returns a pointer to the PDELab's constraints
     *        transformation.
     */
    const ConstraintsTrafo& constraintsTrafo() const
    {
        return *constraintsTrafo_;
    }
#endif // HAVE_DUNE_PDELAB

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const Matrix& matrix() const
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return residual_; }


private:
#if !HAVE_DUNE_PDELAB
    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int nVerts = gridView_().size(dim);

        // allocate raw matrix
        matrix_ = new Matrix(nVerts, nVerts, Matrix::random);

        // find out the global indices of the neighboring vertices of
        // each vertex
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(nVerts);
        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            const Element &elem = *eIt;
            int n = elem.template count<dim>();
            
            // if the element is not in the interior or the process
            // border, all dofs just contain main-diagonal entries
            if (elem.partitionType() != Dune::InteriorEntity &&
                elem.partitionType() != Dune::BorderEntity) 
            {
                for (int i = 0; i < n; ++i) {
                    int globalI = vertexMapper_().map(*eIt, i, dim);
                    neighbors[globalI].insert(globalI);
                }
                continue;
            };

            // loop over all element vertices
            for (int i = 0; i < n - 1; ++i) {
                int globalI = vertexMapper_().map(*eIt, i, dim);
                for (int j = i + 1; j < n; ++j) {
                    int globalJ = vertexMapper_().map(*eIt, j, dim);
                    // make sure that vertex j is in the neighbor set
                    // of vertex i and vice-versa
                    neighbors[globalI].insert(globalJ);
                    neighbors[globalJ].insert(globalI);
                }
            }
        };

        // make vertices neighbors to themselfs
        for (int i = 0; i < nVerts; ++i)
            neighbors[i].insert(i);

        // allocate space for the rows of the matrix
        for (int i = 0; i < nVerts; ++i) {
            matrix_->setrowsize(i, neighbors[i].size());
        }
        matrix_->endrowsizes();

        // fill the rows with indices. each vertex talks to all of its
        // neighbors. (it also talks to itself since vertices are
        // sometimes quite egocentric.)
        for (int i = 0; i < nVerts; ++i) {
            typename NeighborSet::iterator nIt = neighbors[i].begin();
            typename NeighborSet::iterator nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                matrix_->addindex(i, *nIt);
            }
        }
        matrix_->endindices();

#if 0        
        typedef typename Dumux::OverlapFromBCRSMatrix<Matrix> OverlapFromBCRSMatrix;
        typedef typename Dumux::VertexBorderListFromGrid<GridView, VertexMapper> BorderListFromGrid;
        typedef typename OverlapFromBCRSMatrix::BorderList BorderList;
        
        BorderListFromGrid borderListCreator(gridView_(), vertexMapper_());
        int overlapSize = 1;
        OverlapFromBCRSMatrix overlap(*matrix_, 
                                      borderListCreator.borderList(),
                                      overlapSize);
        
        /*
        if (gridView_().comm().rank() == 0) {
            std::cout << "initial border list size: " << borderListCreator.borderList().size() << "\n";
            overlap.printOverlap();
        }
        */
        
/*
        typedef OverlapFromBCRSMatrix Overlap;
        typedef OverlapBCRSMatrix<Matrix, Overlap> OverlapMatrix;
        OverlapMatrix overlapMatrix(*matrix_, overlap);
        overlapMatrix.assignFromNonOverlapping(*matrix_);

        typedef Dune::BlockVector<PrimaryVariables> BlockVector;
        typedef Dumux::OverlapBlockVector<BlockVector, Overlap> OverlapBlockVector;
        BlockVector bv(matrix_->N());
        bv = gridView_().comm().rank();
        OverlapBlockVector overlapVector(bv, overlap);
*/
#endif
    };
#endif

    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // do not do anything if we can re-use the current linearization
        if (reuseMatrix_)
            return;

        // always reset the right hand side.
        residual_ = 0.0;

        if (!enablePartialReassemble) {
            // If partial reassembly of the jacobian is not enabled,
            // we can just reset everything!
            (*matrix_) = 0;

            // reset the parts needed for Jacobian recycling
            if (enableJacobianRecycling) {
                int numVertices = matrix_->N();
                for (int i=0; i < numVertices; ++ i) {
                    storageJacobian_[i] = 0;
                    storageTerm_[i] = 0;
                };
            }
            
            return;
        }

        // reset all entries corrosponding to a red vertex
        for (int rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if (vertexColor_[rowIdx] == Green)
                continue; // the equations for this control volume are
                          // already below the treshold

            // reset the parts needed for Jacobian recycling
            if (enableJacobianRecycling) {
                storageJacobian_[rowIdx] = 0;
                storageTerm_[rowIdx] = 0;
            }

            // set all entries in the row to 0
            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[rowIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                (*colIt) = 0.0;
            }
        };
    }

    // assemble a non-ghost element
    void assembleElement_(const Element &elem)
    {
        if (enablePartialReassemble) {
            int globalElemIdx = model_().elementMapper().map(elem);
            if (elementColor_[globalElemIdx] == Green) {
                ++greenElems_;
                return;
            }
        }

        model_().localJacobian().assemble(elem);

        int numVertices = elem.template count<dim>();
        for (int i=0; i < numVertices; ++ i) {
            int globI = vertexMapper_().map(elem, i, dim);

            // update the right hand side
            if (vertexColor(globI) == Green) {
                continue;
            }

            residual_[globI] += model_().localJacobian().residual(i);
            if (enableJacobianRecycling) {
                // save the flux term and the jacobian of the
                // storage term in case we can reuse the current
                // linearization...
                storageJacobian_[globI] +=
                    model_().localJacobian().storageJacobian(i);

                storageTerm_[globI] +=
                    model_().localJacobian().storageTerm(i);
            }

            // update the jacobian matrix
            for (int j=0; j < numVertices; ++ j) {
                int globJ = vertexMapper_().map(elem, j, dim);
                (*matrix_)[globI][globJ] +=
                    model_().localJacobian().mat(i,j);
            }
        }
    }

    // "assemble" a ghost element
    void assembleGhostElement_(const Element &elem)
    {
        //return; // we just ignore ghosts!
        int n = elem.template count<dim>();
        for (int i=0; i < n; ++i) {
            const VertexPointer vp = elem.template subEntity<dim>(i);
            
            // set main diagonal entries for the vertex
            int vIdx = vertexMapper_().map(*vp);
            typedef typename Matrix::block_type BlockType;
            BlockType &J = (*matrix_)[vIdx][vIdx];
            for (int j = 0; j < BlockType::rows; ++j)
                J[j][j] = 1.0;

            // set residual for the vertex
            residual_[vIdx] = 0;
        }
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
    Matrix *matrix_;
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
    std::vector<Scalar> vertexDelta_;

    int totalElems_;
    int greenElems_;
    
    Scalar nextReassembleAccuracy_;
    Scalar reassembleAccuracy_;
    
#if HAVE_DUNE_PDELAB
    // PDELab stuff
    Constraints *cn_;
    FEM *fem_;
    ScalarGridFunctionSpace *scalarGridFunctionSpace_;
    GridFunctionSpace *gridFunctionSpace_;
    ConstraintsTrafo *constraintsTrafo_;
    LocalOperator *localOperator_;
    GridOperatorSpace *gridOperatorSpace_;
#endif
};

} // namespace Dumux

#endif
