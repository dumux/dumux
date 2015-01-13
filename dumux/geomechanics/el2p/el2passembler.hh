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
 *
 * \brief This file contains an assembler for the Jacobian matrix
 * of the two-phase linear-elastic model based on PDELab.
 */
#ifndef DUMUX_EL2P_ASSEMBLER_HH
#define DUMUX_EL2P_ASSEMBLER_HH

#include <dune/common/version.hh>
#include "el2pproperties.hh"
#include "el2plocaloperator.hh"

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(PressureFEM); //!< Finite element space used for pressure, saturation, ...
NEW_PROP_TAG(DisplacementFEM); //!< Finite element space used for displacement
NEW_PROP_TAG(PressureGridFunctionSpace); //!< Grid function space used for pressure, saturation, ...
NEW_PROP_TAG(DisplacementGridFunctionSpace); //!< Grid function space used for displacement
}

namespace PDELab {

/*!
 * \brief An assembler for the Jacobian matrix
 * of the two-phase linear-elastic model based on PDELab.
 */
template<class TypeTag>
class El2PAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;

#if HAVE_DUNE_PDELAB
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureFEM)) PressureFEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureGridFunctionSpace)) PressureGFS;
    typedef typename PressureGFS::template Child<0>::Type PressureScalarGFS;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DisplacementFEM)) DisplacementFEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DisplacementGridFunctionSpace)) DisplacementGFS;
    typedef typename DisplacementGFS::template Child<0>::Type DisplacementScalarGFS;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Constraints)) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalOperator)) LocalOperator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridOperator)) GridOperator;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    enum{dim = GridView::dimension};
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GridView::template Codim<dim>::Entity Vertex;

    enum {
        enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(ImplicitEnablePartialReassemble)),
        enableJacobianRecycling = GET_PROP_VALUE(TypeTag, PTAG(ImplicitEnableJacobianRecycling)),
    };

    // copying the jacobian assembler is not a good idea
    El2PAssembler(const El2PAssembler &);

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
        Red,

        /*!
         * Vertex/element that needs to be reassembled because a
         * neighboring element/vertex is red
         */
        Yellow,

        /*!
         * Yellow vertex has only non-green neighbor elements.
         *
         * This means that its relative error is below the tolerance,
         * but its defect can be linearized without any additional
         * cost. This is just an "internal" color which is not used
         * ouside of the jacobian assembler.
         */
        Orange,

        /*!
         * Vertex/element that does not need to be reassembled
         */
        Green
    };

    El2PAssembler()
    {
        // set reassemble tolerance to 0, so that if partial
        // reassembly of the jacobian matrix is disabled, the
        // reassemble tolerance is always smaller than the current
        // relative tolerance
        reassembleTolerance_ = 0.0;
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

        constraints_ = Dune::make_shared<Constraints>();

        pressureFEM_ = Dune::make_shared<PressureFEM>(problemPtr_->gridView());
        pressureScalarGFS_ = Dune::make_shared<PressureScalarGFS>(problemPtr_->gridView(), *pressureFEM_, *constraints_);
        pressureGFS_ = Dune::make_shared<PressureGFS>(*pressureScalarGFS_);

        displacementFEM_ = Dune::make_shared<DisplacementFEM>(problemPtr_->gridView());
        displacementScalarGFS_ = Dune::make_shared<DisplacementScalarGFS>(problemPtr_->gridView(), *displacementFEM_, *constraints_);
        displacementGFS_ = Dune::make_shared<DisplacementGFS>(*displacementScalarGFS_);

        gridFunctionSpace_ = Dune::make_shared<GridFunctionSpace>(*pressureGFS_, *displacementGFS_);

        constraintsTrafo_ = Dune::make_shared<ConstraintsTrafo>();

        // initialize the grid operator spaces
        localOperator_ = Dune::make_shared<LocalOperator>(problemPtr_->model());
        gridOperator_ =
            Dune::make_shared<GridOperator>(*gridFunctionSpace_, *constraintsTrafo_,
                                  *gridFunctionSpace_, *constraintsTrafo_, *localOperator_);

        // allocate raw matrix
        matrix_ = Dune::make_shared<JacobianMatrix>(*gridOperator_);

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseMatrix_ = false;

        residual_ = Dune::make_shared<SolutionVector>(*gridFunctionSpace_);

        int numVertices = gridView_().size(dim);
        int numElements = gridView_().size(0);

        totalElems_ = gridView_().comm().sum(numElements);

        // initialize data needed for partial reassembly
        if (enablePartialReassemble) {
            vertexColor_.resize(numVertices);
            vertexDelta_.resize(numVertices);
            elementColor_.resize(numElements);
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
        *matrix_ = 0;
        gridOperator_->jacobian(problemPtr_->model().curSol(), *matrix_);

        *residual_ = 0;
        gridOperator_->residual(problemPtr_->model().curSol(), *residual_);

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
        nextReassembleTolerance_ = 0.0;

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
     * \brief Returns the relative error below which a vertex is
     *        considered to be "green" if partial Jacobian reassembly
     *        is enabled.
     *
     * This returns the _actual_ relative computed seen by
     * computeColors(), not the tolerance which it was given.
     */
    Scalar reassembleTolerance() const
    { return reassembleTolerance_; }

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

        ElementIterator eIt = gridView_().template begin<0>();
        ElementIterator eEndIt = gridView_().template end<0>();

        // mark the red vertices and update the tolerance of the
        // linearization which actually will get achieved
        nextReassembleTolerance_ = 0;
        for (int i = 0; i < vertexColor_.size(); ++i) {
            vertexColor_[i] = Green;
            if (vertexDelta_[i] > relTol) {
                // mark vertex as red if discrepancy is larger than
                // the relative tolerance
                vertexColor_[i] = Red;
            }
            nextReassembleTolerance_ =
                std::max(nextReassembleTolerance_, vertexDelta_[i]);
        };

        // Mark all red elements
        for (; eIt != eEndIt; ++eIt) {
            // find out whether the current element features a red
            // vertex
            bool isRed = false;
            int numVertices = eIt->template count<dim>();
            for (int i=0; i < numVertices; ++i) {
                int globalI = vertexMapper_().map(*eIt, i, dim);
                if (vertexColor_[globalI] == Red) {
                    isRed = true;
                    break;
                }
            };

            // if yes, the element color is also red, else it is not
            // red, i.e. green for the mean time
            int eIdxGlobal = elementMapper_().map(*eIt);
            if (isRed)
                elementColor_[eIdxGlobal] = Red;
            else
                elementColor_[eIdxGlobal] = Green;
        }

        // Mark yellow vertices (as orange for the mean time)
        eIt = gridView_().template begin<0>();
        for (; eIt != eEndIt; ++eIt) {
            int eIdx = this->elementMapper_().map(*eIt);
            if (elementColor_[eIdx] != Red)
                continue; // non-red elements do not tint vertices
                          // yellow!

            int numVertices = eIt->template count<dim>();
            for (int i=0; i < numVertices; ++i) {
                int globalI = vertexMapper_().map(*eIt, i, dim);
                // if a vertex is already red, don't recolor it to
                // yellow!
                if (vertexColor_[globalI] != Red)
                    vertexColor_[globalI] = Orange;
            };
        }

        // Mark yellow elements
        eIt = gridView_().template begin<0>();
        for (; eIt != eEndIt; ++eIt) {
            int eIdx = this->elementMapper_().map(*eIt);
            if (elementColor_[eIdx] == Red) {
                continue; // element is red already!
            }

            // check whether the element features a yellow
            // (resp. orange at this point) vertex
            bool isYellow = false;
            int numVertices = eIt->template count<dim>();
            for (int i=0; i < numVertices; ++i) {
                int globalI = vertexMapper_().map(*eIt, i, dim);
                if (vertexColor_[globalI] == Orange) {
                    isYellow = true;
                    break;
                }
            };

            if (isYellow)
                elementColor_[eIdx] = Yellow;
        }

        // Demote orange vertices to yellow ones if it has at least
        // one green element as a neighbor.
        eIt = gridView_().template begin<0>();
        for (; eIt != eEndIt; ++eIt) {
            int eIdx = this->elementMapper_().map(*eIt);
            if (elementColor_[eIdx] != Green)
                continue; // yellow and red elements do not make
                          // orange vertices yellow!

            int numVertices = eIt->template count<dim>();
            for (int i=0; i < numVertices; ++i) {
                int globalI = vertexMapper_().map(*eIt, i, dim);
                // if a vertex is orange, recolor it to yellow!
                if (vertexColor_[globalI] == Orange)
                    vertexColor_[globalI] = Yellow;
            };
        }

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
     * \param vIdx The local index of the vertex in the element.
     */
    int vertexColor(const Element &element, int vIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int vIdxGlobal = vertexMapper_().subIndex(element, vIdx, dim);
#else 
        int vIdxGlobal = vertexMapper_().map(element, vIdx, dim);
#endif
        return vertexColor_[vIdxGlobal];
    }

    /*!
     * \brief Returns the reassemble color of a vertex
     *
     * \param vIdxGlobal The global index of the vertex.
     */
    int vertexColor(int vIdxGlobal) const
    {
        if (!enablePartialReassemble)
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
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!

        int eIdxGlobal = elementMapper_().map(element);
        return elementColor_[eIdxGlobal];
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
    const JacobianMatrix& matrix() const
    { return *matrix_; }

    /*!
     * \brief Return reference to global Jacobian matrix.
     */
    JacobianMatrix& matrix()
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return *residual_; }


    /*!
     * \brief Return reference to global residual vector.
     */
    SolutionVector& residual()
    { return *residual_; }


private:
#if !HAVE_DUNE_PDELAB
    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int numVertices = gridView_().size(dim);

        // allocate raw matrix
        matrix_ = Dune::make_shared<JacobianMatrix>(numVertices, numVertices, JacobianMatrix::random);

        // find out the global indices of the neighboring vertices of
        // each vertex
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(numVertices);
        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            const Element &element = *eIt;

            // loop over all element vertices
            int n = element.template count<dim>();
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
        for (int i = 0; i < numVertices; ++i)
            neighbors[i].insert(i);

        // allocate space for the rows of the matrix
        for (int i = 0; i < numVertices; ++i) {
            matrix_->setrowsize(i, neighbors[i].size());
        }
        matrix_->endrowsizes();

        // fill the rows with indices. each vertex talks to all of its
        // neighbors. (it also talks to itself since vertices are
        // sometimes quite egocentric.)
        for (int i = 0; i < numVertices; ++i) {
            typename NeighborSet::iterator nIt = neighbors[i].begin();
            typename NeighborSet::iterator nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                matrix_->addindex(i, *nIt);
            }
        }
        matrix_->endindices();
    };
#endif

    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // always reset the right hand side.
        *residual_ = 0.0;

        if (!enablePartialReassemble) {
            // If partial reassembly of the jacobian is not enabled,
            // we can just reset everything!
            (*matrix_) = 0;
            return;
        }

        // reset all entries corrosponding to a red vertex
        for (int rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if (vertexColor_[rowIdx] == Green)
                continue; // the equations for this control volume are
                          // already below the treshold

            // set all entries in the row to 0
            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[rowIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                (*colIt) = 0.0;
            }
        };
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
    Dune::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    Dune::shared_ptr<SolutionVector> residual_;

    // attributes required for jacobian matrix recycling
    bool reuseMatrix_;

    // attributes required for partial jacobian reassembly
    std::vector<EntityColor> vertexColor_;
    std::vector<EntityColor> elementColor_;
    std::vector<Scalar> vertexDelta_;

    int totalElems_;
    int greenElems_;

    Scalar nextReassembleTolerance_;
    Scalar reassembleTolerance_;


    Dune::shared_ptr<Constraints> constraints_;
    Dune::shared_ptr<PressureFEM> pressureFEM_;
    Dune::shared_ptr<DisplacementFEM> displacementFEM_;
    Dune::shared_ptr<PressureScalarGFS> pressureScalarGFS_;
    Dune::shared_ptr<DisplacementScalarGFS> displacementScalarGFS_;
    Dune::shared_ptr<PressureGFS> pressureGFS_;
    Dune::shared_ptr<DisplacementGFS> displacementGFS_;
    Dune::shared_ptr<GridFunctionSpace> gridFunctionSpace_;
    Dune::shared_ptr<ConstraintsTrafo> constraintsTrafo_;
    Dune::shared_ptr<LocalOperator> localOperator_;
    Dune::shared_ptr<GridOperator> gridOperator_;
};

} // namespace PDELab

} // namespace Dumux

#endif
