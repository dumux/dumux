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
 * \brief An assembler for the global Jacobian matrix for models using the box discretization.
 */
#ifndef DUMUX_BOX_ASSEMBLER_HH
#define DUMUX_BOX_ASSEMBLER_HH

#include <dumux/common/exceptions.hh>
#include <dumux/implicit/assembler.hh>
#include <dumux/parallel/vertexhandles.hh>

namespace Dumux {

/*!
 * \ingroup BoxModel
 * \brief An assembler for the global Jacobian matrix for models using the box discretization.
 */
template<class TypeTag>
class BoxAssembler : public ImplicitAssembler<TypeTag>
{
    typedef ImplicitAssembler<TypeTag> ParentType;
    friend class ImplicitAssembler<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum{ dim = GridView::dimension };

public:
    BoxAssembler(): ParentType() {}

private:
    // copying the jacobian assembler is not a good idea
    BoxAssembler(const BoxAssembler &);

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
    void computeColors_(const Scalar relTol)
    {
        if (!this->enablePartialReassemble_())
            return;

        // mark the red vertices and update the tolerance of the
        // linearization which actually will get achieved
        this->nextReassembleAccuracy_ = 0;
        for (unsigned int i = 0; i < this->vertexColor_.size(); ++i) {
            using std::max;
            if (this->delta_[i] > relTol)
                // mark vertex as red if discrepancy is larger than
                // the relative tolerance
                this->vertexColor_[i] = ParentType::Red;
            else
                this->nextReassembleAccuracy_ =
                    max(this->nextReassembleAccuracy_, this->delta_[i]);
        }

        // Mark all red elements
        for (const auto& element : elements(this->gridView_())) {
            // find out whether the current element features a red
            // vertex
            bool isRed = false;

            int numVertices = element.subEntities(dim);

            for (int i=0; i < numVertices; ++i) {
                int globalI = this->vertexMapper_().subIndex(element, i, dim);

                if (this->vertexColor_[globalI] == ParentType::Red) {
                    isRed = true;
                    break;
                }
            }

            // if yes, the element color is also red, else it is not
            // red, i.e. green for the mean time
            int eIdxGlobal = this->elementMapper_().index(element);

            if (isRed)
                this->elementColor_[eIdxGlobal] = ParentType::Red;
            else
                this->elementColor_[eIdxGlobal] = ParentType::Green;
        }

        // Mark yellow vertices (as orange for the mean time)
        for (const auto& element : elements(this->gridView_())) {
            int eIdx = this->elementMapper_().index(element);

            if (this->elementColor_[eIdx] != ParentType::Red)
                continue; // non-red elements do not tint vertices
                          // yellow!

            int numVertices = element.subEntities(dim);

            for (int i = 0; i < numVertices; ++i) {
                int globalI = this->vertexMapper_().subIndex(element, i, dim);

                // if a vertex is already red, don't recolor it to
                // yellow!
                if (this->vertexColor_[globalI] != ParentType::Red) {
                    this->vertexColor_[globalI] = ParentType::Orange;
                }
            }
        }

        // at this point we communicate the yellow vertices to the
        // neighboring processes because a neigbor process may not see
        // the red vertex for yellow border vertices
        typedef typename ParentType::EntityColor EntityColor;
        VertexHandleMin<EntityColor, std::vector<EntityColor>,  VertexMapper>
            minHandle(this->vertexColor_, this->vertexMapper_());
        this->gridView_().communicate(minHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // Mark yellow elements
        for (const auto& element : elements(this->gridView_())) {
            int eIdx = this->elementMapper_().index(element);

            if (this->elementColor_[eIdx] == ParentType::Red) {
                continue; // element is red already!
            }

            // check whether the element features a yellow
            // (resp. orange at this point) vertex
            bool isYellow = false;
            int numVertices = element.subEntities(dim);

            for (int i = 0; i < numVertices; ++i) {
                int globalI = this->vertexMapper_().subIndex(element, i, dim);

                if (this->vertexColor_[globalI] == ParentType::Orange) {
                    isYellow = true;
                    break;
                }
            }

            if (isYellow)
                this->elementColor_[eIdx] = ParentType::Yellow;
        }

        // Demote orange vertices to yellow ones if it has at least
        // one green element as a neighbor.
        for (const auto& element : elements(this->gridView_())) {
            int eIdx = this->elementMapper_().index(element);

            if (this->elementColor_[eIdx] != ParentType::Green)
                continue; // yellow and red elements do not make
                          // orange vertices yellow!

            int numVertices = element.subEntities(dim);

            for (int i = 0; i < numVertices; ++i) {
                int globalI = this->vertexMapper_().subIndex(element, i, dim);

                // if a vertex is orange, recolor it to yellow!
                if (this->vertexColor_[globalI] == ParentType::Orange)
                    this->vertexColor_[globalI] = ParentType::Yellow;
            }
        }

        // demote the border orange vertices
        VertexHandleMax<EntityColor, std::vector<EntityColor>,  VertexMapper>
            maxHandle(this->vertexColor_,
                      this->vertexMapper_());
        this->gridView_().communicate(maxHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // promote the remaining orange vertices to red
        for (unsigned int i=0; i < this->vertexColor_.size(); ++i) {
            // if a vertex is green or yellow don't do anything!
            if (this->vertexColor_[i] == ParentType::Green || this->vertexColor_[i] == ParentType::Yellow)
                continue;

            // make sure the vertex is red (this is a no-op vertices
            // which are already red!)
            this->vertexColor_[i] = ParentType::Red;

            // set the error of this vertex to 0 because the system
            // will be consistently linearized at this vertex
            this->delta_[i] = 0.0;
        }
    }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int numVerticesGlobal = this->gridView_().size(dim);

        // allocate raw matrix
        this->matrix_ = std::make_shared<JacobianMatrix>(numVerticesGlobal, numVerticesGlobal, JacobianMatrix::random);

        // find out the global indices of the neighboring vertices of
        // each vertex
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(numVerticesGlobal);
        for (const auto& element : elements(this->gridView_())) {
            int numVerticesLocal = element.subEntities(dim);

            // if the element is not in the interior or the process
            // border, all dofs just contain main-diagonal entries
            if (element.partitionType() != Dune::InteriorEntity &&
                element.partitionType() != Dune::BorderEntity)
            {
                for (int i = 0; i < numVerticesLocal; ++i) {
                    int globalI = this->vertexMapper_().subIndex(element, i, dim);

                    neighbors[globalI].insert(globalI);
                }
            }
            else
            {
                // loop over all element vertices
                for (int i = 0; i < numVerticesLocal - 1; ++i) {
                    int globalI = this->vertexMapper_().subIndex(element, i, dim);

                    for (int j = i + 1; j < numVerticesLocal; ++j) {
                        int globalJ = this->vertexMapper_().subIndex(element, j, dim);

                        // make sure that vertex j is in the neighbor set
                        // of vertex i and vice-versa
                        neighbors[globalI].insert(globalJ);
                        neighbors[globalJ].insert(globalI);
                    }
                }
            }
        }

        // make vertices neighbors to themselfs
        for (int i = 0; i < numVerticesGlobal; ++i) {
            neighbors[i].insert(i);
        }

        // allocate space for the rows of the matrix
        for (int i = 0; i < numVerticesGlobal; ++i) {
            this->matrix_->setrowsize(i, neighbors[i].size());
        }
        this->matrix_->endrowsizes();

        // fill the rows with indices. each vertex talks to all of its
        // neighbors. (it also talks to itself since vertices are
        // sometimes quite egocentric.)
        for (int i = 0; i < numVerticesGlobal; ++i) {
            typename NeighborSet::iterator nIt = neighbors[i].begin();
            typename NeighborSet::iterator nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                this->matrix_->addindex(i, *nIt);
            }
        }
        this->matrix_->endindices();
    }

    // assemble a non-ghost element
    void assembleElement_(const Element &element)
    {
        if (this->enablePartialReassemble_()) {
            int eIdxGlobal = this->model_().elementMapper().index(element);

            if (this->elementColor_[eIdxGlobal] == ParentType::Green) {
                ++this->greenElems_;

                assembleGreenElement_(element);
                return;
            }
        }

        this->model_().localJacobian().assemble(element);

        int numVerticesLocal = element.subEntities(dim);

        for (int i=0; i < numVerticesLocal; ++ i) {
            int globI = this->vertexMapper_().subIndex(element, i, dim);

            // update the right hand side
            this->residual_[globI] += this->model_().localJacobian().residual(i);
            for (int j = 0; j < this->residual_[globI].dimension; ++j) {
                using std::isfinite;
                if (!isfinite(this->residual_[globI][j])) {
                    DUNE_THROW(NumericalProblem,
                               "residual_[" << globI << "][" << j << "] is not finite");
                }
            }
            if (this->enableJacobianRecycling_()) {
                this->storageTerm_[globI] +=
                    this->model_().localJacobian().storageTerm(i);
            }

            // only update the jacobian matrix for non-green vertices
            if (this->vertexColor(globI) != ParentType::Green) {
                if (this->enableJacobianRecycling_())
                    this->storageJacobian_[globI] +=
                        this->model_().localJacobian().storageJacobian(i);

                // update the jacobian matrix
                for (int j = 0; j < numVerticesLocal; ++ j) {
                    int globJ = this->vertexMapper_().subIndex(element, j, dim);

                    (*this->matrix_)[globI][globJ] +=
                        this->model_().localJacobian().mat(i,j);
                }
            }
        }
    }

    // "assemble" a green element. green elements only get the
    // residual updated, but the jacobian is left alone...
    void assembleGreenElement_(const Element &element)
    {
        this->model_().localResidual().eval(element);

        int numVerticesLocal = element.subEntities(dim);

        for (int i = 0; i < numVerticesLocal; ++ i) {
            int globI = this->vertexMapper_().subIndex(element, i, dim);

            // update the right hand side
            this->residual_[globI] += this->model_().localResidual().residual(i);
            if (this->enableJacobianRecycling_())
                this->storageTerm_[globI] += this->model_().localResidual().storageTerm(i);
        }
    }

    // "assemble" a ghost element
    void assembleGhostElement_(const Element &element)
    {
        int numVerticesLocal = element.subEntities(dim);

        for (int i=0; i < numVerticesLocal; ++i) {
            auto vertex = element.template subEntity<dim>(i);

            if (vertex.partitionType() == Dune::InteriorEntity ||
                vertex.partitionType() == Dune::BorderEntity)
            {
                // do not change the non-ghost vertices
                continue;
            }

            // set main diagonal entries for the vertex
            int vIdx = this->vertexMapper_().index(vertex);

            typedef typename JacobianMatrix::block_type BlockType;
            BlockType &J = (*this->matrix_)[vIdx][vIdx];
            for (int j = 0; j < BlockType::rows; ++j)
                J[j][j] = 1.0;

            // set residual for the vertex
            this->residual_[vIdx] = 0;
        }
    }
};

} // namespace Dumux

#endif
