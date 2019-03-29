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
 * \ingroup Linear
 *
 * \brief Wraps the AMG backend such that it can be used for the elrichards model.
 */
#ifndef DUMUX_ELRICHARDS_AMGBACKEND_HH
#define DUMUX_ELRICHARDS_AMGBACKEND_HH

#include <dumux/linear/amgbackend.hh>

namespace Dumux {

/*!
 * \brief Base class for the ElRichards AMGBackend.
 */
template <class TypeTag, bool isParallel>
class ElRichardsAMGBackendBase : public AMGBackend<TypeTag>
{
    typedef AMGBackend<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    /*!
     * \copydoc AMGBackend::AMGBackend()
     */
    ElRichardsAMGBackendBase(const Problem& problem)
    : ParentType(problem)
    {}
};

/*!
 * \brief Specialization for the parallel setting.
 */
template <class TypeTag>
class ElRichardsAMGBackendBase<TypeTag, true> : public AMGBackend<TypeTag>
{
    typedef AMGBackend<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef typename Dune::BCRSMatrix<MatrixBlock> BlockMatrix;
    typedef typename Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef typename Dune::BlockVector<VectorBlock> BlockVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };

public:
    /*!
     * \copydoc AMGBackend::AMGBackend()
     */
    ElRichardsAMGBackendBase(const Problem& problem)
    : ParentType(problem)
    {
        createBlockMatrixAndVectors_();
    }

    /*!
     * \copydoc AMGBackend::solve()
     */
    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        flatToBlocked_(A, x, b);
        int converged = ParentType::solve(*aBlocked_, *xBlocked_, *bBlocked_);
        blockedToFlat_(x, b);
        return converged;
    }

private:
    void createBlockMatrixAndVectors_()
    {
        int numVertices = this->problem().gridView().size(dim);

        aBlocked_ = std::make_shared<BlockMatrix>(numVertices, numVertices, BlockMatrix::random);
        xBlocked_ = std::make_shared<BlockVector>(numVertices);
        bBlocked_ = std::make_shared<BlockVector>(numVertices);

        // find out the global indices of the neighboring vertices of
        // each vertex
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(numVertices);
        for (const auto& element : elements(this->problem().gridView())) {

            // loop over all element vertices
            int n = element.subEntities(dim);
            for (int i = 0; i < n - 1; ++i) {
                int globalI = this->problem().vertexMapper().subIndex(element, i, dim);
                for (int j = i + 1; j < n; ++j) {
                    int globalJ = this->problem().vertexMapper().subIndex(element, j, dim);
                    // make sure that vertex j is in the neighbor set
                    // of vertex i and vice-versa
                    neighbors[globalI].insert(globalJ);
                    neighbors[globalJ].insert(globalI);
                }
            }
        }

        // make vertices neighbors to themselfs
        for (int i = 0; i < numVertices; ++i)
            neighbors[i].insert(i);

        // allocate space for the rows of the matrix
        for (int i = 0; i < numVertices; ++i) {
            aBlocked_->setrowsize(i, neighbors[i].size());
        }
        aBlocked_->endrowsizes();

        // fill the rows with indices. each vertex talks to all of its
        // neighbors. (it also talks to itself since vertices are
        // sometimes quite egocentric.)
        for (int i = 0; i < numVertices; ++i) {
            auto nIt = neighbors[i].begin();
            const auto& nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                aBlocked_->addindex(i, *nIt);
            }
        }
        aBlocked_->endindices();
    }

    template <class FlatMatrix, class FlatVector>
    void flatToBlocked_(const FlatMatrix& aFlat,
                        const FlatVector& xFlat,
                        const FlatVector& bFlat)
    {
        unsigned numBlocks = xBlocked_->size();
        static const unsigned numMassEq = numEq - dim;
        for (unsigned rowBlockIdx = 0; rowBlockIdx < numBlocks; ++rowBlockIdx)
        {
            for (unsigned rowEqIdx = 0; rowEqIdx < numEq; ++rowEqIdx)
            {
                unsigned rowFlatIdx;
                if (rowEqIdx < numMassEq)
                    rowFlatIdx = rowBlockIdx*numMassEq + rowEqIdx;
                else
                    rowFlatIdx = numBlocks*numMassEq + rowBlockIdx*dim + rowEqIdx - numMassEq;

                (*xBlocked_)[rowBlockIdx][rowEqIdx] = xFlat[rowFlatIdx];
                (*bBlocked_)[rowBlockIdx][rowEqIdx] = bFlat[rowFlatIdx];

                for (auto colBlockIt = (*aBlocked_)[rowBlockIdx].begin();
                     colBlockIt != (*aBlocked_)[rowBlockIdx].end(); ++colBlockIt)
                {
                    unsigned colBlockIdx = colBlockIt.index();
                    auto& aBlock = (*aBlocked_)[rowBlockIdx][colBlockIdx];

                    for (unsigned colEqIdx = 0; colEqIdx < numEq; ++colEqIdx)
                    {
                        unsigned colFlatIdx;
                        if (colEqIdx < numMassEq)
                            colFlatIdx = colBlockIdx*numMassEq + colEqIdx;
                        else
                            colFlatIdx = numBlocks*numMassEq + colBlockIdx*dim + colEqIdx - numMassEq;

                        aBlock[rowEqIdx][colEqIdx] = aFlat[rowFlatIdx][colFlatIdx];
                    }
                }
            }
        }
    }

    template <class FlatVector>
    void blockedToFlat_(FlatVector& xFlat,
                        FlatVector& bFlat)
    {
        unsigned numBlocks = xBlocked_->size();
        static const unsigned numMassEq = numEq - dim;
        for (unsigned rowBlockIdx = 0; rowBlockIdx < numBlocks; ++rowBlockIdx)
        {
            for (unsigned rowEqIdx = 0; rowEqIdx < numEq; ++rowEqIdx)
            {
                unsigned rowFlatIdx;
                if (rowEqIdx < numMassEq)
                    rowFlatIdx = rowBlockIdx*numMassEq + rowEqIdx;
                else
                    rowFlatIdx = numBlocks*numMassEq + rowBlockIdx*dim + rowEqIdx - numMassEq;

                xFlat[rowFlatIdx] = (*xBlocked_)[rowBlockIdx][rowEqIdx];
                bFlat[rowFlatIdx] = (*bBlocked_)[rowBlockIdx][rowEqIdx];
            }
        }
    }

    std::shared_ptr<BlockMatrix> aBlocked_;
    std::shared_ptr<BlockVector> xBlocked_;
    std::shared_ptr<BlockVector> bBlocked_;
};

/*!
 * \brief Wraps the AMG backend such that it can be used for the elrichards model.
 */
template <class TypeTag>
class ElRichardsAMGBackend : public ElRichardsAMGBackendBase<
    TypeTag,
    Dune::Capabilities::canCommunicate<typename GET_PROP_TYPE(TypeTag, Grid),
                                       GET_PROP_TYPE(TypeTag, Grid)::dimension>::v>
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    enum { dofCodim = Grid::dimension };
    enum { isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v };
    typedef ElRichardsAMGBackendBase<TypeTag, isParallel> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    /*!
     * \copydoc AMGBackend::AMGBackend()
     */
    ElRichardsAMGBackend(const Problem& problem)
    : ParentType(problem)
    {}
};

} // namespace Dumux

#endif // DUMUX_ELRichards_AMGBACKEND_HH
