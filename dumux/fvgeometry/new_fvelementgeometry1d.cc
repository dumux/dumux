// $Id$

/*!
 * \brief Specialization of the FVElementGeometry for the 1D-case.
 *
 * The complete public API of this class is inherited from the
 * base class, so this class only fills some gaps.
 */
template <class GridT,
          class ReferenceElementContainerT,
          class ShapeFunctionSetContainerT>
class NewFVElementGeometry<GridT,
                           ReferenceElementContainerT,
                           ShapeFunctionSetContainerT,
                           1>
    : public NewFVElementGeometryBase<NewFVElementGeometry<GridT,
                                                           ReferenceElementContainerT,
                                                           ShapeFunctionSetContainerT,
                                                           1> >
{
    typedef GridT                                           Grid;
    typedef typename Grid::template Codim<0>::Entity        Cell;
    typedef ReferenceElementContainerT                      ReferenceElementContainer;
    typedef typename ReferenceElementContainer::value_type  ReferenceElement;
    typedef ShapeFunctionSetContainerT                      ShapeFunctionSetContainer;
    typedef typename ShapeFunctionSetContainer::value_type  ShapeFunctionSet;
    typedef typename Cell::Geometry                         CellGeometry;

    typedef Dune::IntersectionIteratorGetter<Grid,Dune::LeafTag>      IntersectionIteratorGetter;
    typedef typename IntersectionIteratorGetter::IntersectionIterator IntersectionIterator;

    typedef typename ShapeFunctionSet::ResultType           Scalar;
    typedef typename Grid::ctype                            CoordScalar;

    enum {
        GridDim = Grid::dimension,
        WorldDim = Grid::dimensionworld,

        MaxNodes = 2,
        MaxEdges = 1,
        MaxFaces = 0,
        MaxBoundaryFaces = 2,
        MaxNodesInFace = 1
    };

    typedef NewFVElementGeometry<Grid,
                                 ReferenceElementContainer,
                                 ShapeFunctionSetContainer>  ThisType;
    typedef NewFVElementGeometryBase<ThisType>               ParentType;

    typedef Dune::FieldVector<CoordScalar, GridDim>   LocalCoord;
    typedef Dune::FieldVector<CoordScalar, WorldDim>  WorldCoord;
    typedef FieldMatrix<CoordScalar,GridDim,WorldDim> InverseTransposedJacobian;

    // allow the parent class to access our protected member
    // functions

    /** \todo Please doc me! */

    friend class NewFVElementGeometryBase< ThisType >;

    void updateVolumes_(const Cell &cell)
        {
            this->subContVol_[0].volume = 0.5*this->cellVolume_;
            this->subContVol_[1].volume = 0.5*this->cellVolume_;
        };

    void updateInteriorFaces_(const ShapeFunctionSet &sfs,
                              const ReferenceElement &refElem,
                              const CellGeometry &geometry)
        {
            WorldCoord diffVec;
            LocalCoord temp;

            // fill sub control volume face data
            for (int scvFaceIdx = 0; scvFaceIdx < this->numEdges_; scvFaceIdx++)
            {
                int insideIdx  = refElem.subEntity(scvFaceIdx, GridDim-1, 0, GridDim);
                int outsideIdx = refElem.subEntity(scvFaceIdx, GridDim-1, 1, GridDim);

                this->subContVolFace_[scvFaceIdx].insideIdx_ = insideIdx;
                this->subContVolFace_[scvFaceIdx].outsideIdx_ = outsideIdx;

                // calculate the local integration point and
                // the face normal.
                this->subContVolFace_[scvFaceIdx].ipLocal = 0.5;
                this->subContVolFace_[scvFaceIdx].normal = 1.0;

                // get the global integration point and the Jacobian inverse
                const LocalCoord &ipLocal = this->subContVolFace_[scvFaceIdx].ipLocal;
                this->subContVolFace_[scvFaceIdx].ipGlobal = geometry.global(ipLocal);
                const InverseTransposedJacobian &jacInvT = geometry.jacobianInverseTransposed(ipLocal);

                // calculate the gradients at the centers of the sub control volumes
                for (int node = 0; node < this->numNodes_; node++) {
                    for (int i = 0; i < GridDim; i++)
                        temp[i] = sfs[node].evaluateDerivative(0,
                                                               i,
                                                               ipLocal);
                    // grad[node] = J^-1 * temp
                    jacInvT.mv(temp,
                               this->subContVolFace_[scvFaceIdx].grad[node]);
                }
            }
        }

    void updateBoundaryFace_(IntersectionIterator &it,
                             const ReferenceElement &refElem,
                             const CellGeometry &geometry)
        {
            // fill boundary face data:
            int face = it->numberInSelf();
            int numVerticesOfFace = refElem.size(face, 1, GridDim);
            for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
            {
                int nodeInElement = refElem.subEntity(face, 1, nodeInFace, GridDim);
                int bfIndex = this->boundaryFaceIndex(face, nodeInFace);

                this->boundaryFace_[bfIndex].ipLocal = refElem.position(nodeInElement, GridDim);
                this->boundaryFace_[bfIndex].area = 1.0;

                this->boundaryFace_[bfIndex].ipGlobal =
                    geometry.global(this->boundaryFace_[bfIndex].ipLocal);
            }

        }
};

