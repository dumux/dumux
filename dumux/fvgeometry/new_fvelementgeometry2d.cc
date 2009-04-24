// $Id$

#include "config.h"
#include <dune/grid/common/capabilities.hh>

/*!
 * \brief Specialization of the FVElementGeometry for the 2D-case.
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
                           2>
    : public NewFVElementGeometryBase<GridT,
                                      ReferenceElementContainerT,
                                      ShapeFunctionSetContainerT,
                                      NewFVElementGeometry<GridT,
                                                           ReferenceElementContainerT,
                                                           ShapeFunctionSetContainerT,
                                                           2> >
{

    enum {
        GridDim = Grid::dimension,
        WorldDim = Grid::dimensionworld,

        // TODO/low priority: This should be specific to the actual
        // reference elements used.
        MaxNodes = 4,
        MaxEdges = 4,
        MaxFaces = 0,
        MaxBoundaryFaces = 8,
        MaxNodesInFace = 2
    };

    typedef NewFVElementGeometry<Grid,
                                 ReferenceElementContainer,
                                 ShapeFunctionSetContainer>  ThisType;
    typedef NewFVElementGeometryBase<Grid,
                                     ReferenceElementContainer,
                                     ShapeFunctionSetContainer,
                                     ThisType>               ParentType;

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

    typedef Dune::FieldVector<CoordScalar, GridDim>   LocalCoord;
    typedef Dune::FieldVector<CoordScalar, WorldDim>  WorldCoord;
    typedef FieldMatrix<CoordScalar,GridDim,WorldDim> InverseTransposedJacobian;

    // allow the parent class to access our protected member
    // functions

    /** \todo Please doc me! */

    friend class NewFVElementGeometryBase< ThisType >;

    Scalar quadrilateralArea_(const WorldCoord& p0,
                              const WorldCoord& p1,
                              const WorldCoord& p2,
                              const WorldCoord& p3)
    {
        if (!Capabilities::IsUnstructured<Grid>::v) {
            return this->cellVolume_ / 4.0;
        }
        else
            return 0.5*fabs((p3[0] - p1[0])*(p2[1] - p0[1]) -
                            (p3[1] - p1[1])*(p2[0] - p0[0]));
    }

    void updateVolumes_(const Cell &cell)
    {
        switch (cell.template count<GridDim>()) {
        case 3: // 2D, triangle
            for (int k = 0; k < 3; k++)
                this->subContVol_[k].volume = this->cellVolume_ / 3.0;
            break;

        case 4: // 2D, quadrilinear
            this->subContVol_[0].volume =
                quadrilateralArea_(this->subContVol_[0].global,
                                   this->edgeCoord_[2],
                                   this->elementGlobal_,
                                   this->edgeCoord_[0]);
            this->subContVol_[1].volume =
                quadrilateralArea_(this->subContVol_[1].global,
                                   this->edgeCoord_[1],
                                   this->elementGlobal_,
                                   this->edgeCoord_[2]);
            this->subContVol_[2].volume =
                quadrilateralArea_(this->subContVol_[2].global,
                                   this->edgeCoord_[0],
                                   this->elementGlobal_,
                                   this->edgeCoord_[3]);
            this->subContVol_[3].volume =
                quadrilateralArea_(this->subContVol_[3].global,
                                   this->edgeCoord_[3],
                                   this->elementGlobal_,
                                   this->edgeCoord_[1]);
            break;

        default:
            DUNE_THROW(NotImplemented,
                       "updateVolumes_ dim = "
                       << GridDim
                       << ", numVertices = "
                       << cell.template count<GridDim>());
        }
    };

    void updateInteriorFaces_(const ShapeFunctionSet &sfs,
                              const ReferenceElement &refElem,
                              const CellGeometry &geometry)
    {
        LocalCoord ipLocal;
        WorldCoord diffVec;

        // fill the faces between the sub-control-volumes
        for (int scvFaceIdx = 0; scvFaceIdx < this->numEdges_; scvFaceIdx++)
        {
            int insideIdx = refElem.subEntity(scvFaceIdx, GridDim-1, 0, GridDim);
            int outsideIdx = refElem.subEntity(scvFaceIdx, GridDim-1, 1, GridDim);

            this->subContVolFace_[scvFaceIdx].insideIdx_ = insideIdx;
            this->subContVolFace_[scvFaceIdx].outsideIdx_ = outsideIdx;

            // calculate the local integration point and
            // the face normal.
            //
            // TODO: use a reference element instead of
            // calculating this every time
            ipLocal  = refElem.position(scvFaceIdx, GridDim-1);
            ipLocal += this->elementLocal_;
            ipLocal *= 0.5;

            this->subContVolFace_[scvFaceIdx].ipLocal = ipLocal;

            diffVec = this->elementGlobal_ - this->edgeCoord[scvFaceIdx];
            this->subContVolFace_[scvFaceIdx].normal[0] =  diffVec[1];
            this->subContVolFace_[scvFaceIdx].normal[1] = -diffVec[0];
        }
    }

    void updateBoundaryFace_(IntersectionIterator &it,
                             const ReferenceElement &refElem,
                             const CellGeometry &geometry)
    {
        // fill boundary face data:
        int face = it->indexInInside();
        int numVerticesOfFace = refElem.size(face, 1, GridDim);
        for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
        {
            int nodeInElement = refElem.subEntity(face, 1, nodeInFace, GridDim);
            int bfIndex = this->_boundaryFaceIndex(face, nodeInFace);

            // TODO: use reference element to find the local
            // integration point
            this->boundaryFace_[bfIndex].ipLocal_ = refElem.position(nodeInElement, GridDim);
            this->boundaryFace_[bfIndex].ipLocal_ += refElem.position(face, 1);
            this->boundaryFace_[bfIndex].ipLocal_ *= 0.5;

            this->boundaryFace_[bfIndex].area_ = 0.5*it->geometry().volume();

            this->boundaryFace_[bfIndex].ipGlobal_ =
                geometry.global(this->boundaryFace_[bfIndex].ipLocal_);
        }
    }

};
