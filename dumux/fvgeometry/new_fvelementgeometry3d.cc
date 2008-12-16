// $Id$ 

/*!
 * \brief Specialization of the FVElementGeometry for the 3D-case.
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
                           3>
    : public NewFVElementGeometryBase<NewFVElementGeometry<GridT,
                                                           ReferenceElementContainerT, 
                                                           ShapeFunctionSetContainerT,
                                                           3> >
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
            
        // TODO/low priority: This should be specific to the actual
        // reference elements used.
        MaxNodes = 8,
        MaxEdges = 12,
        MaxFaces = 6,
        MaxBoundaryFaces = 24,
        MaxNodesInFace = 4
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
    friend class NewFVElementGeometryBase< ThisType >;

    void crossProduct_(WorldCoord& c, 
                       const WorldCoord& a,
                       const WorldCoord& b)
        {
            c[0] = a[1]*b[2] - a[2]*b[1];
            c[1] = a[2]*b[0] - a[0]*b[2];
            c[2] = a[0]*b[1] - a[1]*b[0];
        }

    Scalar pyramidVolume_(const WorldCoord& p0,
                          const WorldCoord& p1,
                          const WorldCoord& p2,
                          const WorldCoord& p3,
                          const WorldCoord& p4)
        {
/*
  WorldCoord a = p2 - p0;
  WorldCoord b = p3 - p1;
  WorldCoord h = p4 - p0;
  WorldCoord n = crossProduct(a, b);
  return 1.0/6.0*(n*h);
*/

            WorldCoord a(p2); a -= p0;
            WorldCoord b(p3); b -= p1;
                
            WorldCoord n;
            crossProduct_(n, a, b);
                
            a = p4; a -= p0;

            return 1.0/6.0*(n*a);
        }

    Scalar prismVolume_(const WorldCoord& p0,
                        const WorldCoord& p1,
                        const WorldCoord& p2,
                        const WorldCoord& p3,
                        const WorldCoord& p4,
                        const WorldCoord& p5)
        {
/*
  WorldCoord a = p4 - p0;
  WorldCoord b = p1 - p3;
  WorldCoord c = p1 - p0;
  WorldCoord d = p2 - p0;
  WorldCoord e = p5 - p0;
  WorldCoord m = crossProduct_(a, b);
  WorldCoord n = m + crossProduct_(c, d);

  return fabs(1.0/6.0*(n*e));
*/

            WorldCoord a(p4); a -= p0;
            WorldCoord b(p1); b -= p3;
            WorldCoord m;
            crossProduct_(m, a, b);

            a = p1; a -= p0;
            b = p2; b -= p0;
            WorldCoord n;
            crossProduct_(n, a, b);
            n += m;

            a = p5; a -= p0;

            return fabs(1.0/6.0*(n*a));
        }

    void normalOfQuadrilateral_(WorldCoord &normal, 
                                const WorldCoord& p0,
                                const WorldCoord& p1,
                                const WorldCoord& p2,
                                const WorldCoord& p3)
        {
            WorldCoord a(p2); a -= p0;
            WorldCoord b(p3); b -= p1;
            crossProduct_(normal, a, b);
            normal *= 0.5;
        }

    Scalar quadrilateralArea_(const WorldCoord& p0,
                              const WorldCoord& p1,
                              const WorldCoord& p2,
                              const WorldCoord& p3)
        {
            WorldCoord normal;
            normalOfQuadrilateral3D_(normal, p0, p1, p2, p3);
            return normal.two_norm();
        }

    void getFaceIndices_(int numVertices,
                         int k,
                         int& leftFace,
                         int& rightFace)
        {
            static const int edgeToFaceTet[2][6] = {
                {2, 0, 3, 1, 2, 0},
                {3, 3, 1, 2, 0, 1}
            };
            static const int edgeToFacePyramid[2][8] = {
                {1, 2, 3, 4, 4, 1, 2, 3},
                {0, 0, 0, 0, 1, 2, 3, 4}
            };
            static const int edgeToFacePrism[2][9] = {
                {1, 2, 3, 3, 1, 2, 4, 4, 4},
                {0, 0, 0, 1, 2, 3, 1, 2, 3}
            };
            static const int edgeToFaceHex[2][12] = {
                {0, 2, 3, 1, 4, 1, 0, 5, 2, 4, 5, 3},
                {2, 1, 0, 3, 0, 4, 5, 1, 4, 3, 2, 5}
            };

            switch (numVertices) {
                case 4:
                    leftFace = edgeToFaceTet[0][k];
                    rightFace = edgeToFaceTet[1][k];
                    break;
                case 5:
                    leftFace = edgeToFacePyramid[0][k];
                    rightFace = edgeToFacePyramid[1][k];
                    break;
                case 6:
                    leftFace = edgeToFacePrism[0][k];
                    rightFace = edgeToFacePrism[1][k];
                    break;
                case 8:
                    leftFace = edgeToFaceHex[0][k];
                    rightFace = edgeToFaceHex[1][k];
                    break;
                default:
                    DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for numVertices = " << numVertices);
                    break;
            }
        }

    void getEdgeIndices_(int numVertices,
                         int face,
                         int node,
                         int &leftEdge,
                         int &rightEdge)
        {
            static const int faceAndNodeToLeftEdgeTet[4][4] = {
                {-1,  1,  1,  4},
                { 2, -1,  2,  3},
                { 0,  0, -1,  3},
                { 0,  0,  1, -1}
            };
            static const int faceAndNodeToRightEdgeTet[4][4] = {
                {-1,  4,  5,  5},
                { 3, -1,  5,  5},
                { 3,  4, -1,  4},
                { 2,  1,  2, -1}
            };
            static const int faceAndNodeToLeftEdgePyramid[5][5] = {
                { 3,  0,  1,  2, -1},
                { 0,  0, -1, -1,  4},
                {-1,  1,  1, -1,  5},
                {-1, -1,  2,  2,  6},
                { 3, -1, -1,  3,  4}
            };
            static const int faceAndNodeToRightEdgePyramid[5][5] = {
                { 0,  1,  2,  3, -1},
                { 4,  5, -1, -1,  5},
                {-1,  5,  6, -1,  6},
                {-1, -1,  6,  7,  7},
                { 4, -1, -1,  7,  7}
            };
            static const int faceAndNodeToLeftEdgePrism[5][6] = {
                { 0,  0,  1, -1, -1, -1},
                { 0,  0, -1,  3,  4, -1},
                {-1,  1,  1, -1,  4,  5},
                { 2, -1,  2,  3, -1,  5},
                {-1, -1, -1,  6,  6,  7}
            };
            static const int faceAndNodeToRightEdgePrism[5][6] = {
                { 2,  1,  2, -1, -1, -1},
                { 3,  4, -1,  6,  6, -1},
                {-1,  4,  5, -1,  7,  7},
                { 3, -1,  5,  8, -1,  8},
                {-1, -1, -1,  8,  7,  8}
            };
            static const int faceAndNodeToLeftEdgeHex[6][8] = {
                { 0, -1,  4, -1,  6, -1,  2, -1},
                {-1,  5, -1,  3, -1,  1, -1,  7},
                { 8,  1, -1, -1,  0, 10, -1, -1},
                {-1, -1,  2,  9, -1, -1, 11,  3},
                { 4,  8,  9,  5, -1, -1, -1, -1},
                {-1, -1, -1, -1, 10,  7,  6, 11}
            };
            static const int faceAndNodeToRightEdgeHex[6][8] = {
                { 4, -1,  2, -1,  0, -1,  6, -1},
                {-1,  1, -1,  5, -1,  7, -1,  3},
                { 0,  8, -1, -1, 10,  1, -1, -1},
                {-1, -1,  9,  3, -1, -1,  2, 11},
                { 8,  5,  4,  9, -1, -1, -1, -1},
                {-1, -1, -1, -1,  6, 10, 11,  7}
            };

            switch (numVertices) {
                case 4:
                    leftEdge = faceAndNodeToLeftEdgeTet[face][node];
                    rightEdge = faceAndNodeToRightEdgeTet[face][node];
                    break;
                case 5:
                    leftEdge = faceAndNodeToLeftEdgePyramid[face][node];
                    rightEdge = faceAndNodeToRightEdgePyramid[face][node];
                    break;
                case 6:
                    leftEdge = faceAndNodeToLeftEdgePrism[face][node];
                    rightEdge = faceAndNodeToRightEdgePrism[face][node];
                    break;
                case 8:
                    leftEdge = faceAndNodeToLeftEdgeHex[face][node];
                    rightEdge = faceAndNodeToRightEdgeHex[face][node];
                    break;
                default:
                    DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for numVertices = " << numVertices);
                    break;
            }
        }


        
    Scalar hexahedronVolume_(const WorldCoord &p0,
                             const WorldCoord &p1, 
                             const WorldCoord &p2,
                             const WorldCoord &p3,
                             const WorldCoord &p4, 
                             const WorldCoord &p5,
                             const WorldCoord &p6,
                             const WorldCoord &p7)
        {
            if (!Capabilities::IsUnstructured<Grid>::v)
                return this->cellVolume_ / 8.0;
            else
                return prismVolume_(p0,p1,p2,p4,p5,p6)
                    + prismVolume_(p0,p2,p3,p4,p6,p7);
        }


    void updateVolumes_(const Cell &cell)
        {
            switch (cell.template count<GridDim>())
            {
                case 4: // 3D, tetrahedron
                    for (int k = 0; k < this->numNodes_; k++)
                        this->subContVol_[k].volume = this->cellVolume_ / 4.0;
                    break;
                        
                case 5: // 3D, pyramid
                    this->subContVol_[0].volume = 
                        hexahedronVolume_(this->subContVol_[0].global,
                                          this->edgeCoord_[0],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[3],
                                              
                                          this->edgeCoord_[4],
                                          this->faceCoord_[1],
                                          this->elementGlobal_,
                                          this->faceCoord_[4]);
                        
                    this->subContVol_[1].volume = 
                        hexahedronVolume_(this->subContVol_[1].global,
                                          this->edgeCoord_[1],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[0],
                                              
                                          this->edgeCoord_[5],
                                          this->faceCoord_[2],
                                          this->elementGlobal_,
                                          this->faceCoord_[1]);

                    this->subContVol_[2].volume = 
                        hexahedronVolume_(this->subContVol_[2].global,
                                          this->edgeCoord_[2],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[1],
                                              
                                          this->edgeCoord_[6],
                                          this->faceCoord_[3],
                                          this->elementGlobal_,
                                          this->faceCoord_[2]);
                        
                    this->subContVol_[3].volume =
                        _hexahedronVolume(this->subContVol_[3].global,
                                          this->edgeCoord_[3],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[2],
                                              
                                          this->edgeCoord_[7],
                                          this->faceCoord_[4],
                                          this->elementGlobal_,
                                          this->faceCoord_[3]);
                        
                    this->subContVol_[4].volume = 
                        this.cellVolume_ -
                        this->subContVol_[0].volume -
                        this->subContVol_[1].volume -
                        this->subContVol_[2].volume -
                        this->subContVol_[3].volume;
                    break;
                        
                case 6: // 3D, prism
                    this->subContVol_[0].volume = 
                        hexahedronVolume_(this->subContVol_[0].global,
                                          this->edgeCoord_[0],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[2],
                                              
                                          this->edgeCoord_[3],
                                          this->faceCoord_[1],
                                          this->elementGlobal_,
                                          this->faceCoord_[3]);
                        
                    this->subContVol_[1].volume = 
                        hexahedronVolume_(this->subContVol_[1].global,
                                          this->edgeCoord_[1],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[0],
                                                                        
                                          this->edgeCoord_[4],
                                          this->faceCoord_[2],
                                          this->elementGlobal_,
                                          this->faceCoord_[1]);

                    this->subContVol_[2].volume =
                        hexahedronVolume_(this->subContVol_[2].global,
                                          this->edgeCoord_[2],
                                          this->faceCoord_[0],
                                          this->edgeCoord_[1],
                                                                        
                                          this->edgeCoord_[5],
                                          this->faceCoord_[3],
                                          this->elementGlobal_,
                                          this->faceCoord_[2]);

                    this->subContVol_[3].volume =
                        hexahedronVolume_(this->edgeCoord_[3],
                                          this->faceCoord_[1],
                                          this->elementGlobal_,
                                          this->faceCoord_[3],
                                                                        
                                          this->subContVol_[3].global,
                                          this->edgeCoord_[6],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[8]);

                    this->subContVol_[4].volume =
                        hexahedronVolume_(this->edgeCoord_[4],
                                          this->faceCoord_[2],
                                          this->elementGlobal_,
                                          this->faceCoord_[1],
                                                                        
                                          this->subContVol_[4].global,
                                          this->edgeCoord_[7],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[6]);

                    this->subContVol_[5].volume = 
                        hexahedronVolume_(this->edgeCoord_[5],
                                          this->faceCoord_[3],
                                          this->elementGlobal_,
                                          this->faceCoord_[2],

                                          this->subContVol_[5].global,
                                          this->edgeCoord_[8],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[7]);
                    break;
                        
                case 8: // 3D, hexahedron
                    this->subContVol_[0].volume = 
                        hexahedronVolume_(this->subContVol_[0].global,
                                          this->edgeCoord_[8],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[4],
                                                                        
                                          this->edgeCoord_[0],
                                          this->faceCoord_[2],
                                          this->elementGlobal_,
                                          this->faceCoord_[0]);

                    this->subContVol_[1].volume = 
                        hexahedronVolume_(this->subContVol_[1].global,
                                          this->edgeCoord_[5],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[8],

                                          this->edgeCoord_[1],
                                          this->faceCoord_[1],
                                          this->elementGlobal_,
                                          this->faceCoord_[2]);

                    this->subContVol_[2].volume =
                        hexahedronVolume_(this->subContVol_[2].global,
                                          this->edgeCoord_[4],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[9],

                                          this->edgeCoord_[2],
                                          this->faceCoord_[0],
                                          this->elementGlobal_,
                                          this->faceCoord_[3]);


                    this->subContVol_[3].volume =
                        hexahedronVolume_(this->subContVol_[3].global,
                                          this->edgeCoord_[9],
                                          this->faceCoord_[4],
                                          this->edgeCoord_[5],
                                                                        
                                          this->edgeCoord_[3],
                                          this->faceCoord_[3],
                                          this->elementGlobal_,
                                          this->faceCoord_[1]);


                    this->subContVol_[4].volume = 
                        hexahedronVolume_(this->edgeCoord_[0],
                                          this->faceCoord_[2],
                                          this->elementGlobal_,
                                          this->faceCoord_[0],
                                                                        
                                          this->subContVol_[4].global,
                                          this->edgeCoord_[10],
                                          this->faceCoord_[5],
                                          this->edgeCoord_[6]);
                        

                    this->subContVol_[5].volume = 
                        hexahedronVolume_(this->edgeCoord_[1],
                                          this->faceCoord_[1],
                                          this->elementGlobal_,
                                          this->faceCoord_[2],
                                                                        
                                          this->subContVol_[5].global,
                                          this->edgeCoord_[7],
                                          this->faceCoord_[5],
                                          this->edgeCoord_[10]);

                    this->subContVol_[6].volume = 
                        hexahedronVolume_(this->edgeCoord_[2],
                                          this->faceCoord_[0],
                                          this->elementGlobal_,
                                          this->faceCoord_[3],
                                                                        
                                          this->subContVol_[6].global,
                                          this->edgeCoord_[6],
                                          this->faceCoord_[5],
                                          this->edgeCoord_[11]);

                    this->subContVol_[7].volume = 
                        hexahedronVolume_(this->edgeCoord_[3],
                                          this->faceCoord_[3],
                                          this->elementGlobal_,
                                          this->faceCoord_[1],
                                                                        
                                          this->subContVol_[7].global,
                                          this->edgeCoord_[11],
                                          this->faceCoord_[5],
                                          this->edgeCoord_[7]);
                    break;

                default:
                    DUNE_THROW(NotImplemented, 
                               "updateVolumes_ dim = " 
                               << GridDim
                               << ", numVertices = " 
                               << cell.template count<GridDim>());
            }
        };

    void updateInteriorFaces_(const ReferenceElement &refElem)
        {
            LocalCoord ipLocal;
            WorldCoord diffVec;

            // fill sub control volume face data
            for (int scvFaceIdx = 0; scvFaceIdx < this->numEdges_; scvFaceIdx++) 
            { 
                int insideIdx = refElem.subEntity(scvFaceIdx, GridDim-1,
                                                  0, GridDim);
                int outsideIdx = refElem.subEntity(scvFaceIdx, GridDim-1,
                                                   1, GridDim);
                
                this->subContVolFace_[scvFaceIdx].insideIdx_ = insideIdx;
                this->subContVolFace_[scvFaceIdx].outsideIdx_ = outsideIdx;

                // calculate the local integration point and
                // the face normal.
                //
                // TODO: use a reference FVElementGeometry instead
                // of calculating
                int leftFace;
                int rightFace;
                getFaceIndices_(this->numNodes_, scvFaceIdx, leftFace, rightFace);
                ipLocal  = refElem.position(scvFaceIdx, GridDim-1);
                ipLocal += this->elementLocal_;
                ipLocal += refElem.position(leftFace, 1);
                ipLocal += refElem.position(rightFace, 1);
                ipLocal *= 1 / 4.0;

                this->subContVolFace_[scvFaceIdx].ipLocal_ = ipLocal;
                this->subContVolFace_[scvFaceIdx].normal_ = 
                    normalOfQuadrilateral_(this->edgeCoord_[scvFaceIdx], 
                                           this->faceCoord_[rightFace],
                                           this->elementGlobal_, 
                                           this->faceCoord_[leftFace]);
            }
        }

    void updateBoundaryFace_(IntersectionIterator &it, 
                             const ReferenceElement &refElem,
                             const CellGeometry &geometry)
        {
            LocalCoord ipLocal;
                
            // split finite element face into the sub faces of the
            // sub control volumes
            int face = it->numberInSelf();
            int numVerticesOfFace = refElem.size(face, 1, GridDim);
            for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
            {
                int nodeInElement = refElem.subEntity(face, 1,
                                                      nodeInFace, GridDim);
                int bfIndex = this->_boundaryFaceIndex(face, nodeInFace);
                    
                int leftEdge;
                int rightEdge;

                getEdgeIndices(this->numNodes_,
                               face,
                               nodeInElement,
                               leftEdge, 
                               rightEdge);

                // TODO: define the local integration point on
                // reference FVElementGeometry for each finite
                // element cell type and only transform it to
                // global coordinates here...
                ipLocal  = refElem.position(nodeInElement, GridDim);
                ipLocal += refElem.position(face, 1);
                ipLocal += refElem.position(leftEdge, GridDim-1);
                ipLocal += refElem.position(rightEdge, GridDim-1);
                ipLocal *= 1/4.0;
                this->boundaryFace_[bfIndex].ipLocal = ipLocal;

                this->boundaryFace_[bfIndex].area = 
                    quadrilateralArea_(this->subContVol_[nodeInElement].global,
                                       this->edgeCoord_[rightEdge], 
                                       this->faceCoord_[face],
                                       this->edgeCoord_[leftEdge]);
                    
                this->boundaryFace_[bfIndex].ipGlobal = geometry.global(ipLocal);
            }
        }
};
