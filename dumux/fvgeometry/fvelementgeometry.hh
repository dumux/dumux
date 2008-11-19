// $Id$ 

#ifndef DUNE_FVELEMENTGEOMETRY_HH
#define DUNE_FVELEMENTGEOMETRY_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/utility/intersectiongetter.hh>

namespace Dune
{
    /////////////////////
    // HACK: Functions to initialize the subcontrol volume data
    // structures of FVElementGeomerty.
    //
    // this is a work around to be able to to use specialization for
    // doing this.  it is required since it is not possible to
    // specialize member functions of template classes because of some
    // weird reason I didn't really get...
    template <typename FVElementGeometry, int dim>
    class _FVElemGeomHelper
    {
    public:
        static void fillSubContVolData(FVElementGeometry &eg, int nNodes)
            {
                DUNE_THROW(NotImplemented, "_FVElemGeomHelper::fillSubContVolData dim = " << dim);
            };
    };

    template <typename FVElementGeometry>
    class _FVElemGeomHelper<FVElementGeometry, 1>
    {
    public:
        enum { dim = 1 };
        static void fillSubContVolData(FVElementGeometry &eg, int nNodes)
            {
                // 1D
                eg.subContVol[0].volume = 0.5*eg.cellVolume;
                eg.subContVol[1].volume = 0.5*eg.cellVolume;
            }
    };

    template <typename FVElementGeometry>
    class _FVElemGeomHelper<FVElementGeometry, 2>
    {
    public:
        enum { dim = 2 };

        static void fillSubContVolData(FVElementGeometry &eg, int nNodes)
            {
                switch (nNodes) {
                    case 3: // 2D, triangle
                        eg.subContVol[0].volume = eg.quadrilateralArea(eg.subContVol[0].global, eg.edgeCoord[2], eg.cellGlobal, eg.edgeCoord[1]);
                        eg.subContVol[1].volume = eg.quadrilateralArea(eg.subContVol[1].global, eg.edgeCoord[0], eg.cellGlobal, eg.edgeCoord[2]);
                        eg.subContVol[2].volume = eg.quadrilateralArea(eg.subContVol[2].global, eg.edgeCoord[1], eg.cellGlobal, eg.edgeCoord[0]);
                        break;
                    case 4: // 2D, quadrilinear
                        eg.subContVol[0].volume = eg.quadrilateralArea(eg.subContVol[0].global, eg.edgeCoord[2], eg.cellGlobal, eg.edgeCoord[0]);
                        eg.subContVol[1].volume = eg.quadrilateralArea(eg.subContVol[1].global, eg.edgeCoord[1], eg.cellGlobal, eg.edgeCoord[2]);
                        eg.subContVol[2].volume = eg.quadrilateralArea(eg.subContVol[2].global, eg.edgeCoord[0], eg.cellGlobal, eg.edgeCoord[3]);
                        eg.subContVol[3].volume = eg.quadrilateralArea(eg.subContVol[3].global, eg.edgeCoord[3], eg.cellGlobal, eg.edgeCoord[1]);
                        break;
                    default:
                        DUNE_THROW(NotImplemented, "_FVElemGeomHelper::fillSubContVolData dim = " << dim << ", nNodes = " << nNodes);
                }
            }
    };

    template <typename FVElementGeometry>
    class _FVElemGeomHelper<FVElementGeometry, 3>
    {
    public:
        enum { dim = 3 };

        static void fillSubContVolData(FVElementGeometry &eg, int nNodes)
            {
                switch (nNodes) {
                    case 4: // 3D, tetrahedron
                        for (int k = 0; k < eg.nNodes; k++)
                            eg.subContVol[k].volume = 0.25*eg.cellVolume;
                        break;
                    case 5: // 3D, pyramid
                        eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                                      eg.edgeCoord[0],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[3],
                                                                      eg.edgeCoord[4],
                                                                      eg.faceCoord[1],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[4]);
                        eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                                      eg.edgeCoord[1],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[0],
                                                                      eg.edgeCoord[5],
                                                                      eg.faceCoord[2],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[1]);
                        eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                                      eg.edgeCoord[2],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[1],
                                                                      eg.edgeCoord[6],
                                                                      eg.faceCoord[3],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[2]);
                        eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                                      eg.edgeCoord[3],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[2],
                                                                      eg.edgeCoord[7],
                                                                      eg.faceCoord[4],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[3]);
                        eg.subContVol[4].volume = eg.cellVolume -
                            eg.subContVol[0].volume -
                            eg.subContVol[1].volume -
                            eg.subContVol[2].volume -
                            eg.subContVol[3].volume;
                        break;
                    case 6: // 3D, prism
                        eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                                      eg.edgeCoord[0],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[2],
                                                                      eg.edgeCoord[3],
                                                                      eg.faceCoord[1],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[3]);
                        eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                                      eg.edgeCoord[1],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[0],
                                                                      eg.edgeCoord[4],
                                                                      eg.faceCoord[2],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[1]);
                        eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                                      eg.edgeCoord[2],
                                                                      eg.faceCoord[0],
                                                                      eg.edgeCoord[1],
                                                                      eg.edgeCoord[5],
                                                                      eg.faceCoord[3],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[2]);
                        eg.subContVol[3].volume = eg.hexahedronVolume(eg.edgeCoord[3],
                                                                      eg.faceCoord[1],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[3],
                                                                      eg.subContVol[3].global,
                                                                      eg.edgeCoord[6],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[8]);
                        eg.subContVol[4].volume = eg.hexahedronVolume(eg.edgeCoord[4],
                                                                      eg.faceCoord[2],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[1],
                                                                      eg.subContVol[4].global,
                                                                      eg.edgeCoord[7],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[6]);
                        eg.subContVol[5].volume = eg.hexahedronVolume(eg.edgeCoord[5],
                                                                      eg.faceCoord[3],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[2],
                                                                      eg.subContVol[5].global,
                                                                      eg.edgeCoord[8],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[7]);
                        break;
                    case 8: // 3D, hexahedron
                        eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                                      eg.edgeCoord[8],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[4],
                                                                      eg.edgeCoord[0],
                                                                      eg.faceCoord[2],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[0]);
                        eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                                      eg.edgeCoord[5],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[8],
                                                                      eg.edgeCoord[1],
                                                                      eg.faceCoord[1],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[2]);
                        eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                                      eg.edgeCoord[4],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[9],
                                                                      eg.edgeCoord[2],
                                                                      eg.faceCoord[0],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[3]);
                        eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                                      eg.edgeCoord[9],
                                                                      eg.faceCoord[4],
                                                                      eg.edgeCoord[5],
                                                                      eg.edgeCoord[3],
                                                                      eg.faceCoord[3],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[1]);
                        eg.subContVol[4].volume = eg.hexahedronVolume(eg.edgeCoord[0],
                                                                      eg.faceCoord[2],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[0],
                                                                      eg.subContVol[4].global,
                                                                      eg.edgeCoord[10],
                                                                      eg.faceCoord[5],
                                                                      eg.edgeCoord[6]);
                        eg.subContVol[5].volume = eg.hexahedronVolume(eg.edgeCoord[1],
                                                                      eg.faceCoord[1],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[2],
                                                                      eg.subContVol[5].global,
                                                                      eg.edgeCoord[7],
                                                                      eg.faceCoord[5],
                                                                      eg.edgeCoord[10]);
                        eg.subContVol[6].volume = eg.hexahedronVolume(eg.edgeCoord[2],
                                                                      eg.faceCoord[0],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[3],
                                                                      eg.subContVol[6].global,
                                                                      eg.edgeCoord[6],
                                                                      eg.faceCoord[5],
                                                                      eg.edgeCoord[11]);
                        eg.subContVol[7].volume = eg.hexahedronVolume(eg.edgeCoord[3],
                                                                      eg.faceCoord[3],
                                                                      eg.cellGlobal,
                                                                      eg.faceCoord[1],
                                                                      eg.subContVol[7].global,
                                                                      eg.edgeCoord[11],
                                                                      eg.faceCoord[5],
                                                                      eg.edgeCoord[7]);
                        break;
                    default:
                        DUNE_THROW(NotImplemented, "_FVElemGeomHelper::fillSubContVolData dim = " << dim << ", nNodes = " << nNodes);
                }
            }
    };
    // END HACK
    /////////////////////

    template<class Grid>
    class FVElementGeometry
    {
        friend class _FVElemGeomHelper<FVElementGeometry<Grid>, Grid::dimension>;

        typedef _FVElemGeomHelper<FVElementGeometry<Grid>, Grid::dimension> FVElemGeomHelper;

        enum{dim = Grid::dimension};
        enum{maxNC = (dim < 3 ? 4 : 8)};
        enum{maxNE = (dim < 3 ? 4 : 12)};
        enum{maxNF = (dim < 3 ? 0 : 6)};
        enum{maxCOS = (dim < 3 ? 2 : 4)};
        enum{maxBF = (dim < 3 ? 8 : 24)};
        typedef typename Grid::ctype DT;
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Entity::Geometry Geometry;
        typedef FieldVector<DT,dim> FV;
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        DT quadrilateralArea(const FV& p0, const FV& p1, const FV& p2, const FV& p3)
            {
                return 0.5*fabs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
            }

        void crossProduct(FV& c, const FV& a, const FV& b)
            {
                c[0] = a[1]*b[2] - a[2]*b[1];
                c[1] = a[2]*b[0] - a[0]*b[2];
                c[2] = a[0]*b[1] - a[1]*b[0];
            }

        DT pyramidVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4)
            {
/*
                FV a = p2 - p0;
                FV b = p3 - p1;
                FV h = p4 - p0;
                FV n = crossProduct(a, b);
                return 1.0/6.0*(n*h);
*/

                FV a(p2); a -= p0;
                FV b(p3); b -= p1;
                
                FV n;
                crossProduct(n, a, b);
                
                a = p4; a -= p0;

                return 1.0/6.0*(n*a);
            }

        DT prismVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4, const FV& p5)
            {
/*
                FV a = p4 - p0;
                FV b = p1 - p3;
                FV c = p1 - p0;
                FV d = p2 - p0;
                FV e = p5 - p0;
                FV m = crossProduct(a, b);
                FV n = m + crossProduct(c, d);

                return fabs(1.0/6.0*(n*e));
*/

                FV a(p4); a -= p0;
                FV b(p1); b -= p3;
                FV m;
                crossProduct(m, a, b);

                a = p1; a -= p0;
                b = p2; b -= p0;
                FV n;
                crossProduct(n, a, b);
                n += m;

                a = p5; a -= p0;

                return fabs(1.0/6.0*(n*a));
            }

        DT hexahedronVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3,
                             const FV& p4, const FV& p5, const FV& p6, const FV& p7)
            {
                return prismVolume(p0,p1,p2,p4,p5,p6)
                       + prismVolume(p0,p2,p3,p4,p6,p7);
            }

        void normalOfQuadrilateral3D(FV &normal, const FV& p0, const FV& p1, const FV& p2, const FV& p3)
            {
                FV a(p2); a -= p0;
                FV b(p3); b -= p1;
                crossProduct(normal, a, b);
                normal *= 0.5;
            }

        DT quadrilateralArea3D(const FV& p0, const FV& p1, const FV& p2, const FV& p3)
            {
                FV normal;
                normalOfQuadrilateral3D(normal, p0, p1, p2, p3);
                return normal.two_norm();
            }

        void getFaceIndices(int nNodes, int k, int& leftFace, int& rightFace)
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

                switch (nNodes) {
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
                        DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for nNodes = " << nNodes);
                        break;
                }
            }

        void getEdgeIndices(int nNodes, int face, int node, int& leftEdge, int& rightEdge)
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

                switch (nNodes) {
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
                        DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for nNodes = " << nNodes);
                        break;
                }
            }

    public:
        int boundaryFaceIndex(int face, int nodeInFace) const
            {
                return (face*maxCOS + nodeInFace);
            }

        struct SubControlVolume //!< FV intersected with element
        {
            FV local; //!< local node position
            FV global; //!< global node position
            DT volume; //!< volume of scv
        };

        struct SubControlVolumeFace
        {
            int i,j; //!< scvf seperates corner i and j of elem
            FV ipLocal; //!< integration point in local coords
            FV ipGlobal; //!< integration point in global coords
            FV normal; //!< normal on face at ip pointing to CV j with length equal to |scvf|
            FieldVector<FV, maxNC> grad; //!< derivatives of shape functions at ip
        };

        struct BoundaryFace {
            FV ipLocal; //!< integration point in local coords
            FV ipGlobal; //!< integration point in global coords
            DT area; //!< area of boundary face
        };

        FV cellLocal; //!< local coordinate of cell center
        FV cellGlobal; //!< global coordinate of cell center
        DT cellVolume; //!< cell volume
        SubControlVolume subContVol[maxNC]; //!< data of the sub control volumes
        SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
        BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
        FV edgeCoord[maxNE]; //!< global coordinates of the edge centers
        FV faceCoord[maxNF]; //!< global coordinates of the face centers
        int nNodes; //!< number of nodes
        int nEdges; //!< number of edges
        int nFaces; //!< number of faces (0 in < 3D)

        FVElementGeometry<Grid>()
            {}

        void update(const Entity& e)
            {
                const Geometry& geometry = e.geometry();
                GeometryType gt = geometry.type();

                const typename ReferenceElementContainer<DT,dim>::value_type&
                    referenceElement = ReferenceElements<DT,dim>::general(gt);

                const typename LagrangeShapeFunctionSetContainer<DT,DT,dim>::value_type&
                    sfs=LagrangeShapeFunctions<DT,DT,dim>::general(gt, 1);

                cellVolume = geometry.volume();
                cellLocal = referenceElement.position(0,0);
                cellGlobal = geometry.global(cellLocal);

                nNodes = referenceElement.size(dim);
                nEdges = referenceElement.size(dim-1);
                nFaces = (dim<3)?0:referenceElement.size(1);

                // corners:
                for (int node = 0; node < nNodes; node++) {
                    subContVol[node].local  = referenceElement.position(node, dim);
                    subContVol[node].global = geometry.global(subContVol[node].local);
                }

                // edges:
                for (int edge = 0; edge < nEdges; edge++) {
                    edgeCoord[edge] = geometry.global(referenceElement.position(edge, dim-1));
                }

                // faces:
                for (int face = 0; face < nFaces; face++) {
                    faceCoord[face] = geometry.global(referenceElement.position(face, 1));
                }

                // fill sub control volume data use specialization for this
                // \todo maybe it would be a good idea to migrate everything
                // which is dependend of the grid's dimension to
                // _FVElemGeomHelper in order to benefit from more aggressive
                // compiler optimizations...
                FVElemGeomHelper::fillSubContVolData(*this, nNodes);

                // fill sub control volume face data:
                for (int k = 0; k < nEdges; k++) { // begin loop over edges / sub control volume faces
                    int i = referenceElement.subEntity(k, dim-1, 0, dim);
                    int j = referenceElement.subEntity(k, dim-1, 1, dim);
                    if (nEdges == 4 && (i == 2 || j == 2))
                        std::swap(i, j);
                    subContVolFace[k].i = i;
                    subContVolFace[k].j = j;

                    // calculate the local integration point and
                    // the face normal. note that since dim is a
                    // constant which is known at compile time
                    // the compiler can optimize away all if
                    // cases which don't apply.
                    FV ipLocal;
                    FV diffVec;
                    if (dim==1) {
                        subContVolFace[k].ipLocal = 0.5;
                        subContVolFace[k].normal = 1.0;
                        ipLocal = subContVolFace[k].ipLocal;
                    }
                    else if (dim==2) {
                        ipLocal = referenceElement.position(k, dim-1) + cellLocal;
                        ipLocal *= 0.5;
                        subContVolFace[k].ipLocal = ipLocal;
                        diffVec = cellGlobal - edgeCoord[k];
                        subContVolFace[k].normal[0] = diffVec[1];
                        subContVolFace[k].normal[1] = -diffVec[0];
                    }
                    else if (dim==3) {
                        int leftFace;
                        int rightFace;
                        getFaceIndices(nNodes, k, leftFace, rightFace);
                        ipLocal = referenceElement.position(k, dim-1) + cellLocal
                            + referenceElement.position(leftFace, 1)
                            + referenceElement.position(rightFace, 1);
                        ipLocal *= 0.25;
                        subContVolFace[k].ipLocal = ipLocal;
                        normalOfQuadrilateral3D(subContVolFace[k].normal,
                                                edgeCoord[k], faceCoord[rightFace],
                                                cellGlobal, faceCoord[leftFace]);
                    }

                    // get the global integration point and the Jacobian inverse
                    subContVolFace[k].ipGlobal = geometry.global(ipLocal);
                    FieldMatrix<DT,dim,dim> jacInvT = geometry.jacobianInverseTransposed(ipLocal);


//              std::cout << "SCV Face " << k << ", i = " << i << ", j = " << j
//                          << ", ipLocal = " << ipLocal << ", ipGlobal = " << subContVolFace[k].ipGlobal << ", normal = " << subContVolFace[k].normal
//                          << std::endl;

                    // calculate the shape function gradients
                    for (int node = 0; node < nNodes; node++) {
                        FV grad(0),temp;
                        for (int l = 0; l < dim; l++)
                            temp[l] = sfs[node].evaluateDerivative(0, l, subContVolFace[k].ipLocal);
                        jacInvT.umv(temp, grad);
                        subContVolFace[k].grad[node] = grad;
                    }
                } // end loop over edges / sub control volume faces

                // fill boundary face data:
                IntersectionIterator endit = e.ileafend();
                for (IntersectionIterator it = e.ileafbegin(); it != endit; ++it)
                    if (it->boundary())
                    {
                        int face = it->numberInSelf();
                        int nNodesOfFace = referenceElement.size(face, 1, dim);
                        for (int nodeInFace = 0; nodeInFace < nNodesOfFace; nodeInFace++)
                        {
                            int nodeInElement = referenceElement.subEntity(face, 1, nodeInFace, dim);
                            int bfIndex = boundaryFaceIndex(face, nodeInFace);
                            switch (dim) {
                                case 1:
                                    boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim);
                                    boundaryFace[bfIndex].area = 1.0;
                                    break;
                                case 2:
                                    boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim)
                                        + referenceElement.position(face, 1);
                                    boundaryFace[bfIndex].ipLocal *= 0.5;
                                    boundaryFace[bfIndex].area = 0.5*it->intersectionGlobal().volume();
                                    break;
                                case 3:
                                    int leftEdge;
                                    int rightEdge;
                                    getEdgeIndices(nNodes, face, nodeInElement, leftEdge, rightEdge);
                                    boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim)
                                        + referenceElement.position(face, 1)
                                        + referenceElement.position(leftEdge, dim-1)
                                        + referenceElement.position(rightEdge, dim-1);
                                    boundaryFace[bfIndex].ipLocal *= 0.25;
                                    boundaryFace[bfIndex].area = quadrilateralArea3D(subContVol[nodeInElement].global,
                                                                                     edgeCoord[rightEdge], faceCoord[face], edgeCoord[leftEdge]);
                                    break;
                                default:
                                    DUNE_THROW(NotImplemented, "FVElementGeometry for dim = " << dim);
                            }
                            boundaryFace[bfIndex].ipGlobal = geometry.global(boundaryFace[bfIndex].ipLocal);

                            //                    std::cout << "boundary face " << face << ", node = " << nodeInElement << ", ipLocal = "
                            //                        << boundaryFace[bfIndex].ipLocal << ", ipGlobal = " << boundaryFace[bfIndex].ipGlobal
                            //                        << ", area = " << boundaryFace[bfIndex].area << std::endl;

                        }
                    }
            }

    };

}


#endif

