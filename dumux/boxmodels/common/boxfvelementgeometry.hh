// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 */
#ifndef DUMUX_BOX_FV_ELEMENTGEOMETRY_HH
#define DUMUX_BOX_FV_ELEMENTGEOMETRY_HH

#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
};

/////////////////////
// HACK: Functions to initialize the subcontrol volume data
// structures of BoxFVElementGeomerty.
//
// this is a work around to be able to to use specialization for
// doing this.  it is required since it is not possible to
// specialize member functions of template classes because of some
// weird reason I didn't really get...

//! \cond INTERNAL

template <typename BoxFVElementGeometry, int dim>
class _BoxFVElemGeomHelper
{
public:
    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices)
    {
        DUNE_THROW(Dune::NotImplemented, "_BoxFVElemGeomHelper::fillSubContVolData dim = " << dim);
    };
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 1>
{
public:
    enum { dim = 1 };
    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices)
    {
        // 1D
        eg.subContVol[0].volume = 0.5*eg.elementVolume;
        eg.subContVol[1].volume = 0.5*eg.elementVolume;
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 2>
{
public:
    enum { dim = 2 };

    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices)
    {
        switch (numVertices) {
        case 3: // 2D, triangle
            eg.subContVol[0].volume = eg.elementVolume/3;
            eg.subContVol[1].volume = eg.elementVolume/3;
            eg.subContVol[2].volume = eg.elementVolume/3;
            break;
        case 4: // 2D, quadrilinear
            eg.subContVol[0].volume = eg.quadrilateralArea(eg.subContVol[0].global, eg.edgeCoord[2], eg.elementGlobal, eg.edgeCoord[0]);
            eg.subContVol[1].volume = eg.quadrilateralArea(eg.subContVol[1].global, eg.edgeCoord[1], eg.elementGlobal, eg.edgeCoord[2]);
            eg.subContVol[2].volume = eg.quadrilateralArea(eg.subContVol[2].global, eg.edgeCoord[0], eg.elementGlobal, eg.edgeCoord[3]);
            eg.subContVol[3].volume = eg.quadrilateralArea(eg.subContVol[3].global, eg.edgeCoord[3], eg.elementGlobal, eg.edgeCoord[1]);
            break;
        default:
            DUNE_THROW(Dune::NotImplemented, "_BoxFVElemGeomHelper::fillSubContVolData dim = " << dim << ", numVertices = " << numVertices);
        }
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 3>
{
public:
    enum { dim = 3 };

    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices)
    {
        switch (numVertices) {
        case 4: // 3D, tetrahedron
            for (int k = 0; k < eg.numVertices; k++)
                eg.subContVol[k].volume = eg.elementVolume / 4.0;
            break;
        case 5: // 3D, pyramid
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          eg.edgeCoord[2],
                                                          eg.faceCoord[0],
                                                          eg.edgeCoord[0],
                                                          eg.edgeCoord[4],
                                                          eg.faceCoord[3],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[1]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          eg.edgeCoord[1],
                                                          eg.faceCoord[0],
                                                          eg.edgeCoord[2],
                                                          eg.edgeCoord[5],
                                                          eg.faceCoord[2],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[3]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          eg.edgeCoord[0],
                                                          eg.faceCoord[0],
                                                          eg.edgeCoord[3],
                                                          eg.edgeCoord[6],
                                                          eg.faceCoord[1],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[4]);
            eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                          eg.edgeCoord[3],
                                                          eg.faceCoord[0],
                                                          eg.edgeCoord[1],
                                                          eg.edgeCoord[7],
                                                          eg.faceCoord[4],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[2]);
            eg.subContVol[4].volume = eg.elementVolume -
                eg.subContVol[0].volume -
                eg.subContVol[1].volume -
                eg.subContVol[2].volume -
                eg.subContVol[3].volume;
            break;
        case 6: // 3D, prism
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          eg.edgeCoord[3],
                                                          eg.faceCoord[3],
                                                          eg.edgeCoord[4],
                                                          eg.edgeCoord[0],
                                                          eg.faceCoord[0],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[1]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          eg.edgeCoord[5],
                                                          eg.faceCoord[3],
                                                          eg.edgeCoord[3],
                                                          eg.edgeCoord[1],
                                                          eg.faceCoord[2],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[0]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          eg.edgeCoord[4],
                                                          eg.faceCoord[3],
                                                          eg.edgeCoord[5],
                                                          eg.edgeCoord[2],
                                                          eg.faceCoord[1],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[2]);
            eg.subContVol[3].volume = eg.hexahedronVolume(eg.edgeCoord[0],
                                                          eg.faceCoord[0],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[1],
                                                          eg.subContVol[3].global,
                                                          eg.edgeCoord[6],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[7]);
            eg.subContVol[4].volume = eg.hexahedronVolume(eg.edgeCoord[1],
                                                          eg.faceCoord[2],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[0],
                                                          eg.subContVol[4].global,
                                                          eg.edgeCoord[8],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[6]);
            eg.subContVol[5].volume = eg.hexahedronVolume(eg.edgeCoord[2],
                                                          eg.faceCoord[1],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[2],
                                                          eg.subContVol[5].global,
                                                          eg.edgeCoord[7],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[8]);
            break;
        case 8: // 3D, hexahedron
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          eg.edgeCoord[6],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[4],
                                                          eg.edgeCoord[0],
                                                          eg.faceCoord[2],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[0]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          eg.edgeCoord[5],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[6],
                                                          eg.edgeCoord[1],
                                                          eg.faceCoord[1],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[2]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          eg.edgeCoord[4],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[7],
                                                          eg.edgeCoord[2],
                                                          eg.faceCoord[0],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[3]);
            eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                          eg.edgeCoord[7],
                                                          eg.faceCoord[4],
                                                          eg.edgeCoord[5],
                                                          eg.edgeCoord[3],
                                                          eg.faceCoord[3],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[1]);
            eg.subContVol[4].volume = eg.hexahedronVolume(eg.edgeCoord[0],
                                                          eg.faceCoord[2],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[0],
                                                          eg.subContVol[4].global,
                                                          eg.edgeCoord[10],
                                                          eg.faceCoord[5],
                                                          eg.edgeCoord[8]);
            eg.subContVol[5].volume = eg.hexahedronVolume(eg.edgeCoord[1],
                                                          eg.faceCoord[1],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[2],
                                                          eg.subContVol[5].global,
                                                          eg.edgeCoord[9],
                                                          eg.faceCoord[5],
                                                          eg.edgeCoord[10]);
            eg.subContVol[6].volume = eg.hexahedronVolume(eg.edgeCoord[2],
                                                          eg.faceCoord[0],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[3],
                                                          eg.subContVol[6].global,
                                                          eg.edgeCoord[8],
                                                          eg.faceCoord[5],
                                                          eg.edgeCoord[11]);
            eg.subContVol[7].volume = eg.hexahedronVolume(eg.edgeCoord[3],
                                                          eg.faceCoord[3],
                                                          eg.elementGlobal,
                                                          eg.faceCoord[1],
                                                          eg.subContVol[7].global,
                                                          eg.edgeCoord[11],
                                                          eg.faceCoord[5],
                                                          eg.edgeCoord[9]);
            break;
        default:
            DUNE_THROW(Dune::NotImplemented, "_BoxFVElemGeomHelper::fillSubContVolData dim = " << dim << ", numVertices = " << numVertices);
        }
    }
};

//! \endcond

// END HACK
/////////////////////

/*!
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 *
 * The box scheme is a vertex centered finite volume approach. This
 * means that each vertex corrosponds to a control volume which
 * intersects each of the vertex' neighboring elements. If only
 * looking at a single element of the primary grid (which is what this
 * class does), the element is subdivided into multiple fragments of
 * control volumes called sub-control volumes. Each of the element's
 * vertices corrosponds to exactly one sub-control volume in this
 * scenario.
 *
 * For the box methods the sub-control volumes are constructed by
 * connecting the element's center with each edge of the element.
 */
template<class TypeTag>
class BoxFVElementGeometry
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef BoxFVElementGeometry<TypeTag>   ThisType;

    /** \todo Please doc me! */
    friend class _BoxFVElemGeomHelper<ThisType, Grid::dimension>;

    typedef _BoxFVElemGeomHelper<ThisType, Grid::dimension> BoxFVElemGeomHelper;

    enum{dim = Grid::dimension};
    enum{maxNC = (dim < 3 ? 4 : 8)};
    enum{maxNE = (dim < 3 ? 4 : 12)};
    enum{maxNF = (dim < 3 ? 1 : 6)};
    enum{maxCOS = (dim < 3 ? 2 : 4)};
    enum{maxBF = (dim < 3 ? 8 : 24)};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;
  //typedef Scalar CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dim> FV;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;

    Scalar quadrilateralArea(const FV& p0, const FV& p1, const FV& p2, const FV& p3)
    {
        return 0.5*fabs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
    }

    void crossProduct(FV& c, const FV& a, const FV& b)
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    Scalar pyramidVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4)
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

    Scalar prismVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4, const FV& p5)
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

    Scalar hexahedronVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3,
                             const FV& p4, const FV& p5, const FV& p6, const FV& p7)
    {
        return
            prismVolume(p0,p1,p2,p4,p5,p6)
            + prismVolume(p0,p2,p3,p4,p6,p7);
    }

    void normalOfQuadrilateral3D(FV &normal, const FV& p0, const FV& p1, const FV& p2, const FV& p3)
    {
        FV a(p2); a -= p0;
        FV b(p3); b -= p1;
        crossProduct(normal, a, b);
        normal *= 0.5;
    }

    Scalar quadrilateralArea3D(const FV& p0, const FV& p1, const FV& p2, const FV& p3)
    {
        FV normal;
        normalOfQuadrilateral3D(normal, p0, p1, p2, p3);
        return normal.two_norm();
    }

    void getFaceIndices(int numVertices, int k, int& leftFace, int& rightFace)
    {
        static const int edgeToFaceTet[2][6] = {
            {1, 0, 3, 2, 1, 3},
            {0, 2, 0, 1, 3, 2}
        };
        static const int edgeToFacePyramid[2][8] = {
            {1, 2, 3, 4, 1, 3, 4, 2},
            {0, 0, 0, 0, 3, 2, 1, 4}
        };
        static const int edgeToFacePrism[2][9] = {
            {1, 0, 2, 0, 1, 2, 4, 4, 4},
            {0, 2, 1, 3, 3, 3, 0, 1, 2}
        };
        static const int edgeToFaceHex[2][12] = {
            {0, 2, 3, 1, 4, 1, 2, 4, 0, 5, 5, 3},
            {2, 1, 0, 3, 0, 4, 4, 3, 5, 1, 2, 5}
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
            DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry :: getFaceIndices for numVertices = " << numVertices);
            break;
        }
    }

    void getEdgeIndices(int numVertices, int face, int vert, int& leftEdge, int& rightEdge)
    {
        static const int faceAndVertexToLeftEdgeTet[4][4] = {
                { 0, 0, 2, -1},
                { 0, 0, -1, 3},
                { 1, -1, 1, 3},
                {-1, 2, 2, 4}
        };
        static const int faceAndVertexToRightEdgeTet[4][4] = {
                { 1, 2, 1, -1},
                { 3, 4, -1, 4},
                { 3, -1, 5, 5},
                {-1, 4, 5, 5}
        };
        static const int faceAndVertexToLeftEdgePyramid[5][5] = {
            { 0, 2, 3, 1, -1},
            { 0, -1, 0, -1, 4},
            {-1, 1, -1, 1, 5},
            { 2, 2, -1, -1, 4},
            {-1, -1, 3, 3, 7}
        };
        static const int faceAndVertexToRightEdgePyramid[5][5] = {
            { 2, 1, 0, 3, -1},
            { 4, -1, 6, -1, 6},
            {-1, 5, -1, 7, 7},
            { 4, 5, -1, -1, 5},
            {-1, -1, 6, 7, 6}
        };
        static const int faceAndVertexToLeftEdgePrism[5][6] = {
            { 3, 3, -1, 0, 1, -1},
            { 4, -1, 4, 0, -1, 2},
            {-1, 5, 5, -1, 1, 2},
            { 3, 3, 5, -1, -1, -1},
            {-1, -1, -1, 6, 6, 8}
        };
        static const int faceAndVertexToRightEdgePrism[5][6] = {
            { 0, 1, -1, 6, 6, -1},
            { 0, -1, 2, 7, -1, 7},
            {-1, 1, 2, -1, 8, 8},
            { 4, 5, 4, -1, -1, -1},
            {-1, -1, -1, 7, 8, 7}
        };
        static const int faceAndVertexToLeftEdgeHex[6][8] = {
            { 0, -1, 4, -1, 8, -1, 2, -1},
            {-1, 5, -1, 3, -1, 1, -1, 9},
            { 6, 1, -1, -1, 0, 10, -1, -1},
            {-1, -1, 2, 7, -1, -1, 11, 3},
            { 4, 6, 7, 5, -1, -1, -1, -1},
            {-1, -1, -1, -1, 10, 9, 8, 11}
        };
        static const int faceAndVertexToRightEdgeHex[6][8] = {
            { 4, -1, 2, -1, 0, -1, 8, -1},
            {-1, 1, -1, 5, -1, 9, -1, 3},
            { 0, 6, -1, -1, 10, 1, -1, -1},
            {-1, -1, 7, 3, -1, -1, 2, 11},
            { 6, 5, 4, 7, -1, -1, -1, -1},
            {-1, -1, -1, -1, 8, 10, 11, 9}
        };

        switch (numVertices) {
        case 4:
            leftEdge = faceAndVertexToLeftEdgeTet[face][vert];
            rightEdge = faceAndVertexToRightEdgeTet[face][vert];
            break;
        case 5:
            leftEdge = faceAndVertexToLeftEdgePyramid[face][vert];
            rightEdge = faceAndVertexToRightEdgePyramid[face][vert];
            break;
        case 6:
            leftEdge = faceAndVertexToLeftEdgePrism[face][vert];
            rightEdge = faceAndVertexToRightEdgePrism[face][vert];
            break;
        case 8:
            leftEdge = faceAndVertexToLeftEdgeHex[face][vert];
            rightEdge = faceAndVertexToRightEdgeHex[face][vert];
            break;
        default:
            DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry :: getFaceIndices for numVertices = " << numVertices);
            break;
        }
    }

public:
    int boundaryFaceIndex(int face, int vertInFace) const
    {
        return (face*maxCOS + vertInFace);
    }

    struct SubControlVolume //! FV intersected with element
    {
        FV local; //!< local vert position
        FV global; //!< global vert position
        FV localCenter; //!< local position of scv center
        Scalar volume; //!< volume of scv
        Dune::FieldVector<FV, maxNC> grad; //! derivative of shape function associated with the sub control volume
        Dune::FieldVector<FV, maxNC> gradCenter; //! derivative of shape function at the center of the sub control volume
        Dune::FieldVector<Scalar, maxNC> shapeValue; //! value of shape function associated with the sub control volume
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
        FV ipLocal; //!< integration point in local coords
        FV ipGlobal; //!< integration point in global coords
        FV normal; //!< normal on face at ip pointing to CV j with length equal to |scvf|
        Dune::FieldVector<FV, maxNC> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNC> shapeValue; //!< value of shape functions at ip
    };

    struct BoundaryFace //! boundary face of a sub control volume
    {
        FV ipLocal; //!< integration point in local coords
        FV ipGlobal; //!< integration point in global coords
        Scalar area; //!< area of boundary face
        FV normal; //!< normal on face at ip pointing to CV j with length equal to |scvf|
        Dune::FieldVector<FV, maxNC> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNC> shapeValue; //!< value of shape functions at ip
    };

    FV elementLocal; //!< local coordinate of element center
    FV elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[maxNC]; //!< data of the sub control volumes
    SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    FV edgeCoord[maxNE]; //!< global coordinates of the edge centers
    FV faceCoord[maxNF]; //!< global coordinates of the face centers
    int numVertices; //!< number of verts
    int numEdges; //!< number of edges
    int numFaces; //!< number of faces (0 in < 3D)

    const LocalFiniteElementCache feCache_;
    bool computeGradientAtScvCenters;

    BoxFVElementGeometry()
    {
        computeGradientAtScvCenters = false;
    }

    BoxFVElementGeometry(bool computeGradientAtCenters)
    {
        computeGradientAtScvCenters = computeGradientAtCenters;
    };

    void update(const GridView& gridView, const Element& e)
    {
        const Geometry& geometry = e.geometry();
        Dune::GeometryType gt = geometry.type();

        const typename Dune::GenericReferenceElementContainer<CoordScalar,dim>::value_type&
            referenceElement = Dune::GenericReferenceElements<CoordScalar,dim>::general(gt);

        const LocalFiniteElement
            &localFiniteElement = feCache_.get(gt);

        elementVolume = geometry.volume();
        elementLocal = referenceElement.position(0,0);
        elementGlobal = geometry.global(elementLocal);

        numVertices = referenceElement.size(dim);
        numEdges = referenceElement.size(dim-1);
        numFaces = (dim<3)?0:referenceElement.size(1);

        // corners:
        for (int vert = 0; vert < numVertices; vert++) {
            subContVol[vert].local = referenceElement.position(vert, dim);
            subContVol[vert].global = geometry.global(subContVol[vert].local);
            subContVol[vert].inner = true;
        }

        // edges:
        for (int edge = 0; edge < numEdges; edge++) {
            edgeCoord[edge] = geometry.global(referenceElement.position(edge, dim-1));
        }

        // faces:
        for (int face = 0; face < numFaces; face++) {
            faceCoord[face] = geometry.global(referenceElement.position(face, 1));
        }

        // fill sub control volume data use specialization for this
        // \todo maybe it would be a good idea to migrate everything
        // which is dependend of the grid's dimension to
        // _BoxFVElemGeomHelper in order to benefit from more aggressive
        // compiler optimizations...
        BoxFVElemGeomHelper::fillSubContVolData(*this, numVertices);

        // fill sub control volume face data:
        for (int k = 0; k < numEdges; k++) { // begin loop over edges / sub control volume faces
            int i = referenceElement.subEntity(k, dim-1, 0, dim);
            int j = referenceElement.subEntity(k, dim-1, 1, dim);
            if (numEdges == 4 && (i == 2 || j == 2))
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
                ipLocal = referenceElement.position(k, dim-1) + elementLocal;
                ipLocal *= 0.5;
                subContVolFace[k].ipLocal = ipLocal;
                diffVec = elementGlobal - edgeCoord[k];
                subContVolFace[k].normal[0] = diffVec[1];
                subContVolFace[k].normal[1] = -diffVec[0];

                diffVec = subContVol[j].global;
                diffVec -= subContVol[i].global;
                // make sure the normal points to the right direction
                if (subContVolFace[k].normal * diffVec < 0)
                    subContVolFace[k].normal *= -1;

            }
            else if (dim==3) {
                int leftFace;
                int rightFace;
                getFaceIndices(numVertices, k, leftFace, rightFace);
                ipLocal = referenceElement.position(k, dim-1) + elementLocal
                    + referenceElement.position(leftFace, 1)
                    + referenceElement.position(rightFace, 1);
                ipLocal *= 0.25;
                subContVolFace[k].ipLocal = ipLocal;
                normalOfQuadrilateral3D(subContVolFace[k].normal,
                                        edgeCoord[k], faceCoord[rightFace],
                                        elementGlobal, faceCoord[leftFace]);
            }

            // get the global integration point and the Jacobian inverse
            subContVolFace[k].ipGlobal = geometry.global(ipLocal);
            Dune::FieldMatrix<CoordScalar,dim,dim> jacInvT = geometry.jacobianInverseTransposed(ipLocal);


            //              std::cout << "SCV Face " << k << ", i = " << i << ", j = " << j
            //                          << ", ipLocal = " << ipLocal << ", ipGlobal = " << subContVolFace[k].ipGlobal << ", normal = " << subContVolFace[k].normal
            //                          << std::endl;

            // calculate the shape function gradients
            //typedef Dune::FieldVector< Dune::FieldVector< CoordScalar, dim >, 1 > ShapeJacobian;
            typedef Dune::FieldVector< CoordScalar, 1 > ShapeValue;
            std::vector<ShapeJacobian> localJac;
            std::vector<ShapeValue>    shapeVal;
            localFiniteElement.localBasis().evaluateJacobian(subContVolFace[k].ipLocal, localJac);
            localFiniteElement.localBasis().evaluateFunction(subContVolFace[k].ipLocal, shapeVal);
            for (int vert = 0; vert < numVertices; vert++) {
                jacInvT.mv(localJac[vert][0], subContVolFace[k].grad[vert]);
                subContVolFace[k].shapeValue[vert] = Scalar(shapeVal[vert]);
            }
        } // end loop over edges / sub control volume faces

        // fill boundary face data:
        IntersectionIterator endit = gridView.iend(e);
        for (IntersectionIterator it = gridView.ibegin(e); it != endit; ++it)
            if (it->boundary())
            {
                int face = it->indexInInside();
                int numVerticesOfFace = referenceElement.size(face, 1, dim);
                for (int vertInFace = 0; vertInFace < numVerticesOfFace; vertInFace++)
                {
                    int vertInElement = referenceElement.subEntity(face, 1, vertInFace, dim);
                    int bfIdx = boundaryFaceIndex(face, vertInFace);
                    subContVol[vertInElement].inner = false;
                    switch ((short) dim) {
                    case 1:
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim);
                        boundaryFace[bfIdx].area = 1.0;
                        break;
                    case 2:
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim)
                            + referenceElement.position(face, 1);
                        boundaryFace[bfIdx].ipLocal *= 0.5;
                        boundaryFace[bfIdx].area = 0.5*it->geometry().volume();
                        break;
                    case 3:
                        int leftEdge;
                        int rightEdge;
                        getEdgeIndices(numVertices, face, vertInElement, leftEdge, rightEdge);
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim)
                            + referenceElement.position(face, 1)
                            + referenceElement.position(leftEdge, dim-1)
                            + referenceElement.position(rightEdge, dim-1);
                        boundaryFace[bfIdx].ipLocal *= 0.25;
                        boundaryFace[bfIdx].area = quadrilateralArea3D(subContVol[vertInElement].global,
                                                                       edgeCoord[rightEdge], faceCoord[face], edgeCoord[leftEdge]);
                        break;
                    default:
                        DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry for dim = " << dim);
                    }
                    boundaryFace[bfIdx].ipGlobal = geometry.global(boundaryFace[bfIdx].ipLocal);

                    // ASSUME constant normal
                    Dune::FieldVector<CoordScalar, dim-1> localDimM1(0);
                    boundaryFace[bfIdx].normal = it->unitOuterNormal(localDimM1);
                    boundaryFace[bfIdx].normal *= boundaryFace[bfIdx].area;

                    typedef Dune::FieldVector< CoordScalar, 1 > ShapeValue;
                    std::vector<ShapeJacobian> localJac;
                    std::vector<ShapeValue>    shapeVal;
                    localFiniteElement.localBasis().evaluateJacobian(boundaryFace[bfIdx].ipLocal, localJac);
                    localFiniteElement.localBasis().evaluateFunction(boundaryFace[bfIdx].ipLocal, shapeVal);

                    Dune::FieldMatrix<CoordScalar,dim,dim> jacInvT = geometry.jacobianInverseTransposed(boundaryFace[bfIdx].ipLocal);
                    for (int vert = 0; vert < numVertices; vert++)
                    {
                        jacInvT.mv(localJac[vert][0], boundaryFace[bfIdx].grad[vert]);
                        boundaryFace[bfIdx].shapeValue[vert] = Scalar(shapeVal[vert]);
                    }

                    //                    std::cout << "boundary face " << face << ", vert = " << vertInElement << ", ipLocal = "
                    //                        << boundaryFace[bfIdx].ipLocal << ", ipGlobal = " << boundaryFace[bfIdx].ipGlobal
                    //                        << ", area = " << boundaryFace[bfIdx].area << std::endl;

                }
            }


            // calculate gradients at the center of the scv
            for (int scvIdx = 0; scvIdx < numVertices; scvIdx++)
                if (dim == 2)
                {
//                    if (!subContVol[scvIdx].inner)
                    {
                        switch (scvIdx)
                        {
                        case 0:
                            if (numVertices == 4) {
                                subContVol[scvIdx].localCenter[0] = 0.25;
                                subContVol[scvIdx].localCenter[1] = 0.25;
                            }
                            else {
                                subContVol[scvIdx].localCenter[0] = 1.0/6.0;
                                subContVol[scvIdx].localCenter[1] = 1.0/6.0;
                            }
                            break;
                        case 1:
                            if (numVertices == 4) {
                                subContVol[scvIdx].localCenter[0] = 0.75;
                                subContVol[scvIdx].localCenter[1] = 0.25;
                            }
                            else {
                                subContVol[scvIdx].localCenter[0] = 4.0/6.0;
                                subContVol[scvIdx].localCenter[1] = 1.0/6.0;
                            }
                            break;
                        case 2:
                            if (numVertices == 4) {
                                subContVol[scvIdx].localCenter[0] = 0.25;
                                subContVol[scvIdx].localCenter[1] = 0.75;
                            }
                            else {
                                subContVol[scvIdx].localCenter[0] = 1.0/6.0;
                                subContVol[scvIdx].localCenter[1] = 4.0/6.0;
                            }
                            break;
                        case 3:
                            subContVol[scvIdx].localCenter[0] = 0.75;
                            subContVol[scvIdx].localCenter[1] = 0.75;
                            break;
                        }
                    }

                    std::vector<ShapeJacobian> localJac;
                    localFiniteElement.localBasis().evaluateJacobian(subContVol[scvIdx].localCenter, localJac);

                    Dune::FieldMatrix<CoordScalar,dim,dim> jacInvT =
                            geometry.jacobianInverseTransposed(subContVol[scvIdx].localCenter);
                    for (int vert = 0; vert < numVertices; vert++)
                        jacInvT.mv(localJac[vert][0], subContVol[scvIdx].gradCenter[vert]);
                }
    }
};

}

#endif

