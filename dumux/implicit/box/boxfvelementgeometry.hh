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
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 */
#ifndef DUMUX_BOX_FV_ELEMENTGEOMETRY_HH
#define DUMUX_BOX_FV_ELEMENTGEOMETRY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/propertysystem.hh>
#include "boxproperties.hh"

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(ImplicitUseTwoPointFlux);
}

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
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 1>
{
public:
    enum { dim = 1 };
    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices, GlobalPosition *edgeCoord, GlobalPosition *faceCoord)
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

    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices, GlobalPosition *edgeCoord, GlobalPosition *faceCoord)
    {
        switch (numVertices) {
        case 3: // 2D, triangle
            eg.subContVol[0].volume = eg.elementVolume/3;
            eg.subContVol[1].volume = eg.elementVolume/3;
            eg.subContVol[2].volume = eg.elementVolume/3;
            break;
        case 4: // 2D, quadrilinear
            eg.subContVol[0].volume = eg.quadrilateralArea(eg.subContVol[0].global, edgeCoord[2], eg.elementGlobal, edgeCoord[0]);
            eg.subContVol[1].volume = eg.quadrilateralArea(eg.subContVol[1].global, edgeCoord[1], eg.elementGlobal, edgeCoord[2]);
            eg.subContVol[2].volume = eg.quadrilateralArea(eg.subContVol[2].global, edgeCoord[0], eg.elementGlobal, edgeCoord[3]);
            eg.subContVol[3].volume = eg.quadrilateralArea(eg.subContVol[3].global, edgeCoord[3], eg.elementGlobal, edgeCoord[1]);
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

    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &eg, int numVertices, GlobalPosition *edgeCoord, GlobalPosition *faceCoord)
    {
        switch (numVertices) {
        case 4: // 3D, tetrahedron
            for (int k = 0; k < eg.numScv; k++)
                eg.subContVol[k].volume = eg.elementVolume / 4.0;
            break;
        case 5: // 3D, pyramid
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          edgeCoord[2],
                                                          faceCoord[0],
                                                          edgeCoord[0],
                                                          edgeCoord[4],
                                                          faceCoord[3],
                                                          eg.elementGlobal,
                                                          faceCoord[1]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          edgeCoord[1],
                                                          faceCoord[0],
                                                          edgeCoord[2],
                                                          edgeCoord[5],
                                                          faceCoord[2],
                                                          eg.elementGlobal,
                                                          faceCoord[3]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          edgeCoord[0],
                                                          faceCoord[0],
                                                          edgeCoord[3],
                                                          edgeCoord[6],
                                                          faceCoord[1],
                                                          eg.elementGlobal,
                                                          faceCoord[4]);
            eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                          edgeCoord[3],
                                                          faceCoord[0],
                                                          edgeCoord[1],
                                                          edgeCoord[7],
                                                          faceCoord[4],
                                                          eg.elementGlobal,
                                                          faceCoord[2]);
            eg.subContVol[4].volume = eg.elementVolume -
                eg.subContVol[0].volume -
                eg.subContVol[1].volume -
                eg.subContVol[2].volume -
                eg.subContVol[3].volume;
            break;
        case 6: // 3D, prism
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          edgeCoord[3],
                                                          faceCoord[3],
                                                          edgeCoord[4],
                                                          edgeCoord[0],
                                                          faceCoord[0],
                                                          eg.elementGlobal,
                                                          faceCoord[1]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          edgeCoord[5],
                                                          faceCoord[3],
                                                          edgeCoord[3],
                                                          edgeCoord[1],
                                                          faceCoord[2],
                                                          eg.elementGlobal,
                                                          faceCoord[0]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          edgeCoord[4],
                                                          faceCoord[3],
                                                          edgeCoord[5],
                                                          edgeCoord[2],
                                                          faceCoord[1],
                                                          eg.elementGlobal,
                                                          faceCoord[2]);
            eg.subContVol[3].volume = eg.hexahedronVolume(edgeCoord[0],
                                                          faceCoord[0],
                                                          eg.elementGlobal,
                                                          faceCoord[1],
                                                          eg.subContVol[3].global,
                                                          edgeCoord[6],
                                                          faceCoord[4],
                                                          edgeCoord[7]);
            eg.subContVol[4].volume = eg.hexahedronVolume(edgeCoord[1],
                                                          faceCoord[2],
                                                          eg.elementGlobal,
                                                          faceCoord[0],
                                                          eg.subContVol[4].global,
                                                          edgeCoord[8],
                                                          faceCoord[4],
                                                          edgeCoord[6]);
            eg.subContVol[5].volume = eg.hexahedronVolume(edgeCoord[2],
                                                          faceCoord[1],
                                                          eg.elementGlobal,
                                                          faceCoord[2],
                                                          eg.subContVol[5].global,
                                                          edgeCoord[7],
                                                          faceCoord[4],
                                                          edgeCoord[8]);
            break;
        case 8: // 3D, hexahedron
            eg.subContVol[0].volume = eg.hexahedronVolume(eg.subContVol[0].global,
                                                          edgeCoord[6],
                                                          faceCoord[4],
                                                          edgeCoord[4],
                                                          edgeCoord[0],
                                                          faceCoord[2],
                                                          eg.elementGlobal,
                                                          faceCoord[0]);
            eg.subContVol[1].volume = eg.hexahedronVolume(eg.subContVol[1].global,
                                                          edgeCoord[5],
                                                          faceCoord[4],
                                                          edgeCoord[6],
                                                          edgeCoord[1],
                                                          faceCoord[1],
                                                          eg.elementGlobal,
                                                          faceCoord[2]);
            eg.subContVol[2].volume = eg.hexahedronVolume(eg.subContVol[2].global,
                                                          edgeCoord[4],
                                                          faceCoord[4],
                                                          edgeCoord[7],
                                                          edgeCoord[2],
                                                          faceCoord[0],
                                                          eg.elementGlobal,
                                                          faceCoord[3]);
            eg.subContVol[3].volume = eg.hexahedronVolume(eg.subContVol[3].global,
                                                          edgeCoord[7],
                                                          faceCoord[4],
                                                          edgeCoord[5],
                                                          edgeCoord[3],
                                                          faceCoord[3],
                                                          eg.elementGlobal,
                                                          faceCoord[1]);
            eg.subContVol[4].volume = eg.hexahedronVolume(edgeCoord[0],
                                                          faceCoord[2],
                                                          eg.elementGlobal,
                                                          faceCoord[0],
                                                          eg.subContVol[4].global,
                                                          edgeCoord[10],
                                                          faceCoord[5],
                                                          edgeCoord[8]);
            eg.subContVol[5].volume = eg.hexahedronVolume(edgeCoord[1],
                                                          faceCoord[1],
                                                          eg.elementGlobal,
                                                          faceCoord[2],
                                                          eg.subContVol[5].global,
                                                          edgeCoord[9],
                                                          faceCoord[5],
                                                          edgeCoord[10]);
            eg.subContVol[6].volume = eg.hexahedronVolume(edgeCoord[2],
                                                          faceCoord[0],
                                                          eg.elementGlobal,
                                                          faceCoord[3],
                                                          eg.subContVol[6].global,
                                                          edgeCoord[8],
                                                          faceCoord[5],
                                                          edgeCoord[11]);
            eg.subContVol[7].volume = eg.hexahedronVolume(edgeCoord[3],
                                                          faceCoord[3],
                                                          eg.elementGlobal,
                                                          faceCoord[1],
                                                          eg.subContVol[7].global,
                                                          edgeCoord[11],
                                                          faceCoord[5],
                                                          edgeCoord[9]);
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
 * \ingroup BoxModel
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
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};

    typedef BoxFVElementGeometry<TypeTag>   ThisType;

    /** \todo Please doc me! */
    friend class _BoxFVElemGeomHelper<ThisType, dim>;

    typedef _BoxFVElemGeomHelper<ThisType, dim> BoxFVElemGeomHelper;

    enum{maxNC = (dim < 3 ? 4 : 8)};
    enum{maxNE = (dim < 3 ? 4 : 12)};
    enum{maxNF = (dim < 3 ? 1 : 6)};
    enum{maxCOS = (dim < 3 ? 2 : 4)};
    enum{maxBF = (dim < 3 ? 8 : 24)};
    enum{maxNFAP = (dim < 3 ? 4 : 8)};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Geometry::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;

    Scalar quadrilateralArea(const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2,
                             const GlobalPosition& p3)
    {
        return 0.5*fabs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
    }

    void crossProduct(GlobalPosition& c, const GlobalPosition& a, const GlobalPosition& b)
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    Scalar pyramidVolume (const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2,
                          const GlobalPosition& p3, const GlobalPosition& p4)
    {
        GlobalPosition a(p2); a -= p0;
        GlobalPosition b(p3); b -= p1;

        GlobalPosition n;
        crossProduct(n, a, b);

        a = p4; a -= p0;

        return 1.0/6.0*(n*a);
    }

    Scalar prismVolume (const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2,
                        const GlobalPosition& p3, const GlobalPosition& p4, const GlobalPosition& p5)
    {
        GlobalPosition a(p4);
        for (int k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];
        GlobalPosition b(p1);
        for (int k = 0; k < dimWorld; ++k)
            b[k] -= p3[k];
        GlobalPosition m;
        crossProduct(m, a, b);

        a = p1;
        for (int k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];
        b = p2;
        for (int k = 0; k < dimWorld; ++k)
            b[k] -= p0[k];
        GlobalPosition n;
        crossProduct(n, a, b);
        n += m;

        a = p5;
        for (int k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];

        return fabs(1.0/6.0*(n*a));
    }

    Scalar hexahedronVolume (const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2, const GlobalPosition& p3,
                             const GlobalPosition& p4, const GlobalPosition& p5, const GlobalPosition& p6, const GlobalPosition& p7)
    {
        return
            prismVolume(p0,p1,p2,p4,p5,p6)
            + prismVolume(p0,p2,p3,p4,p6,p7);
    }

    void normalOfQuadrilateral3D(GlobalPosition &normal, const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2,
                                 const GlobalPosition& p3)
    {
        GlobalPosition a(p2);
        for (int k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];
        GlobalPosition b(p3);
        for (int k = 0; k < dimWorld; ++k)
            b[k] -= p1[k];

        crossProduct(normal, a, b);
        normal *= 0.5;
    }

    Scalar quadrilateralArea3D(const GlobalPosition& p0, const GlobalPosition& p1, const GlobalPosition& p2, const GlobalPosition& p3)
    {
        GlobalPosition normal;
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
            {1, 0, 2, 0, 3, 2, 4, 1, 4},
            {0, 2, 1, 3, 1, 3, 0, 4, 2}
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

    void getEdgeIndices(int numVertices, int fIdx, int vIdx, int& leftEdge, int& rightEdge)
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
            leftEdge = faceAndVertexToLeftEdgeTet[fIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgeTet[fIdx][vIdx];
            break;
        case 5:
            leftEdge = faceAndVertexToLeftEdgePyramid[fIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgePyramid[fIdx][vIdx];
            break;
        case 6:
            leftEdge = faceAndVertexToLeftEdgePrism[fIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgePrism[fIdx][vIdx];
            break;
        case 8:
            leftEdge = faceAndVertexToLeftEdgeHex[fIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgeHex[fIdx][vIdx];
            break;
        default:
            DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry :: getFaceIndices for numVertices = " << numVertices);
            break;
        }
    }

public:
    int boundaryFaceIndex(const int fIdx, const int vertInFace) const
    {
        return (fIdx*maxCOS + vertInFace);
    }

    struct SubControlVolume //! FV intersected with element
    {
        LocalPosition local; //!< local position
        GlobalPosition global; //!< global position
        LocalPosition localCenter; //!< local position of scv center
        Scalar volume; //!< volume of scv
        Dune::FieldVector<GlobalPosition, maxNC> grad; //! derivative of shape function associated with the sub control volume
        Dune::FieldVector<GlobalPosition, maxNC> gradCenter; //! derivative of shape function at the center of the sub control volume
        Dune::FieldVector<Scalar, maxNC> shapeValue; //! value of shape function associated with the sub control volume
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
        LocalPosition ipLocal; //!< integration point in local coords
        GlobalPosition ipGlobal; //!< integration point in global coords
        GlobalPosition normal; //!< normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        Scalar area; //!< area of face
        Dune::FieldVector<GlobalPosition, maxNC> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNC> shapeValue; //!< value of shape functions at ip
        Dune::FieldVector<int, maxNFAP> fapIndices; //!< indices w.r.t.neighbors of the flux approximation points
        unsigned numFap; //!< number of flux approximation points
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[maxNC]; //!< data of the sub control volumes
    SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    int numScv; //!< number of subcontrol volumes
    int numScvf; //!< number of inner-domain subcontrolvolume faces 
    int numNeighbors; //!< needed for compatibility with cc models
    std::vector<ElementPointer> neighbors; //!< needed for compatibility with cc models

    const LocalFiniteElementCache feCache_;

    void updateInner(const Element& element) //!< needed for compatibility with cc models
    {}

    void update(const GridView& gridView, const Element& element)
    {
        const Geometry& geometry = element.geometry();
        Dune::GeometryType geomType = geometry.type();

        const ReferenceElement &referenceElement = ReferenceElements::general(geomType);

        const LocalFiniteElement
            &localFiniteElement = feCache_.get(geomType);

        elementVolume = geometry.volume();
        elementLocal = referenceElement.position(0,0);
        elementGlobal = geometry.global(elementLocal);

        numScv = referenceElement.size(dim);
        numScvf = referenceElement.size(dim-1);
        numNeighbors = 0;

        // subcontrol volumes:
        for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
            subContVol[scvIdx].local = referenceElement.position(scvIdx, dim);
            subContVol[scvIdx].global = geometry.global(subContVol[scvIdx].local);
            subContVol[scvIdx].inner = true;
        }

        // edges:
        GlobalPosition *edgeCoordinates = new GlobalPosition[numScvf];
        for (int edge = 0; edge < numScvf; edge++) {
            edgeCoordinates[edge] = geometry.global(referenceElement.position(edge, dim-1));
        }

        // faces:
        int elementFaces = (dim<3)?0:referenceElement.size(1);
        GlobalPosition *faceCoordinates = new GlobalPosition[elementFaces];
        for (int fIdx = 0; fIdx < elementFaces; fIdx++) {
            faceCoordinates[fIdx] = geometry.global(referenceElement.position(fIdx, 1));
        }

        // fill sub control volume data use specialization for this
        // \todo maybe it would be a good idea to migrate everything
        // which is dependend of the grid's dimension to
        // _BoxFVElemGeomHelper in order to benefit from more aggressive
        // compiler optimizations...
        BoxFVElemGeomHelper::fillSubContVolData(*this, numScv, edgeCoordinates, faceCoordinates);

        // fill sub control volume face data:
        for (int k = 0; k < numScvf; k++) { // begin loop over edges / sub control volume faces
            int i = referenceElement.subEntity(k, dim-1, 0, dim);
            int j = referenceElement.subEntity(k, dim-1, 1, dim);
            if (numScvf == 4 && (i == 2 || j == 2))
                std::swap(i, j);
            subContVolFace[k].i = i;
            subContVolFace[k].j = j;

            // calculate the local integration point and
            // the face normal. note that since dim is a
            // constant which is known at compile time
            // the compiler can optimize away all if
            // cases which don't apply.
            LocalPosition ipLocal;
            GlobalPosition diffVec;
            if (dim==1) {
                subContVolFace[k].ipLocal = 0.5;
                subContVolFace[k].normal = 1.0;
                ipLocal = subContVolFace[k].ipLocal;
            }
            else if (dim==2) {
                ipLocal = referenceElement.position(k, dim-1) + elementLocal;
                ipLocal *= 0.5;
                subContVolFace[k].ipLocal = ipLocal;
                diffVec = elementGlobal - edgeCoordinates[k];
                subContVolFace[k].normal[0] = diffVec[1];
                subContVolFace[k].normal[1] = -diffVec[0];

                diffVec = subContVol[j].global;
                for (int m = 0; m < dimWorld; ++m)
                    diffVec[m] -= subContVol[i].global[m];
                // make sure the normal points to the right direction
                if (subContVolFace[k].normal * diffVec < 0)
                    subContVolFace[k].normal *= -1.0;
            }
            else if (dim==3) {
                int leftFace;
                int rightFace;
                getFaceIndices(numScv, k, leftFace, rightFace);
                ipLocal = referenceElement.position(k, dim-1) + elementLocal
                    + referenceElement.position(leftFace, 1)
                    + referenceElement.position(rightFace, 1);
                ipLocal *= 0.25;
                subContVolFace[k].ipLocal = ipLocal;
                normalOfQuadrilateral3D(subContVolFace[k].normal,
                                        edgeCoordinates[k], faceCoordinates[rightFace],
                                        elementGlobal, faceCoordinates[leftFace]);
            }

            subContVolFace[k].area = subContVolFace[k].normal.two_norm();

            bool useTwoPointFlux
                = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, UseTwoPointFlux);

            if (useTwoPointFlux)
            {
                GlobalPosition distVec = subContVol[i].global;
                distVec -= subContVol[j].global;
                distVec /= distVec.two_norm2();

                // gradients using a two-point flux approximation
                subContVolFace[k].numFap = 2;
                for (unsigned int idx = 0; idx < subContVolFace[k].numFap; idx++)
                {
                    subContVolFace[k].grad[idx] = distVec;
                    subContVolFace[k].shapeValue[idx] = 0.5;
                }
                subContVolFace[k].grad[1] *= -1.0;

                subContVolFace[k].fapIndices[0] = subContVolFace[k].i;
                subContVolFace[k].fapIndices[1] = subContVolFace[k].j;
            }
            else
            {
                // get the global integration point and the Jacobian inverse
                subContVolFace[k].ipGlobal = geometry.global(ipLocal);
                JacobianInverseTransposed jacInvT = 
                    geometry.jacobianInverseTransposed(ipLocal);

                // calculate the shape function gradients
                //typedef Dune::FieldVector< Dune::FieldVector< CoordScalar, dim >, 1 > ShapeJacobian;
                typedef Dune::FieldVector< Scalar, 1 > ShapeValue;
                std::vector<ShapeJacobian> localJac;
                std::vector<ShapeValue>    shapeVal;
                localFiniteElement.localBasis().evaluateJacobian(subContVolFace[k].ipLocal, localJac);
                localFiniteElement.localBasis().evaluateFunction(subContVolFace[k].ipLocal, shapeVal);
                subContVolFace[k].numFap = numScv;
                for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
                    jacInvT.mv(localJac[scvIdx][0], subContVolFace[k].grad[scvIdx]);
                    subContVolFace[k].shapeValue[scvIdx] = Scalar(shapeVal[scvIdx]);
                    subContVolFace[k].fapIndices[scvIdx] = scvIdx;
                }
            }
        } // end loop over edges / sub control volume faces

        // fill boundary face data:
        IntersectionIterator isEndIt = gridView.iend(element);
        for (IntersectionIterator isIt = gridView.ibegin(element); isIt != isEndIt; ++isIt)
            if (isIt->boundary())
            {
                int fIdx = isIt->indexInInside();
                int numVerticesOfFace = referenceElement.size(fIdx, 1, dim);
                for (int vertInFace = 0; vertInFace < numVerticesOfFace; vertInFace++)
                {
                    int vertInElement = referenceElement.subEntity(fIdx, 1, vertInFace, dim);
                    int bfIdx = boundaryFaceIndex(fIdx, vertInFace);
                    subContVol[vertInElement].inner = false;
                    switch ((short) dim) {
                    case 1:
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim);
                        boundaryFace[bfIdx].area = 1.0;
                        break;
                    case 2:
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim)
                            + referenceElement.position(fIdx, 1);
                        boundaryFace[bfIdx].ipLocal *= 0.5;
                        boundaryFace[bfIdx].area = 0.5*isIt->geometry().volume();
                        break;
                    case 3:
                        int leftEdge;
                        int rightEdge;
                        getEdgeIndices(numScv, fIdx, vertInElement, leftEdge, rightEdge);
                        boundaryFace[bfIdx].ipLocal = referenceElement.position(vertInElement, dim)
                            + referenceElement.position(fIdx, 1)
                            + referenceElement.position(leftEdge, dim-1)
                            + referenceElement.position(rightEdge, dim-1);
                        boundaryFace[bfIdx].ipLocal *= 0.25;
                        boundaryFace[bfIdx].area = quadrilateralArea3D(subContVol[vertInElement].global,
                                                                       edgeCoordinates[rightEdge], faceCoordinates[fIdx],
                                                                       edgeCoordinates[leftEdge]);
                        break;
                    default:
                        DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry for dim = " << dim);
                    }
                    boundaryFace[bfIdx].ipGlobal = geometry.global(boundaryFace[bfIdx].ipLocal);
                    boundaryFace[bfIdx].i = vertInElement;
                    boundaryFace[bfIdx].j = vertInElement;

                    // ASSUME constant normal
                    Dune::FieldVector<CoordScalar, dim-1> localDimM1(0);
                    boundaryFace[bfIdx].normal = isIt->unitOuterNormal(localDimM1);
                    boundaryFace[bfIdx].normal *= boundaryFace[bfIdx].area;

                    typedef Dune::FieldVector< Scalar, 1 > ShapeValue;
                    std::vector<ShapeJacobian> localJac;
                    std::vector<ShapeValue>    shapeVal;
                    localFiniteElement.localBasis().evaluateJacobian(boundaryFace[bfIdx].ipLocal, localJac);
                    localFiniteElement.localBasis().evaluateFunction(boundaryFace[bfIdx].ipLocal, shapeVal);

                    JacobianInverseTransposed jacInvT = 
                        geometry.jacobianInverseTransposed(boundaryFace[bfIdx].ipLocal);
                    boundaryFace[bfIdx].numFap = numScv;
                    for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
                    {
                        jacInvT.mv(localJac[scvIdx][0], boundaryFace[bfIdx].grad[scvIdx]);
                        boundaryFace[bfIdx].shapeValue[scvIdx] = Scalar(shapeVal[scvIdx]);
                        boundaryFace[bfIdx].fapIndices[scvIdx] = scvIdx;
                    }
                }
            }

        bool evalGradientsAtSCVCenter = GET_PROP_VALUE(TypeTag, EvalGradientsAtSCVCenter);
        if(evalGradientsAtSCVCenter)
        {
            // calculate gradients at the center of the scv
            for (int scvIdx = 0; scvIdx < numScv; scvIdx++){
                if (dim == 2)
                {
                    switch (scvIdx)
                    {
                    case 0:
                        if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 0.25;
                            subContVol[scvIdx].localCenter[1] = 0.25;
                        }
                        else {
                            subContVol[scvIdx].localCenter[0] = 1.0/6.0;
                            subContVol[scvIdx].localCenter[1] = 1.0/6.0;
                        }
                        break;
                    case 1:
                        if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 0.75;
                            subContVol[scvIdx].localCenter[1] = 0.25;
                        }
                        else {
                            subContVol[scvIdx].localCenter[0] = 4.0/6.0;
                            subContVol[scvIdx].localCenter[1] = 1.0/6.0;
                        }
                        break;
                    case 2:
                        if (numScv == 4) {
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

                else if (dim == 3)
                {
                    switch (scvIdx)
                    {
                    case 0:
                        if (numScv == 8) {
                            subContVol[scvIdx].localCenter[0] = 0.25;
                            subContVol[scvIdx].localCenter[1] = 0.25;
                            subContVol[scvIdx].localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[1] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 1:
                        if (numScv == 8) {
                            subContVol[scvIdx].localCenter[0] = 0.75;
                            subContVol[scvIdx].localCenter[1] = 0.25;
                            subContVol[scvIdx].localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 7.0/16.0;
                            subContVol[scvIdx].localCenter[1] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 2:
                        if (numScv == 8) {
                            subContVol[scvIdx].localCenter[0] = 0.25;
                            subContVol[scvIdx].localCenter[1] = 0.75;
                            subContVol[scvIdx].localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[1] = 7.0/16.0;
                            subContVol[scvIdx].localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 3:
                        if (numScv == 8) {
                            subContVol[scvIdx].localCenter[0] = 0.75;
                            subContVol[scvIdx].localCenter[1] = 0.75;
                            subContVol[scvIdx].localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subContVol[scvIdx].localCenter[0] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[1] = 3.0/16.0;
                            subContVol[scvIdx].localCenter[2] = 7.0/16.0;
                        }
                        break;
                    case 4:
                        subContVol[scvIdx].localCenter[0] = 0.25;
                        subContVol[scvIdx].localCenter[1] = 0.25;
                        subContVol[scvIdx].localCenter[2] = 0.75;
                        break;
                    case 5:
                        subContVol[scvIdx].localCenter[0] = 0.75;
                        subContVol[scvIdx].localCenter[1] = 0.25;
                        subContVol[scvIdx].localCenter[2] = 0.75;
                        break;
                    case 6:
                        subContVol[scvIdx].localCenter[0] = 0.25;
                        subContVol[scvIdx].localCenter[1] = 0.75;
                        subContVol[scvIdx].localCenter[2] = 0.75;
                        break;
                    case 7:
                        subContVol[scvIdx].localCenter[0] = 0.75;
                        subContVol[scvIdx].localCenter[1] = 0.75;
                        subContVol[scvIdx].localCenter[2] = 0.75;
                        break;
                    }
                }
                std::vector<ShapeJacobian> localJac;
                localFiniteElement.localBasis().evaluateJacobian(subContVol[scvIdx].localCenter, localJac);

                JacobianInverseTransposed jacInvT =
                    geometry.jacobianInverseTransposed(subContVol[scvIdx].localCenter);
                for (int vIdx = 0; vIdx < numScv; vIdx++)
                    jacInvT.mv(localJac[vIdx][0], subContVol[scvIdx].gradCenter[vIdx]);
            }
        }

        delete[] edgeCoordinates;
        delete[] faceCoordinates;
    }
};
}

#endif

