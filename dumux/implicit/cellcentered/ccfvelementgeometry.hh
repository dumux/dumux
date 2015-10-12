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
 *        the cell-centered fv scheme.
 */
#ifndef DUMUX_CC_FV_ELEMENTGEOMETRY_HH
#define DUMUX_CC_FV_ELEMENTGEOMETRY_HH

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
}

/*!
 * \ingroup CCModel
 * \brief Represents the finite volume geometry of a single element in
 *        the cell-centered fv scheme.
 */
template<class TypeTag>
class CCFVElementGeometry
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};

    enum{maxNFAP = 2}; //! maximum number of flux approximation points (two-point flux)
    enum{maxNE = 2*dim*(1 << (dim - 1))}; //! maximum number of neighbors (works for one hanging node per face)
    enum{maxBF = 2*dim}; //! maximum number of boundary faces
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    struct SubControlVolume //! FV intersected with element
    {
        LocalPosition local; //!< local position
        GlobalPosition global; //!< global position
        Scalar volume; //!< volume of scv
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
        LocalPosition ipLocal; //!< integration point in local coords
        GlobalPosition ipGlobal; //!< integration point in global coords
        GlobalPosition normal; //!< normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        Scalar area; //!< area of face
        Dune::FieldVector<GlobalPosition, maxNFAP> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNFAP> shapeValue; //!< value of shape functions at ip
        Dune::FieldVector<int, maxNFAP> fapIndices; //!< indices w.r.t.neighbors of the flux approximation points
        unsigned numFap; //!< number of flux approximation points
        unsigned fIdx; //!< the index (w.r.t. the element) of the face (codim 1 entity) that the scvf is part of
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[1]; //!< data of the sub control volumes
    SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    int numScv; //!< number of subcontrol volumes
    int numScvf; //!< number of inner-domain subcontrolvolume faces
    int numNeighbors; //!< number of neighboring elements including the element itself
    std::vector<Element> neighbors; //!< stores the neighboring elements

    void updateInner(const Element& element)
    {
        const Geometry geometry = element.geometry();

        elementVolume = geometry.volume();
        elementGlobal = geometry.center();
        elementLocal = geometry.local(elementGlobal);

        numScv = 1;
        numScvf = 0;

        subContVol[0].local = elementLocal;
        subContVol[0].global = elementGlobal;
        subContVol[0].inner = true;
        subContVol[0].volume = elementVolume;

        // initialize neighbors list with self:
        numNeighbors = 1;
        neighbors.clear();
        neighbors.reserve(maxNE);
        neighbors.push_back(element);
    }

    void update(const GridView& gridView, const Element& element)
    {
        updateInner(element);

        const Geometry geometry = element.geometry();

        bool onBoundary = false;

        // fill neighbor information and control volume face data:
        IntersectionIterator isEndIt = gridView.iend(element);
        for (IntersectionIterator isIt = gridView.ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto isGeometry = isIt->geometry();

            // neighbor information and inner cvf data:
            if (isIt->neighbor())
            {
                numNeighbors++;
                neighbors.push_back(isIt->outside());

                int scvfIdx = numNeighbors - 2;
                SubControlVolumeFace& scvFace = subContVolFace[scvfIdx];

                scvFace.i = 0;
                scvFace.j = scvfIdx + 1;

                scvFace.ipGlobal = isGeometry.center();
                scvFace.ipLocal =  geometry.local(scvFace.ipGlobal);
                scvFace.normal = isIt->centerUnitOuterNormal();
                Scalar volume = isGeometry.volume();
                scvFace.normal *= volume;
                scvFace.area = volume;

                GlobalPosition distVec = elementGlobal
                                       - neighbors[scvfIdx+1].geometry().center();
                distVec /= distVec.two_norm2();

                // gradients using a two-point flux approximation
                scvFace.numFap = 2;
                for (unsigned int fapIdx = 0; fapIdx < scvFace.numFap; fapIdx++)
                {
                    scvFace.grad[fapIdx] = distVec;
                    scvFace.shapeValue[fapIdx] = 0.5;
                }
                scvFace.grad[1] *= -1.0;

                scvFace.fapIndices[0] = scvFace.i;
                scvFace.fapIndices[1] = scvFace.j;

                scvFace.fIdx = isIt->indexInInside();
            }

            // boundary cvf data
            if (isIt->boundary())
            {
                onBoundary = true;
                int bfIdx = isIt->indexInInside();
                SubControlVolumeFace& bFace = boundaryFace[bfIdx];

                bFace.ipGlobal = isGeometry.center();
                bFace.ipLocal =  geometry.local(bFace.ipGlobal);
                bFace.normal = isIt->centerUnitOuterNormal();
                Scalar volume = isGeometry.volume();
                bFace.normal *= volume;
                bFace.area = volume;
                bFace.i = 0;
                bFace.j = 0;

                GlobalPosition distVec = elementGlobal - bFace.ipGlobal;
                distVec /= distVec.two_norm2();

                // gradients using a two-point flux approximation
                bFace.numFap = 2;
                for (unsigned int fapIdx = 0; fapIdx < bFace.numFap; fapIdx++)
                {
                    bFace.grad[fapIdx] = distVec;
                    bFace.shapeValue[fapIdx] = 0.5;
                }
                bFace.grad[1] *= -1.0;

                bFace.fapIndices[0] = bFace.i;
                bFace.fapIndices[1] = bFace.j;
            }
        }

        // set the number of inner-domain subcontrolvolume faces
        numScvf = numNeighbors - 1;

        // treat elements on the boundary
        if (onBoundary)
        {
            for (int bfIdx = 0; bfIdx < element.subEntities(1); bfIdx++)
            {
                SubControlVolumeFace& bFace = boundaryFace[bfIdx];
                bFace.j = numNeighbors + bfIdx;
                bFace.fapIndices[1] = bFace.j;
                neighbors.push_back(element);
            }
        }
    }

    /*!
     * \brief For compatibilty with the box element geometry
     */
    int boundaryFaceIndex(const int fIdx, const int vIdxInFace) const
    {
        return fIdx;
    }
};

}

#endif

