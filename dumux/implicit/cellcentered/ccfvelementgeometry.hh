// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Bernd Flemisch, Andreas Lauser               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Represents the finite volume geometry of a single element in
 *        the cell-centered fv scheme.
 */
#ifndef DUMUX_CC_FV_ELEMENTGEOMETRY_HH
#define DUMUX_CC_FV_ELEMENTGEOMETRY_HH

#include <dune/common/version.hh>
#include <dune/grid/common/intersectioniterator.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_COMMON, 2, 2, 0)
// dune 2.2
#include <dune/geometry/referenceelements.hh>
#else
// dune 2.1
#include <dune/geometry/referenceelements.hh>
#endif

#include <dune/localfunctions/lagrange/pqkfactory.hh>

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

    enum{maxNFAP = 2};
    enum{maxNE = (dim < 3 ? 4 : 12)};
    enum{maxNF = (dim < 3 ? 1 : 6)};
    enum{maxCOS = (dim < 3 ? 2 : 4)};
    enum{maxBF = (dim < 3 ? 8 : 24)};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    struct SubControlVolume //! FV intersected with element
    {
        LocalPosition local; //!< local vert position
        GlobalPosition global; //!< global vert position
        Scalar volume; //!< volume of scv
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
        LocalPosition ipLocal; //!< integration point in local coords
        GlobalPosition ipGlobal; //!< integration point in global coords
        Vector normal; //!< normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        Scalar area; //!< area of face
        Dune::FieldVector<Vector, maxNFAP> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNFAP> shapeValue; //!< value of shape functions at ip
        Dune::FieldVector<int, maxNFAP> fapIndices; //!< indices w.r.t.neighbors of the flux approximation points
        unsigned faceIdx;
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[1]; //!< data of the sub control volumes
    SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    int numVertices; //!< number of verts
    int numEdges; //!< number of edges
    int numFaces; //!< number of faces (0 in < 3D)
    int numScv; //!< number of subcontrol volumes
    int numSCV; //!< \deprecated number of subcontrol volumes
    int numNeighbors; //!< number of neighboring elements including the element itself
    int numFap; //!< number of flux approximation points
    int numFAP; //!< \deprecated number of flux approximation points
    std::vector<ElementPointer> neighbors; //!< stores pointers for the neighboring elements
    
    void updateInner(const Element& element)
    {
        const Geometry& geometry = element.geometry();

        elementVolume = geometry.volume();
        elementGlobal = geometry.center();
        elementLocal = geometry.local(elementGlobal);

        numVertices = element.template count<dim>();
        numEdges = element.template count<dim-1>();
        numFaces = (dim<3)? 0 : element.template count<1>();
        numScv = 1;
        numSCV = numScv;
        numFap = 2;
        numFAP = numFap;
        
        subContVol[0].local = elementLocal;
        subContVol[0].global = elementGlobal;
        subContVol[0].inner = true;
        subContVol[0].volume = elementVolume;

        // initialize neighbors list with self: 
        numNeighbors = 1;
        neighbors.clear();
        neighbors.reserve(maxNE);
        ElementPointer elementPointer(element);
        neighbors.push_back(elementPointer);
    }
    
    void update(const GridView& gridView, const Element& element)
    {
        updateInner(element);
        
        const Geometry& geometry = element.geometry();

        // fill neighbor information and control volume face data:
        IntersectionIterator endit = gridView.iend(element);
        for (IntersectionIterator it = gridView.ibegin(element); it != endit; ++it)
        {
            // neighbor information and inner cvf data:
            if (it->neighbor())
            {
                numNeighbors++;
                ElementPointer elementPointer(it->outside());
                neighbors.push_back(elementPointer);
                
                int k = numNeighbors - 2;
                
                subContVolFace[k].i = 0;
                subContVolFace[k].j = k+1;
                
                subContVolFace[k].ipGlobal = it->geometry().center();
                subContVolFace[k].ipLocal =  geometry.local(subContVolFace[k].ipGlobal);
                subContVolFace[k].normal = it->centerUnitOuterNormal();
                subContVolFace[k].normal *= it->geometry().volume();
                subContVolFace[k].area = it->geometry().volume();

                GlobalPosition distVec = elementGlobal;
                distVec -= neighbors[k+1]->geometry().center();
                distVec /= distVec.two_norm2();
                
                // gradients using a two-point flux approximation
                for (int idx = 0; idx < 2; idx++)
                {
                    subContVolFace[k].grad[idx] = distVec;
                    subContVolFace[k].shapeValue[idx] = 0.5;
                }
                subContVolFace[k].grad[1] *= -1.0;
                
                subContVolFace[k].fapIndices[0] = subContVolFace[k].i;
                subContVolFace[k].fapIndices[1] = subContVolFace[k].j;
                
                subContVolFace[k].faceIdx = it->indexInInside();
            }

            // boundary cvf data
            if (it->boundary())
            {
                int bfIdx = it->indexInInside();
                boundaryFace[bfIdx].ipGlobal = it->geometry().center();
                boundaryFace[bfIdx].ipLocal =  geometry.local(boundaryFace[bfIdx].ipGlobal);
                boundaryFace[bfIdx].normal = it->centerUnitOuterNormal();
                boundaryFace[bfIdx].normal *= it->geometry().volume();
                boundaryFace[bfIdx].area = it->geometry().volume();
                boundaryFace[bfIdx].i = 0;
                boundaryFace[bfIdx].j = 0;

                GlobalPosition distVec = elementGlobal;
                distVec -= boundaryFace[bfIdx].ipGlobal;
                distVec /= 2.0*distVec.two_norm2();
                
                // gradients using a two-point flux approximation
                for (int idx = 0; idx < 2; idx++)
                {
                    boundaryFace[bfIdx].grad[idx] = distVec;
                    boundaryFace[bfIdx].shapeValue[idx] = 0.5;
                }
                boundaryFace[bfIdx].grad[0] *= -1.0;
                
                boundaryFace[bfIdx].fapIndices[0] = boundaryFace[bfIdx].i;
                boundaryFace[bfIdx].fapIndices[1] = boundaryFace[bfIdx].j;
            }
        }
        
        for (int nIdx = 0; nIdx < numNeighbors-1; nIdx++)
        {
            switch (subContVolFace[nIdx].faceIdx)
            {
                case 0:
                    boundaryFace[1].j = nIdx+1;
                    boundaryFace[1].fapIndices[1] = boundaryFace[1].j;
                    break;
                case 1:
                    boundaryFace[0].j = nIdx+1;
                    boundaryFace[0].fapIndices[1] = boundaryFace[0].j;
                    break;
                case 2:
                    boundaryFace[3].j = nIdx+1;
                    boundaryFace[3].fapIndices[1] = boundaryFace[3].j;
                    break;
                case 3:
                    boundaryFace[2].j = nIdx+1;
                    boundaryFace[2].fapIndices[1] = boundaryFace[2].j;
                    break;
                case 4:
                    boundaryFace[5].j = nIdx+1;
                    boundaryFace[5].fapIndices[1] = boundaryFace[5].j;
                    break;
                case 5:
                    boundaryFace[4].j = nIdx+1;
                    boundaryFace[4].fapIndices[1] = boundaryFace[4].j;
                    break;
            }
        }
    }
};

}

#endif

