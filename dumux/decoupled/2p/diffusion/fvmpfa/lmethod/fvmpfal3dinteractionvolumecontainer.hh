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
#ifndef DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_HH
#define DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_HH

// dumux environment
#include <dumux/decoupled/common/pressureproperties.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include <dumux/decoupled/common/fv/mpfa/mpfalinteractionvolume3d.hh>

/**
 * @file
 * @brief  Interactionvolume container for 3-d MPFA L-method
 */

namespace Dumux
{

bool sort_compare(const std::vector<int>& entryI, const std::vector<int>& entryJ)
{
    return (entryI[1] < entryJ[1]);
}

//! \ingroup diffusion
/*! Interactionvolume container for 3-d MPFA L-method
 */
template<class TypeTag>
class FvMpfaL3dInteractionVolumeContainer
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolumeContainer) Implementation;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;
#else
    typedef typename Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dim>::Entity Vertex;
    typedef typename Element::Geometry ElementGeometry;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    typedef typename Intersection::Geometry IntersectionGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

    enum
        {
            pressureEqIdx = Indices::pressureEqIdx,
        };

    enum
        {
            innerEdgeFace = 2, innerSideFace = 1
        };
    enum
        {
            realFaceArea = 0, fluxFaceArea = 1
        };

public:
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolume) InteractionVolume;

private:
    typedef std::vector<InteractionVolume> GlobalInteractionVolumeVector;
    typedef std::vector<Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2 * dim> > FaceAreaVector;
protected:
    void storeSubVolumeElements(const Element& element, std::vector < std::vector<int> >& elemVertMap);
    void storeIntersectionInfo(const Element& element, std::vector < std::vector<int> >& elemVertMap);
    void storeInnerInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex, bool sameLevel = true);
    void storeBoundaryInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex);
private:
    void storeInteractionVolumeInfo();
public:

    void update()
    {
        interactionVolumes_.clear();
        realFluxFaceArea_.clear();

        realFluxFaceArea_.resize(problem_.gridView().size(dim), Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2 * dim>(Dune::FieldVector<Scalar, 2>(0.0)));
        interactionVolumes_.resize(problem_.gridView().size(dim));

        asImp_().storeInteractionVolumeInfo();
    }

    void initialize(bool solveTwice = true)
    {
        interactionVolumes_.clear();
        realFluxFaceArea_.clear();

        realFluxFaceArea_.resize(problem_.gridView().size(dim), Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2 * dim>(Dune::FieldVector<Scalar, 2>(0.0)));
        interactionVolumes_.resize(problem_.gridView().size(dim));

        asImp_().storeInteractionVolumeInfo();

        return;
    }

    InteractionVolume& interactionVolume(int vertexIdx)
    {
        return interactionVolumes_[vertexIdx];
    }

    InteractionVolume& interactionVolume(int vertexIdx) const
    {
        return interactionVolumes_[vertexIdx];
    }

    GlobalInteractionVolumeVector& interactionVolumesGlobal()
    {
        return interactionVolumes_;
    }

    GlobalInteractionVolumeVector& interactionVolumesGlobal() const
    {
        return interactionVolumes_;
    }

    Scalar faceAreaFactor(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = getRealFaceArea(interactionVolume, elemGlobalIdx, elemLocalIdx, localFaceIdx);
        factor /= getRealFluxFaceArea(interactionVolume, elemGlobalIdx, elemLocalIdx, localFaceIdx);

        return factor;
    }

    Scalar getRealFluxFaceArea(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][interactionVolume.getIndexOnElement(elemLocalIdx, localFaceIdx)][fluxFaceArea];

        return factor;
    }

    Scalar getRealFaceArea(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][interactionVolume.getIndexOnElement(elemLocalIdx, localFaceIdx)][realFaceArea];

        return factor;
    }

    FvMpfaL3dInteractionVolumeContainer(Problem& problem) :
        problem_(problem)
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }
    }

protected:
    void addRealFluxFaceArea_(Scalar faceArea, int globalIdx, int faceIdx)
    {
        realFluxFaceArea_[globalIdx][faceIdx][fluxFaceArea] += faceArea;
    }
    void addRealFaceArea_(Scalar faceArea, int globalIdx, int faceIdx)
    {
        realFluxFaceArea_[globalIdx][faceIdx][realFaceArea] += faceArea;
    }

    Problem& problem_;

    GlobalInteractionVolumeVector interactionVolumes_;
    FaceAreaVector realFluxFaceArea_;
private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements(const Element& element, std::vector < std::vector<int> >& elemVertMap)
{
    int globalIdx = problem_.variables().index(element);

    int globalVertIdx = problem_.variables().vertexMapper().map(element, 0, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 7);
    elemVertMap[globalVertIdx][7] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 1, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 6);
    elemVertMap[globalVertIdx][6] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 2, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 5);
    elemVertMap[globalVertIdx][5] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 3, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 4);
    elemVertMap[globalVertIdx][4] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 4, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 3);
    elemVertMap[globalVertIdx][3] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 5, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 2);
    elemVertMap[globalVertIdx][2] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 6, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 1);
    elemVertMap[globalVertIdx][1] = globalIdx;

    globalVertIdx = problem_.variables().vertexMapper().map(element, 7, dim);
    interactionVolumes_[globalVertIdx].setSubVolumeElement(element, 0);
    elemVertMap[globalVertIdx][0] = globalIdx;
}

template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeIntersectionInfo(const Element& element, std::vector < std::vector<int> >& elemVertMap)
{
    BoundaryTypes bcType;

    int globalIdx = problem_.variables().index(element);

    const ElementGeometry& geometry = element.geometry();

    const ReferenceElement& referenceElement = ReferenceElements::general(geometry.type());

    int levelI = element.level();

    // run through all intersections
    IntersectionIterator isIt = problem_.gridView().ibegin(element);
    IntersectionIterator isItEnd = problem_.gridView().iend(element);
    for (; isIt != isItEnd; ++isIt)
    {
        int indexInInside = isIt->indexInInside();

        DimVector normal = isIt->centerUnitOuterNormal();

        const IntersectionGeometry& isGeometry = isIt->geometry();

        Scalar faceVol = isGeometry.volume();

        const DimVector& globalPosFace = isGeometry.center();

        bool takeIntersection = true;
        if (isIt->neighbor())
        {
            ElementPointer outside = isIt->outside();
            int globalIdxJ = problem_.variables().index(*outside);

            if (levelI == outside.level() && globalIdx > globalIdxJ)
                takeIntersection = false;
            if (levelI < outside.level())
                takeIntersection = false;
        }

        if (takeIntersection)
        {
            addRealFaceArea_(faceVol, globalIdx, indexInInside);
            if (isIt->neighbor())
            {
                int globalIdxJ = problem_.variables().index(*(isIt->outside()));
                addRealFaceArea_(faceVol, globalIdxJ, isIt->indexInOutside());
            }

            for (int i = 0; i < isGeometry.corners(); i++)
            {
                int localVertIdx = referenceElement.subEntity(indexInInside, 1, i, dim);

                int globalVertIdx = problem_.variables().vertexMapper().map(element, localVertIdx, dim);

                InteractionVolume& interactionVolume = interactionVolumes_[globalVertIdx];

                if (elemVertMap[globalVertIdx][0] == globalIdx)
                {
                    if (indexInInside == 1)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 0, 0);
                        interactionVolume.setNormal(normal, 0, 0);
                        interactionVolume.setFacePosition(globalPosFace, 0);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 1, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 1, 1);
                        }
                    }
                    else if (indexInInside == 3)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 0, 1);
                        interactionVolume.setNormal(normal, 0, 1);
                        interactionVolume.setFacePosition(globalPosFace, 3);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 2, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 2, 0);
                        }
                    }
                    else if (indexInInside == 5)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 0, 2);
                        interactionVolume.setNormal(normal, 0, 2);
                        interactionVolume.setFacePosition(globalPosFace, 8);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 4, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 4, 0);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][1] == globalIdx)
                {
                    if (indexInInside == 3)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 1, 0);
                        interactionVolume.setNormal(normal, 1, 0);
                        interactionVolume.setFacePosition(globalPosFace, 1);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 3, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 3, 1);
                        }
                    }
                    else if (indexInInside == 0)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 1, 1);
                        interactionVolume.setNormal(normal, 1, 1);
                        interactionVolume.setFacePosition(globalPosFace, 0);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 0, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 0, 0);
                        }
                    }
                    else if (indexInInside == 5)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 1, 2);
                        interactionVolume.setNormal(normal, 1, 2);
                        interactionVolume.setFacePosition(globalPosFace, 9);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 5, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 5, 0);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][2] == globalIdx)
                {
                    if (indexInInside == 2)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 2, 0);
                        interactionVolume.setNormal(normal, 2, 0);
                        interactionVolume.setFacePosition(globalPosFace, 3);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 0, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 0, 1);
                        }
                    }
                    else if (indexInInside == 1)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 2, 1);
                        interactionVolume.setNormal(normal, 2, 1);
                        interactionVolume.setFacePosition(globalPosFace, 2);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 3, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 3, 0);
                        }
                    }
                    else if (indexInInside == 5)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 2, 2);
                        interactionVolume.setNormal(normal, 2, 2);
                        interactionVolume.setFacePosition(globalPosFace, 11);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 6, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 6, 0);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][3] == globalIdx)
                {
                    if (indexInInside == 0)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 3, 0);
                        interactionVolume.setNormal(normal, 3, 0);
                        interactionVolume.setFacePosition(globalPosFace, 2);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 2, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 2, 1);
                        }
                    }
                    else if (indexInInside == 2)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 3, 1);
                        interactionVolume.setNormal(normal, 3, 1);
                        interactionVolume.setFacePosition(globalPosFace, 1);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 1, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 1, 0);
                        }
                    }
                    else if (indexInInside == 5)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 3, 2);
                        interactionVolume.setNormal(normal, 3, 2);
                        interactionVolume.setFacePosition(globalPosFace, 10);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 7, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 7, 0);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][4] == globalIdx)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 4, 0);
                        interactionVolume.setNormal(normal, 4, 0);
                        interactionVolume.setFacePosition(globalPosFace, 8);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 0, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 0, 2);
                        }
                    }
                    else if (indexInInside == 1)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 4, 1);
                        interactionVolume.setNormal(normal, 4, 1);
                        interactionVolume.setFacePosition(globalPosFace, 4);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 5, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 5, 2);
                        }
                    }
                    else if (indexInInside == 3)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 4, 2);
                        interactionVolume.setNormal(normal, 4, 2);
                        interactionVolume.setFacePosition(globalPosFace, 7);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 6, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 6, 1);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][5] == globalIdx)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 5, 0);
                        interactionVolume.setNormal(normal, 5, 0);
                        interactionVolume.setFacePosition(globalPosFace, 9);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 1, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 1, 2);
                        }
                    }
                    else if (indexInInside == 3)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 5, 1);
                        interactionVolume.setNormal(normal, 5, 1);
                        interactionVolume.setFacePosition(globalPosFace, 5);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 7, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 7, 2);
                        }

                    }
                    else if (indexInInside == 0)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 5, 2);
                        interactionVolume.setNormal(normal, 5, 2);
                        interactionVolume.setFacePosition(globalPosFace, 4);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 4, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 4, 1);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][6] == globalIdx)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 6, 0);
                        interactionVolume.setNormal(normal, 6, 0);
                        interactionVolume.setFacePosition(globalPosFace, 11);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 2, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 2, 2);
                        }
                    }
                    else if (indexInInside == 2)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 6, 1);
                        interactionVolume.setNormal(normal, 6, 1);
                        interactionVolume.setFacePosition(globalPosFace, 7);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 4, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 4, 2);
                        }
                    }
                    else if (indexInInside == 1)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 6, 2);
                        interactionVolume.setNormal(normal, 6, 2);
                        interactionVolume.setFacePosition(globalPosFace, 6);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 7, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 7, 1);
                        }
                    }
                }
                if (elemVertMap[globalVertIdx][7] == globalIdx)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 7, 0);
                        interactionVolume.setNormal(normal, 7, 0);
                        interactionVolume.setFacePosition(globalPosFace, 10);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 3, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 3, 2);
                        }
                    }
                    else if (indexInInside == 0)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 7, 1);
                        interactionVolume.setNormal(normal, 7, 1);
                        interactionVolume.setFacePosition(globalPosFace, 6);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 6, 2);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 6, 2);
                        }
                    }
                    else if (indexInInside == 2)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 7, 2);
                        interactionVolume.setNormal(normal, 7, 2);
                        interactionVolume.setFacePosition(globalPosFace, 5);

                        if (isIt->neighbor())
                        {
                            int indexInOutside = isIt->indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 5, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 5, 1);
                        }
                    }
                }
                if (isIt->boundary())
                {
                    if (elemVertMap[globalVertIdx][0] == globalIdx)
                    {
                        if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 0);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 0);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 0);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 3);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 3);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 3);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 8);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 8);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 8);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][1] == globalIdx)
                    {
                        if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 1);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 1);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 1);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 0);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 0);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 0);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 9);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 9);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 9);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][2] == globalIdx)
                    {
                        if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 3);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 3);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 3);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 2);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 2);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 2);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 11);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 11);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 11);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][3] == globalIdx)
                    {
                        if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 2);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 2);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 2);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 1);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 1);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 1);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 10);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 10);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 10);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][4] == globalIdx)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 8);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 8);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 8);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 4);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 4);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 4);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 7);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 7);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 7);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][5] == globalIdx)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 9);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 9);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 9);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 5);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 5);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 5);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 4);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 4);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 4);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][6] == globalIdx)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 11);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 11);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 11);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 7);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 7);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 7);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 6);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 6);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 6);
                            }
                        }
                    }
                    if (elemVertMap[globalVertIdx][7] == globalIdx)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 10);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 10);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 10);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 6);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 6);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 6);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, *isIt);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 5);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, *isIt);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 5);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, *isIt);
                                interactionVolume.setDirichletCondition(boundValues, 5);
                            }
                        }
                    }
                }
            }
        }

    }
}

template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeInnerInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex, bool sameLevel)
{
    const DimVector& centerPos = vertex.geometry().center();

    interactionVolume.setCenterPosition(centerPos);

    if (sameLevel)
    {
        ElementPointer& elementPointer1 = interactionVolume.getSubVolumeElement(0);
        ElementPointer& elementPointer8 = interactionVolume.getSubVolumeElement(7);

        const ElementGeometry& geometry1 = elementPointer1->geometry();
        const ElementGeometry& geometry8 = elementPointer8->geometry();

        const ReferenceElement& referenceElement = ReferenceElements::general(geometry1.type());

        DimVector edgeCoord(geometry1.global(referenceElement.position(9, dim - 1)));
        interactionVolume.setEdgePosition(edgeCoord, 2);
        edgeCoord = geometry1.global(referenceElement.position(3, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 0);
        edgeCoord = geometry1.global(referenceElement.position(11, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 5);

        edgeCoord = geometry8.global(referenceElement.position(4, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 4);
        edgeCoord = geometry8.global(referenceElement.position(6, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 3);
        edgeCoord = geometry8.global(referenceElement.position(0, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 1);
    }

    DimVector edgeCoord1(interactionVolume.getEdgePosition(0));
    DimVector edgeCoord2(interactionVolume.getEdgePosition(1));
    DimVector edgeCoord3(interactionVolume.getEdgePosition(2));
    DimVector edgeCoord4(interactionVolume.getEdgePosition(3));
    DimVector edgeCoord5(interactionVolume.getEdgePosition(4));
    DimVector edgeCoord6(interactionVolume.getEdgePosition(5));

    DimVector crossProductVector1(0);
    DimVector crossProductVector2(0);

    GlobalPosition globalPosFace = interactionVolume.getFacePosition(0);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord1 - edgeCoord3;
    Scalar faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 0);

    globalPosFace = interactionVolume.getFacePosition(1);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord1 - edgeCoord4;

    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;

    interactionVolume.setFaceArea(faceArea, 1);

    globalPosFace = interactionVolume.getFacePosition(2);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord1 - edgeCoord5;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 2);

    globalPosFace = interactionVolume.getFacePosition(3);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord1 - edgeCoord6;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 3);

    globalPosFace = interactionVolume.getFacePosition(4);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord2 - edgeCoord3;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 4);

    globalPosFace = interactionVolume.getFacePosition(5);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord2 - edgeCoord4;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 5);

    globalPosFace = interactionVolume.getFacePosition(6);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord2 - edgeCoord5;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 6);

    globalPosFace = interactionVolume.getFacePosition(7);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord2 - edgeCoord6;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 7);

    globalPosFace = interactionVolume.getFacePosition(8);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord3 - edgeCoord6;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 8);

    globalPosFace = interactionVolume.getFacePosition(9);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord3 - edgeCoord4;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 9);

    globalPosFace = interactionVolume.getFacePosition(10);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord4 - edgeCoord5;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 10);

    globalPosFace = interactionVolume.getFacePosition(11);
    crossProductVector1 = centerPos - globalPosFace;
    crossProductVector2 = edgeCoord5 - edgeCoord6;
    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
    interactionVolume.setFaceArea(faceArea, 11);
}

template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeBoundaryInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex)
{
    const DimVector& centerPos = vertex.geometry().center();

    interactionVolume.setCenterPosition(centerPos);

    //corner
    switch (interactionVolume.getElementNumber())
    {
    case 1:
        {
            if (interactionVolume.hasSubVolumeElement(0))
            {
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(10);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(0);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 0, 0)/4.0,0);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 0, 1)/4.0,3);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 0, 2)/4.0,8);
            }
            if (interactionVolume.hasSubVolumeElement(1))
            {
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(10);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(1);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 1, 0)/4.0,1);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 1, 1)/4.0,0);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 1, 2)/4.0,9);
            }
            if (interactionVolume.hasSubVolumeElement(2))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(10);

                ElementPointer& element = interactionVolume.getSubVolumeElement(2);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 2, 0)/4.0,3);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 2, 1)/4.0,2);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 2, 2)/4.0,11);
            }
            if (interactionVolume.hasSubVolumeElement(3))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(3);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 3, 0)/4.0,2);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 3, 1)/4.0,1);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 3, 2)/4.0,10);
            }
            if (interactionVolume.hasSubVolumeElement(4))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(10);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(4);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 4, 0)/4.0,8);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 4, 1)/4.0,4);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 4, 2)/4.0,7);
            }
            if (interactionVolume.hasSubVolumeElement(5))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(6);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(10);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(5);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 5, 0)/4.0,9);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 5, 1)/4.0,5);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 5, 2)/4.0,4);
            }
            if (interactionVolume.hasSubVolumeElement(6))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(5);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(10);

                ElementPointer& element = interactionVolume.getSubVolumeElement(6);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 6, 0)/4.0,11);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 6, 1)/4.0,7);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 6, 2)/4.0,6);
            }
            if (interactionVolume.hasSubVolumeElement(7))
            {
                interactionVolume.setOutsideFace(0);
                interactionVolume.setOutsideFace(1);
                interactionVolume.setOutsideFace(2);
                interactionVolume.setOutsideFace(3);
                interactionVolume.setOutsideFace(4);
                interactionVolume.setOutsideFace(7);
                interactionVolume.setOutsideFace(8);
                interactionVolume.setOutsideFace(9);
                interactionVolume.setOutsideFace(11);

                ElementPointer& element = interactionVolume.getSubVolumeElement(7);
                int globalIdx = problem_.variables().index(*element);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 7, 0)/4.0,10);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 7, 1)/4.0,6);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx, 7, 2)/4.0,5);
            }
            break;
        }
    case 2:
        {
            // edge
            if (interactionVolume.hasSubVolumeElement(0))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 1)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 2)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 0)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 2)/4.0,9);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(2))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(2);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 0)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 2)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 1)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 2)/4.0,11);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(4))
                {
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(4);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 0)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 1)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 1)/4.0,4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 2)/4.0,7);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(7))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(7);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(5))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(5);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 0)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 1)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 5, 0)/4.0,9);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 5, 2)/4.0,4);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(6))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(9);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(6);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 0)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 2)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 6, 0)/4.0,11);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 6, 1)/4.0,7);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(3))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(3);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 1)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 2)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 3, 0)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 3, 1)/4.0,1);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(5))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(5);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(10);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 5, 1)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 5, 2)/4.0,4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 0)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 1)/4.0,0);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(4))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(4);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 5, 0)/4.0,9);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 5, 1)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 0)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 2)/4.0,7);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(6))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(6);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(4))
                {
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(4);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 6, 0)/4.0,11);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 6, 2)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 0)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 4, 1)/4.0,4);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(2))
                {
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(8);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(2);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 6, 1)/4.0,7);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 6, 2)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 0)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 1)/4.0,2);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(3))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(3);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 3, 0)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 3, 2)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 1)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 2)/4.0,9);

                    return;
                }
                if (interactionVolume.hasSubVolumeElement(2))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(9);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(2);
                    int globalIdx2 = problem_.variables().index(*element2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 3, 1)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 3, 2)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 0)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 2)/4.0,11);

                    return;
                }
            }
            break;
        }
    case 4:
        {
            //side
            if (interactionVolume.hasSubVolumeElement(0))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
                    int globalIdx2 = problem_.variables().index(*element2);
                    if (interactionVolume.hasSubVolumeElement(2) && interactionVolume.hasSubVolumeElement(3))
                    {
                        interactionVolume.setOutsideFace(4);
                        interactionVolume.setOutsideFace(5);
                        interactionVolume.setOutsideFace(6);
                        interactionVolume.setOutsideFace(7);

                        ElementPointer& element3 = interactionVolume.getSubVolumeElement(2);
                        int globalIdx3 = problem_.variables().index(*element3);
                        ElementPointer& element4 = interactionVolume.getSubVolumeElement(3);
                        int globalIdx4 = problem_.variables().index(*element4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 2)/4.0, 8);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 2)/4.0, 9);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 2, 2)/4.0, 10);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 3, 2)/4.0, 11);

                        return;
                    }
                    if (interactionVolume.hasSubVolumeElement(4) && interactionVolume.hasSubVolumeElement(5))
                    {
                        interactionVolume.setOutsideFace(2);
                        interactionVolume.setOutsideFace(6);
                        interactionVolume.setOutsideFace(10);
                        interactionVolume.setOutsideFace(11);

                        ElementPointer& element3 = interactionVolume.getSubVolumeElement(4);
                        int globalIdx3 = problem_.variables().index(*element3);
                        ElementPointer& element4 = interactionVolume.getSubVolumeElement(5);
                        int globalIdx4 = problem_.variables().index(*element4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 1)/4.0, 3);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 1, 0)/4.0, 1);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 4, 2)/4.0, 7);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 5, 1)/4.0, 5);

                        return;
                    }
                }
                if (interactionVolume.hasSubVolumeElement(2) && interactionVolume.hasSubVolumeElement(4)
                    && interactionVolume.hasSubVolumeElement(6))
                {
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(2);
                    int globalIdx2 = problem_.variables().index(*element2);
                    ElementPointer& element3 = interactionVolume.getSubVolumeElement(4);
                    int globalIdx3 = problem_.variables().index(*element3);
                    ElementPointer& element4 = interactionVolume.getSubVolumeElement(6);
                    int globalIdx4 = problem_.variables().index(*element4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 0, 0)/4.0, 0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 2, 1)/4.0, 2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 4, 1)/4.0, 4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 6, 2)/4.0, 6);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(7))
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(7);
                int globalIdx1 = problem_.variables().index(*element1);
                if (interactionVolume.hasSubVolumeElement(5))
                {
                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(5);
                    int globalIdx2 = problem_.variables().index(*element2);
                    if (interactionVolume.hasSubVolumeElement(1) && interactionVolume.hasSubVolumeElement(3))
                    {
                        interactionVolume.setOutsideFace(3);
                        interactionVolume.setOutsideFace(7);
                        interactionVolume.setOutsideFace(8);
                        interactionVolume.setOutsideFace(11);

                        ElementPointer& element3 = interactionVolume.getSubVolumeElement(1);
                        int globalIdx3 = problem_.variables().index(*element3);
                        ElementPointer& element4 = interactionVolume.getSubVolumeElement(3);
                        int globalIdx4 = problem_.variables().index(*element4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 1)/4.0, 6);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 5, 2)/4.0, 4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 1, 1)/4.0, 0);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 3, 0)/4.0, 2);

                        return;
                    }
                    if (interactionVolume.hasSubVolumeElement(4) && interactionVolume.hasSubVolumeElement(6))
                    {
                        interactionVolume.setOutsideFace(0);
                        interactionVolume.setOutsideFace(1);
                        interactionVolume.setOutsideFace(2);
                        interactionVolume.setOutsideFace(3);

                        ElementPointer& element3 = interactionVolume.getSubVolumeElement(4);
                        int globalIdx3 = problem_.variables().index(*element3);
                        ElementPointer& element4 = interactionVolume.getSubVolumeElement(6);
                        int globalIdx4 = problem_.variables().index(*element4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 0)/4.0, 10);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 5, 0)/4.0, 9);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 4, 0)/4.0, 8);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 6, 0)/4.0, 11);

                        return;
                    }
                }
                if (interactionVolume.hasSubVolumeElement(6) && interactionVolume.hasSubVolumeElement(2)
                    && interactionVolume.hasSubVolumeElement(3))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(9);

                    ElementPointer& element2 = interactionVolume.getSubVolumeElement(6);
                    int globalIdx2 = problem_.variables().index(*element2);
                    ElementPointer& element3 = interactionVolume.getSubVolumeElement(2);
                    int globalIdx3 = problem_.variables().index(*element3);
                    ElementPointer& element4 = interactionVolume.getSubVolumeElement(3);
                    int globalIdx4 = problem_.variables().index(*element4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx1, 7, 2)/4.0, 5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx2, 6, 1)/4.0, 7);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx3, 2, 0)/4.0, 3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, globalIdx4, 3, 1)/4.0, 1);

                    return;
                }
            }
            break;
        }
    default:
        DUNE_THROW(Dune::NotImplemented, "Boundary shape not implemented");
        break;
    }
}


// only for 3-D general quadrilateral
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeInteractionVolumeInfo()
{
    std::vector < std::vector<int> > elemVertMap(problem_.gridView().size(dim), std::vector<int>(8, -1));

    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        ElementPointer ePtr = *eIt;
        storeSubVolumeElements(*ePtr, elemVertMap);
    }

    for  (int i = 0; i < interactionVolumes_.size(); i++)
        if (interactionVolumes_[i].getElementNumber() == 0)
            interactionVolumes_[i].printInteractionVolumeInfo();

    // run through all elements

    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get common geometry information for the following computation

        ElementPointer ePtr = *eIt;
        storeIntersectionInfo(*ePtr, elemVertMap);
    }

    // run through all vertices
    VertexIterator vItEnd = problem_.gridView().template end<dim>();
    for (VertexIterator vIt = problem_.gridView().template begin<dim>(); vIt != vItEnd; ++vIt)
    {
        int globalVertIdx = problem_.variables().index(*vIt);

        InteractionVolume& interactionVolume = interactionVolumes_[globalVertIdx];

        if (interactionVolume.getElementNumber() == 8)
        {
            storeInnerInteractionVolume(interactionVolume, *vIt);
        }
        else if (interactionVolume.isBoundaryInteractionVolume())
        {
            storeBoundaryInteractionVolume(interactionVolume, *vIt);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,"Interaction volume is no boundary volume but consists of less than 8 elements");
        }

        if (!interactionVolume.isBoundaryInteractionVolume())
        {
            ElementPointer& elementPointer1 = interactionVolume.getSubVolumeElement(0);
            ElementPointer& elementPointer2 = interactionVolume.getSubVolumeElement(1);
            ElementPointer& elementPointer3 = interactionVolume.getSubVolumeElement(2);
            ElementPointer& elementPointer4 = interactionVolume.getSubVolumeElement(3);
            ElementPointer& elementPointer5 = interactionVolume.getSubVolumeElement(4);
            ElementPointer& elementPointer6 = interactionVolume.getSubVolumeElement(5);
            ElementPointer& elementPointer7 = interactionVolume.getSubVolumeElement(6);
            ElementPointer& elementPointer8 = interactionVolume.getSubVolumeElement(7);

            int globalIdx1 = problem_.variables().index(*elementPointer1);
            int globalIdx2 = problem_.variables().index(*elementPointer2);
            int globalIdx3 = problem_.variables().index(*elementPointer3);
            int globalIdx4 = problem_.variables().index(*elementPointer4);
            int globalIdx5 = problem_.variables().index(*elementPointer5);
            int globalIdx6 = problem_.variables().index(*elementPointer6);
            int globalIdx7 = problem_.variables().index(*elementPointer7);
            int globalIdx8 = problem_.variables().index(*elementPointer8);

            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 0), globalIdx1, interactionVolume.getIndexOnElement(0, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 1), globalIdx1, interactionVolume.getIndexOnElement(0, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 2), globalIdx1, interactionVolume.getIndexOnElement(0, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 0), globalIdx2, interactionVolume.getIndexOnElement(1, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 1), globalIdx2, interactionVolume.getIndexOnElement(1, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 2), globalIdx2, interactionVolume.getIndexOnElement(1, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 0), globalIdx3, interactionVolume.getIndexOnElement(2, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 1), globalIdx3, interactionVolume.getIndexOnElement(2, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 2), globalIdx3, interactionVolume.getIndexOnElement(2, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 0), globalIdx4, interactionVolume.getIndexOnElement(3, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 1), globalIdx4, interactionVolume.getIndexOnElement(3, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 2), globalIdx4, interactionVolume.getIndexOnElement(3, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 0), globalIdx5, interactionVolume.getIndexOnElement(4, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 1), globalIdx5, interactionVolume.getIndexOnElement(4, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 2), globalIdx5, interactionVolume.getIndexOnElement(4, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 0), globalIdx6, interactionVolume.getIndexOnElement(5, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 1), globalIdx6, interactionVolume.getIndexOnElement(5, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 2), globalIdx6, interactionVolume.getIndexOnElement(5, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 0), globalIdx7, interactionVolume.getIndexOnElement(6, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 1), globalIdx7, interactionVolume.getIndexOnElement(6, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 2), globalIdx7, interactionVolume.getIndexOnElement(6, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 0), globalIdx8, interactionVolume.getIndexOnElement(7, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 1), globalIdx8, interactionVolume.getIndexOnElement(7, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 2), globalIdx8, interactionVolume.getIndexOnElement(7, 2));
        }

    }
}
} // end of Dune namespace
#endif
