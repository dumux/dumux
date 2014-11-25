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
#ifndef DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_ADAPTIVE_HH
#define DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_ADAPTIVE_HH

// dumux environment
#include <dumux/decoupled/common/pressureproperties.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include <dumux/decoupled/common/fv/mpfa/mpfalinteractionvolume3dadaptive.hh>
#include "fvmpfal3dinteractionvolumecontainer.hh"

/**
 * @file
 * @brief  Interactionvolume container for 3-d MPFA L-method on an h-adaptive grid
 */

namespace Dumux
{

//! \ingroup FVPressure2p mpfa
/*! \brief Interactionvolume container for 3-d MPFA L-method on an h-adaptive grid
 *
 * Container class which stores MPFA-interaction-volume information for each vertex of a DUNE grid.
 * Each <tt>InteractionVolume</tt> object stores the information which is necessary to calculate MPFA transmissibility matrices:
 *
 * - relationship and orientation of the elements around a vertex (see doc/docextra/3dmpfa)
 * - geometric information, such as element/face/edge positions, normals, ...
 *
 */
template<class TypeTag>
class FvMpfaL3dInteractionVolumeContainerAdaptive: public FvMpfaL3dInteractionVolumeContainer<TypeTag>
{
    friend class FvMpfaL3dInteractionVolumeContainer<TypeTag>;
    typedef FvMpfaL3dInteractionVolumeContainer<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dim>::Entity Vertex;
    typedef typename Element::Geometry ElementGeometry;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef std::vector<Dune::FieldVector<std::set<int>, 2*dim> > FaceVerticesVector;

    enum
        {
            pressureEqIdx = Indices::pressureEqIdx,
        };

    typedef IndexTranslatorAdaptive IndexTranslator;

public:
    //! Type for storing an MPFA-interaction-volume.
    //! (Usually of type Dumux::FvMpfaL3dInteractionVolume or Dumux::FvMpfaL3dInteractionVolumeAdaptive)
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolume) InteractionVolume;

private:

    void storeHangingNodeInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex);
    void storeInnerInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex);

protected:
    void storeInteractionVolumeInfo();

public:

    /*! \brief Returns the set of vertices on an element face
     *
     * The DUNE reference elements does not allow to access hanging nodes from a given element face.
     * However, if a flux through a entire element face has to be calculated, e.g. if single fluxes
     * have to be updated in an implicit treatment of the transport equation, it is necessary to get
     * the complete set of vertices on a face: 4 corners + all hanging nodes.
     */
    std::set<int>& faceVerticeIndices(int eIdxGlobal, int fIdx)
    {
        return faceVertices_[eIdxGlobal][fIdx];
    }

    //! Constructs a FvMpfaL3dInteractionVolumeContainerAdaptive object
    /**
     * \param problem A problem class object
     */
    FvMpfaL3dInteractionVolumeContainerAdaptive(Problem& problem) :
        ParentType(problem), problem_(problem)
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }
    }

private:
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolumeContainer) Implementation;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

    Problem& problem_;
    FaceVerticesVector faceVertices_;
};

/*! \brief Stores additional information which can be constructed for interaction
 *         volumes of non-boundary vertices.
 *
 * Stores additional information which can be constructed for interaction volumes of
 * non-boundary vertices:
 *
 *  - edge coordinates (coordinates of edge-continuity-points)
 *  - flux face areas
 *
 * Assumes a local storage following the DUNE reference element index, which is
 * performed by the function
 * Dumux::FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements
 * (const Element& element, std::vector < std::vector<int> >& elemVertMap).
 *
 * In the case of an adaptive grids with hanging nodes it is important to notice,
 * that the smaller cell of two intersecting cells of different grid level determines
 * the geometric information of the flux face (size, position of the center, etc.).
 * This has to be stored correctly for both sides (elements) of an intersection.
 *
 * \param interactionVolume An interaction volume object
 * \param vertex The vertex (level dim entity) for which the interaction volume is stored
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainerAdaptive<TypeTag>::storeInnerInteractionVolume(InteractionVolume& interactionVolume,
                                                                                       const Vertex& vertex)
{
    // if cells of different grid level appear, the geometric information is constructed going from the
    // coarsest level to the finest level. This ensures that coarser information ins always overwritten by
    // finer information.
    if (!interactionVolume.sameLevel())
    {
        // sort the local element indices according to the grid level
        std::vector < std::vector<int> > levelIdx(8, std::vector<int>(2));
        for (int i = 0; i < 8; i++)
        {
            levelIdx[i][0] = i;
            levelIdx[i][1] = interactionVolume.getSubVolumeElement(i)->level();
        }

        std::sort(levelIdx.begin(), levelIdx.end(), sort_compare);

        // Generate and store the geometric information going from the coarsest to the finest level.
        // For the calculation we take advantage from the fact that the ordering inside the interaction volume
        // with respect to the DUNE reference element is known due to the storage process of the elements in
        // Dumux::FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements
        //        (const Element& element, std::vector < std::vector<int> >& elemVertMap)
        for (int i = 0; i < 8; i++)
        {
            int idx = levelIdx[i][0];
            ElementPointer& elementPointer = interactionVolume.getSubVolumeElement(idx);

            const ElementGeometry& geometry = elementPointer->geometry();

            const ReferenceElement& referenceElement = ReferenceElements::general(geometry.type());

            switch (idx)
            {
            case 0:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(9, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 2);
                    edgeCoord = geometry.global(referenceElement.position(3, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 0);
                    edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 5);

                    break;
                }
            case 1:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(2, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 0);
                    edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 2);
                    edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 3);

                    break;
                }
            case 2:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(1, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 0);
                    edgeCoord = geometry.global(referenceElement.position(9, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 4);
                    edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 5);

                    break;
                }
            case 3:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(0, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 0);
                    edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 4);
                    edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 3);

                    break;
                }
            case 4:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(3, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 1);
                    edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 2);
                    edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 5);

                    break;
                }
            case 5:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(2, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 1);
                    edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 2);
                    edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 3);

                    break;
                }
            case 6:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(1, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 1);
                    edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 4);
                    edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 5);

                    break;
                }
            case 7:
                {
                    DimVector edgeCoord(geometry.global(referenceElement.position(0, dim - 1)));
                    interactionVolume.setEdgePosition(edgeCoord, 1);
                    edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 4);
                    edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                    interactionVolume.setEdgePosition(edgeCoord, 3);

                    break;
                }
            }
        }
        ParentType::storeInnerInteractionVolume(interactionVolume, vertex, false);
    }
    else
    {
        ParentType::storeInnerInteractionVolume(interactionVolume, vertex, true);
    }
}

/*! \brief Stores additional information which can be constructed for interaction volumes around hanging nodes.
 *
 *  Stores additional information which can be constructed for interaction volumes around hanging nodes.
 *
 *  - missing cells: As hanging nodes cannot be accessed from a cell face using the DUNE reference element,
 *                   the attached coarser cells do not appear in the interaction volume object after
 *                   execution of the function Dumux::FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements
 *                                             (const Element& element, std::vector < std::vector<int> >& elemVertMap).
 *                   We take advantage of this fact because it allows us to identify the type of hanging-node interaction volume.
 *                   If, for example, only two cells are stored, we know that the interaction volume is of type 5 according to
 *                   Wolff 2013: http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/
 *                   As first step, the orientation of the stored elements in the local interaction volume indices
 *                   is rotated such that cells 1 and 2 according to the local index (see doc/docextra/3dmpfa)
 *                   are always of the smallest existing level. Thus, the relationship between the grid elements
 *                   in a hanging-node-interaction-volume is always known if the type of interaction volume is known
 *                   (1-7, see Wolff 2013: http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/).
 *                   As a second step, the missing cells are added.
 *  - edge coordinates (coordinates of edge-continuity-points)
 *  - flux face areas
 *
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainerAdaptive<TypeTag>::storeHangingNodeInteractionVolume(InteractionVolume& interactionVolume,
                                                                                             const Vertex& vertex)
{
    const DimVector& centerPos = vertex.geometry().center();

    interactionVolume.setCenterPosition(centerPos);

    // sort the local element indices according to the grid level
    std::vector < std::vector<int> > levelIdx(8, std::vector<int>(2));
    for (int i = 0; i < 8; i++)
    {
        levelIdx[i][0] = i;
        if (interactionVolume.hasSubVolumeElement(i))
            levelIdx[i][1] = interactionVolume.getSubVolumeElement(i)->level();
        else
            levelIdx[i][1] = -1;
    }

    std::sort(levelIdx.begin(), levelIdx.end(), sort_compare);

    // Generate and store the geometric information going from the coarsest to the finest level.
    // For the calculation we take advantage from the fact that the ordering inside the interaction volume
    // with respect to the DUNE reference element is known due to the storage process of the elements in
    // Dumux::FvMpfaL3dInteractionVolumeContainer<TypeTag>::
    // storeSubVolumeElements(const Element& element, std::vector < std::vector<int> >& elemVertMap)
    for (int i = 0; i < 8; i++)
    {
        if (levelIdx[i][1] < 0)
            continue;

        int idx = levelIdx[i][0];

        ElementPointer& elementPointer = interactionVolume.getSubVolumeElement(idx);

        const ElementGeometry& geometry = elementPointer->geometry();

        const ReferenceElement& referenceElement = ReferenceElements::general(geometry.type());

        switch (idx)
        {
        case 0:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(9, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 2);
                edgeCoord = geometry.global(referenceElement.position(3, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 0);
                edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 5);

                break;
            }
        case 1:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(2, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 0);
                edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 2);
                edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 3);

                break;
            }
        case 2:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(1, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 0);
                edgeCoord = geometry.global(referenceElement.position(9, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 4);
                edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 5);

                break;
            }
        case 3:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(0, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 0);
                edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 4);
                edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 3);

                break;
            }
        case 4:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(3, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 1);
                edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 2);
                edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 5);

                break;
            }
        case 5:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(2, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 1);
                edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 2);
                edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 3);

                break;
            }
        case 6:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(1, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 1);
                edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 4);
                edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 5);

                break;
            }
        case 7:
            {
                DimVector edgeCoord(geometry.global(referenceElement.position(0, dim - 1)));
                interactionVolume.setEdgePosition(edgeCoord, 1);
                edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 4);
                edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                interactionVolume.setEdgePosition(edgeCoord, 3);

                break;
            }
        }
    }

    // Choose the type of the hanging-node-interaction-volume depending on the number of stored fine
    // elements (see dissertation M. Wolff, http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/)
    // and add missing elements and geometric information
    switch (interactionVolume.getElementNumber())
    {
    // hanging-node interaction volume of type 5 or 7
    case 2:
        {
            InteractionVolume hangingNodeVolume;

            std::vector<int> elemIdxOld;
            for (int i = 0; i < InteractionVolume::subVolumeTotalNum; i++)
            {
                if (interactionVolume.hasSubVolumeElement(i))
                    elemIdxOld.push_back(i);
            }

            int zeroFaceIdx = IndexTranslator::getFaceIndexFromElements(elemIdxOld[0], elemIdxOld[1]);

            std::vector<int> elemIdxNew(2);
            elemIdxNew[0] = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elemIdxOld[0]);
            elemIdxNew[1] = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elemIdxOld[1]);


            hangingNodeVolume.setSubVolumeElement(interactionVolume.getSubVolumeElement(elemIdxOld[0]),
                                                  elemIdxNew[0]);
            hangingNodeVolume.setSubVolumeElement(interactionVolume.getSubVolumeElement(elemIdxOld[1]),
                                                  elemIdxNew[1]);

            for (int elem = 0; elem < InteractionVolume::subVolumeTotalNum; elem++)
            {
                for (int i = 0; i < 3; i++)
                {
                    int faceIdxOld = IndexTranslator::getFaceIndexFromSubVolume(elem, i);
                    int faceIdxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, faceIdxOld);

                    for (int j = 0; j < 3; j++)
                    {
                        int elemNew = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elem);
                        int faceIdxNewTest = IndexTranslator::getFaceIndexFromSubVolume(elemNew, j);

                        if (faceIdxNew == faceIdxNewTest)
                        {
                            hangingNodeVolume.setNormal(interactionVolume.getNormal(elem, i),
                                                        elemNew, j);
                            hangingNodeVolume.setIndexOnElement(
                                                                interactionVolume.getIndexOnElement(elem, i), elemNew, j);
                        }
                    }
                }
            }

            for (int i = 0; i < InteractionVolume::fluxFacesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setFacePosition(interactionVolume.getFacePosition(i), idxNew);
            }

            for (int i = 0; i < InteractionVolume::fluxEdgesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewEdgeIdxFromOldFaceIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setEdgePosition(interactionVolume.getEdgePosition(i), idxNew);
            }

            interactionVolume = hangingNodeVolume;

            interactionVolume.setHangingNodeType(InteractionVolume::twoSmallCells);

            interactionVolume.setCenterPosition(centerPos);

            ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);

            IntersectionIterator isIt = problem_.gridView().ibegin(*element1);
            IntersectionIterator isEndIt = problem_.gridView().iend(*element1);
            for (; isIt != isEndIt; ++isIt)
            {
                int idxInInside = isIt->indexInInside();

                if (idxInInside == interactionVolume.getIndexOnElement(0, 2))
                {
                    if (isIt->neighbor())
                    {
                        ElementPointer outside = isIt->outside();
                        if (element1->level() > outside->level())
                        {
                            interactionVolume.setSubVolumeElement(outside, 4);
                            interactionVolume.setSubVolumeElement(outside, 5);
                        }
                    }
                }
                else if (idxInInside == interactionVolume.getIndexOnElement(0, 1))
                {
                    if (isIt->neighbor())
                    {
                        ElementPointer outside = isIt->outside();
                        if (element1->level() > outside->level())
                        {
                            interactionVolume.setSubVolumeElement(outside, 2);
                            interactionVolume.setSubVolumeElement(outside, 3);
                        }
                    }
                }
            }
            ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
            ElementPointer& element45 = interactionVolume.getSubVolumeElement(4);
            ElementPointer& element23 = interactionVolume.getSubVolumeElement(2);

            IntersectionIterator isIt1 = problem_.gridView().ibegin(*element45);
            IntersectionIterator isIt1End = problem_.gridView().iend(*element45);
            for (; isIt1 != isIt1End; ++isIt1)
            {
                if (isIt1->neighbor())
                {
                    ElementPointer element45Outside = isIt1->outside();

                    IntersectionIterator isIt2 = problem_.gridView().ibegin(*element23);
                    IntersectionIterator isIt2End = problem_.gridView().iend(*element23);
                    for (; isIt2 != isIt2End; ++isIt2)
                    {
                        if (isIt2->neighbor())
                        {
                            ElementPointer element23Outside = isIt2->outside();

                            if (element45Outside == element23Outside && element45Outside != element1
                                && element45Outside != element2)
                            {
                                interactionVolume.setSubVolumeElement(element45Outside, 6);
                                interactionVolume.setSubVolumeElement(element45Outside, 7);
                                DimVector normal = isIt2->centerUnitOuterNormal();
                                interactionVolume.setNormal(normal, 2, 2);
                                interactionVolume.setNormal(normal, 3, 2);
                                normal *= -1;
                                interactionVolume.setNormal(normal, 6, 0);
                                interactionVolume.setNormal(normal, 7, 0);

                                GlobalPosition globalPosFace = isIt2->geometry().center();
                                interactionVolume.setFacePosition(globalPosFace, 10);
                                interactionVolume.setFacePosition(globalPosFace, 11);
                                interactionVolume.setEdgePosition(centerPos, 4);

                                Scalar faceArea = isIt2->geometry().volume()/4.0;

                                interactionVolume.setFaceArea(faceArea, 10);
                                interactionVolume.setFaceArea(faceArea, 11);

                                normal = isIt1->centerUnitOuterNormal();
                                interactionVolume.setNormal(normal, 4, 2);
                                interactionVolume.setNormal(normal, 5, 1);
                                normal *= -1;
                                interactionVolume.setNormal(normal, 6, 1);
                                interactionVolume.setNormal(normal, 7, 2);

                                globalPosFace = isIt1->geometry().center();
                                interactionVolume.setFacePosition(globalPosFace, 5);
                                interactionVolume.setFacePosition(globalPosFace, 7);
                                interactionVolume.setEdgePosition(centerPos, 1);

                                faceArea = isIt1->geometry().volume()/4.0;

                                interactionVolume.setFaceArea(faceArea, 5);
                                interactionVolume.setFaceArea(faceArea, 7);
                            }
                        }
                    }
                }
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

            globalPosFace = interactionVolume.getFacePosition(3);
            crossProductVector1 = centerPos - globalPosFace;
            crossProductVector2 = edgeCoord1 - edgeCoord6;
            faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
            interactionVolume.setFaceArea(faceArea, 3);

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

            break;
        }
        //hanging-node interaction volume of type 1, 3 or 4
    case 4:
        {
            InteractionVolume hangingNodeVolume;

            std::vector<int> elemIdxOld;
            for (int i = 0; i < InteractionVolume::subVolumeTotalNum; i++)
            {
                if (interactionVolume.hasSubVolumeElement(i))
                    elemIdxOld.push_back(i);
            }

            std::set<int> zeroFaceIdxVec;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    int fIdx = IndexTranslator::getFaceIndexFromElements(elemIdxOld[i], elemIdxOld[j]);
                    if (fIdx >= 0)
                        zeroFaceIdxVec.insert(fIdx);
                }
            }

            std::vector<int> elemIdxNew(4);
            int zeroFaceIdx = 0;

            if (zeroFaceIdxVec.size() == 2)
            {
                hangingNodeVolume.setHangingNodeType(InteractionVolume::fourSmallCellsDiag);

                if (zeroFaceIdxVec.find(0) == zeroFaceIdxVec.end())
                    zeroFaceIdx = *zeroFaceIdxVec.begin();

                for (int i = 0; i < 4; i++)
                {
                    elemIdxNew[i] = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elemIdxOld[i]);
                    hangingNodeVolume.setSubVolumeElement(interactionVolume.getSubVolumeElement(elemIdxOld[i]),
                                                          elemIdxNew[i]);
                }
            }
            else if (zeroFaceIdxVec.size() == 4)
            {
                std::set<int>::iterator it = zeroFaceIdxVec.begin();
                for (; it != zeroFaceIdxVec.end(); ++it)
                {
                    zeroFaceIdx = *it;

                    bool isFace = true;
                    for (int i = 0; i < 4; i++)
                    {
                        elemIdxNew[i] = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx,
                                                                                        elemIdxOld[i]);
                        if (elemIdxNew[i] == 4 || elemIdxNew[i] == 5 || elemIdxNew[i] == 6 || elemIdxNew[i] == 7)
                        {
                            isFace = false;
                            break;
                        }
                    }
                    if (isFace)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            hangingNodeVolume.setSubVolumeElement(interactionVolume.getSubVolumeElement(elemIdxOld[i]), elemIdxNew[i]);
                        }

                        break;
                    }
                }
            }
            else
            {
                std::set<int>::iterator it = zeroFaceIdxVec.begin();
                std::cout<<"zeroFaceIdxVec = ";
                for(;it != zeroFaceIdxVec.end(); ++it)
                {
                    std::cout<<" "<<*it;
                }
                std::cout<<"\n";
            }

            for (int elem = 0; elem < InteractionVolume::subVolumeTotalNum; elem++)
            {
                for (int i = 0; i < 3; i++)
                {
                    int faceIdxOld = IndexTranslator::getFaceIndexFromSubVolume(elem, i);
                    int faceIdxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, faceIdxOld);

                    for (int j = 0; j < 3; j++)
                    {
                        int elemNew = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elem);
                        int faceIdxNewTest = IndexTranslator::getFaceIndexFromSubVolume(elemNew, j);

                        if (faceIdxNew == faceIdxNewTest)
                        {
                            hangingNodeVolume.setNormal(interactionVolume.getNormal(elem, i),
                                                        elemNew, j);
                            hangingNodeVolume.setIndexOnElement(interactionVolume.getIndexOnElement(elem, i), elemNew, j);
                        }
                    }
                }
            }

            for (int i = 0; i < InteractionVolume::fluxFacesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setFacePosition(interactionVolume.getFacePosition(i), idxNew);
            }

            for (int i = 0; i < InteractionVolume::fluxEdgesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewEdgeIdxFromOldFaceIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setEdgePosition(interactionVolume.getEdgePosition(i), idxNew);
            }

            interactionVolume = hangingNodeVolume;

            interactionVolume.setCenterPosition(centerPos);

            if (zeroFaceIdxVec.size() == 4)
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);
                ElementPointer& element4 = interactionVolume.getSubVolumeElement(3);

                ElementPointer outside1 = element1;
                ElementPointer outside4 = element4;

                IntersectionIterator isIt1 = problem_.gridView().ibegin(*element1);
                IntersectionIterator isItEnd1 = problem_.gridView().iend(*element1);
                for (; isIt1 != isItEnd1; ++isIt1)
                {
                    if (isIt1->neighbor())
                    {
                        if (isIt1->indexInInside() == interactionVolume.getIndexOnElement(0, 2))
                        {
                            outside1 = isIt1->outside();
                            break;
                        }
                    }
                }
                IntersectionIterator isIt4 = problem_.gridView().ibegin(*element4);
                IntersectionIterator isItEnd4 = problem_.gridView().iend(*element4);
                for (; isIt4 != isItEnd4; ++isIt4)
                {
                    if (isIt4->neighbor())
                    {
                        if (isIt4->indexInInside() == interactionVolume.getIndexOnElement(3, 2))
                        {
                            outside4 = isIt4->outside();
                            break;
                        }
                    }
                }
                if (outside1 != outside4)
                {
                    interactionVolume.setHangingNodeType(InteractionVolume::fourSmallCellsEdge);
                }
                else if (outside1 == outside4)
                {
                    interactionVolume.setHangingNodeType(InteractionVolume::fourSmallCellsFace);
                }
            }

            if (interactionVolume.getHangingNodeType() == InteractionVolume::fourSmallCellsFace)
            {
                ElementPointer& element = interactionVolume.getSubVolumeElement(0);

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

                IntersectionIterator isIt = problem_.gridView().ibegin(*element);
                IntersectionIterator isEndIt = problem_.gridView().iend(*element);
                for (; isIt != isEndIt; ++isIt)
                {
                    if (isIt->indexInInside() == interactionVolume.getIndexOnElement(0, 2))
                    {
                        if (isIt->neighbor())
                        {
                            ElementPointer outside = isIt->outside();
                            if (element->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 4);
                                interactionVolume.setSubVolumeElement(outside, 5);
                                interactionVolume.setSubVolumeElement(outside, 6);
                                interactionVolume.setSubVolumeElement(outside, 7);

                                break;
                            }
                        }
                    }

                }
            }
            else if (interactionVolume.getHangingNodeType() == InteractionVolume::fourSmallCellsEdge)
            {
                ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);
                ElementPointer& element2 = interactionVolume.getSubVolumeElement(1);
                ElementPointer& element3 = interactionVolume.getSubVolumeElement(2);
                ElementPointer& element4 = interactionVolume.getSubVolumeElement(3);

                IntersectionIterator isIt1 = problem_.gridView().ibegin(*element1);
                IntersectionIterator isItEnd1 = problem_.gridView().iend(*element1);
                for (; isIt1 != isItEnd1; ++isIt1)
                {
                    if (isIt1->neighbor())
                    {
                        if (isIt1->indexInInside() == interactionVolume.getIndexOnElement(0, 2))
                        {
                            ElementPointer outside = isIt1->outside();
                            if (element1->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 4);
                            }
                        }
                    }
                }
                IntersectionIterator isIt2 = problem_.gridView().ibegin(*element2);
                IntersectionIterator isItEnd2 = problem_.gridView().iend(*element2);
                for (; isIt2 != isItEnd2; ++isIt2)
                {
                    if (isIt2->neighbor())
                    {
                        if (isIt2->indexInInside() == interactionVolume.getIndexOnElement(1, 2))
                        {
                            ElementPointer outside = isIt2->outside();
                            if (element2->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 5);

                                break;
                            }
                        }
                    }
                }
                IntersectionIterator isIt3 = problem_.gridView().ibegin(*element3);
                IntersectionIterator isItEnd3 = problem_.gridView().iend(*element3);
                for (; isIt3 != isItEnd3; ++isIt3)
                {
                    if (isIt3->neighbor())
                    {
                        if (isIt3->indexInInside() == interactionVolume.getIndexOnElement(2, 2))
                        {
                            ElementPointer outside = isIt3->outside();
                            if (element3->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 6);
                                break;
                            }
                        }
                    }
                }
                IntersectionIterator isIt4 = problem_.gridView().ibegin(*element4);
                IntersectionIterator isItEnd4 = problem_.gridView().iend(*element4);
                for (; isIt4 != isItEnd4; ++isIt4)
                {
                    if (isIt4->neighbor())
                    {
                        if (isIt4->indexInInside() == interactionVolume.getIndexOnElement(3, 2))
                        {
                            ElementPointer outside = isIt4->outside();
                            if (element4->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 7);

                                break;
                            }
                        }
                    }
                }

                DimVector edgeCoord1(interactionVolume.getEdgePosition(0));
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

                ElementPointer& element5 = interactionVolume.getSubVolumeElement(4);
                ElementPointer& element6 = interactionVolume.getSubVolumeElement(5);
                ElementPointer& element7 = interactionVolume.getSubVolumeElement(6);
                ElementPointer& element8 = interactionVolume.getSubVolumeElement(7);

                if (element5 == element6)
                {
                    interactionVolume.setFacePosition(element5->geometry().center(), 4);
                    interactionVolume.setFacePosition(element7->geometry().center(), 6);

                    IntersectionIterator isIt = problem_.gridView().ibegin(*element5);
                    IntersectionIterator isEndIt = problem_.gridView().iend(*element5);
                    for (; isIt != isEndIt; ++isIt)
                    {
                        if (isIt->neighbor())
                        {
                            ElementPointer outside = isIt->outside();

                            if (outside == element7 || outside == element8)
                            {
                                int indexInInside = isIt->indexInInside();
                                interactionVolume.setIndexOnElement(indexInInside, 4, 2);
                                interactionVolume.setIndexOnElement(indexInInside, 5, 1);
                                DimVector normal = isIt->centerUnitOuterNormal();
                                interactionVolume.setNormal(normal, 4, 2);
                                interactionVolume.setNormal(normal, 5, 1);
                                GlobalPosition globalPosFace(isIt->geometry().center());
                                interactionVolume.setFacePosition(globalPosFace, 5);
                                interactionVolume.setFacePosition(globalPosFace, 7);
                                int indexInOutside = isIt->indexInOutside();
                                interactionVolume.setIndexOnElement(indexInOutside, 6, 1);
                                interactionVolume.setIndexOnElement(indexInOutside, 7, 2);
                                normal *= -1;
                                interactionVolume.setNormal(normal, 6, 1);
                                interactionVolume.setNormal(normal, 7, 2);

                                break;
                            }
                        }
                    }

                    DimVector edgeCoord2(interactionVolume.getFacePosition(7));

                    globalPosFace = interactionVolume.getFacePosition(5);
                    crossProductVector1 = centerPos - globalPosFace;
                    crossProductVector2 = edgeCoord2 - edgeCoord4;
                    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
                    interactionVolume.setFaceArea(faceArea, 5);
                    //                        interactionVolume.setFaceArea(0.0, 5);
                    globalPosFace = interactionVolume.getFacePosition(7);
                    crossProductVector1 = centerPos - globalPosFace;
                    crossProductVector2 = edgeCoord2 - edgeCoord6;
                    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
                    interactionVolume.setFaceArea(faceArea, 7);
                    //                        interactionVolume.setFaceArea(0.0, 7);
                }
                else if (element5 == element7)
                {
                    interactionVolume.setFacePosition(element6->geometry().center(), 5);
                    interactionVolume.setFacePosition(element5->geometry().center(), 7);
                    interactionVolume.setFacePosition(globalPosFace, 6);

                    IntersectionIterator isIt = problem_.gridView().ibegin(*element5);
                    IntersectionIterator isEndIt = problem_.gridView().iend(*element5);
                    for (; isIt != isEndIt; ++isIt)
                    {
                        if (isIt->neighbor())
                        {
                            ElementPointer outside = isIt->outside();

                            if (outside == element6 || outside == element8)
                            {
                                int indexInInside = isIt->indexInInside();
                                interactionVolume.setIndexOnElement(indexInInside, 4, 1);
                                interactionVolume.setIndexOnElement(indexInInside, 6, 2);
                                DimVector normal = isIt->centerUnitOuterNormal();
                                interactionVolume.setNormal(normal, 4, 1);
                                interactionVolume.setNormal(normal, 6, 2);
                                GlobalPosition globalPosFace(isIt->geometry().center());
                                interactionVolume.setFacePosition(globalPosFace, 4);
                                interactionVolume.setFacePosition(globalPosFace, 6);
                                int indexInOutside = isIt->indexInOutside();
                                interactionVolume.setIndexOnElement(indexInOutside, 5, 2);
                                interactionVolume.setIndexOnElement(indexInOutside, 7, 1);
                                normal *= -1;
                                interactionVolume.setNormal(normal, 5, 2);
                                interactionVolume.setNormal(normal, 7, 1);

                                break;
                            }
                        }
                    }
                    DimVector edgeCoord2(interactionVolume.getFacePosition(4));

                    globalPosFace = interactionVolume.getFacePosition(4);
                    crossProductVector1 = centerPos - globalPosFace;
                    crossProductVector2 = edgeCoord2 - edgeCoord3;
                    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
                    interactionVolume.setFaceArea(faceArea, 4);

                    globalPosFace = interactionVolume.getFacePosition(6);
                    crossProductVector1 = centerPos - globalPosFace;
                    crossProductVector2 = edgeCoord2 - edgeCoord5;
                    faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
                    interactionVolume.setFaceArea(faceArea, 6);
                    //                        interactionVolume.setFaceArea(0.0, 4);
                    //                        interactionVolume.setFaceArea(0.0, 6);
                }

                const ElementGeometry& geometry = element5->geometry();

                const ReferenceElement& referenceElement = ReferenceElements::general(geometry.type());

                int oldSubVolumElemIdx = IndexTranslator::getOldElemIdxFromNewFaceIdxto0(zeroFaceIdx, 4);
                int oldEdgeIdx = IndexTranslator::getOldEdgeIdxFromNewFaceIdxto0(zeroFaceIdx, 1);

                DimVector edgeCoord(0.0);
                switch (oldSubVolumElemIdx)
                {
                case 0:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 2:
                            edgeCoord = geometry.global(referenceElement.position(9, dim - 1));
                            break;
                        case 0:
                            edgeCoord = geometry.global(referenceElement.position(3, dim - 1));
                            break;
                        case 5:
                            edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                            break;
                        }

                        break;
                    }
                case 1:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 0:
                            edgeCoord = geometry.global(referenceElement.position(2, dim - 1));
                            break;
                        case 2:
                            edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                            break;
                        case 3:
                            edgeCoord = geometry.global(referenceElement.position(11, dim - 1));
                            break;
                        }

                        break;
                    }
                case 2:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 0:
                            edgeCoord = geometry.global(referenceElement.position(1, dim - 1));
                            break;
                        case 4:
                            edgeCoord = geometry.global(referenceElement.position(9, dim - 1));
                            break;
                        case 5:
                            edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                            break;
                        }

                        break;
                    }
                case 3:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 0:
                            edgeCoord = geometry.global(referenceElement.position(0, dim - 1));
                            break;
                        case 4:
                            edgeCoord = geometry.global(referenceElement.position(8, dim - 1));
                            break;
                        case 3:
                            edgeCoord = geometry.global(referenceElement.position(10, dim - 1));
                            break;
                        }

                        break;
                    }
                case 4:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 1:
                            edgeCoord = geometry.global(referenceElement.position(3, dim - 1));
                            break;
                        case 2:
                            edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                            break;
                        case 5:
                            edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                            break;
                        }

                        break;
                    }
                case 5:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 1:
                            edgeCoord = geometry.global(referenceElement.position(2, dim - 1));
                            break;
                        case 2:
                            edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                            break;
                        case 3:
                            edgeCoord = geometry.global(referenceElement.position(7, dim - 1));
                            break;
                        }

                        break;
                    }
                case 6:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 1:
                            edgeCoord = geometry.global(referenceElement.position(1, dim - 1));
                            break;
                        case 4:
                            edgeCoord = geometry.global(referenceElement.position(5, dim - 1));
                            break;
                        case 5:
                            edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                            break;
                        }

                        break;
                    }
                case 7:
                    {
                        switch (oldEdgeIdx)
                        {
                        case 1:
                            edgeCoord = geometry.global(referenceElement.position(0, dim - 1));
                            break;
                        case 4:
                            edgeCoord = geometry.global(referenceElement.position(4, dim - 1));
                            break;
                        case 3:
                            edgeCoord = geometry.global(referenceElement.position(6, dim - 1));
                            break;
                        }

                        break;
                    }
                }

                interactionVolume.setEdgePosition(edgeCoord, 1);
            }
            else if (interactionVolume.getHangingNodeType() == InteractionVolume::fourSmallCellsDiag)
            {
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

                globalPosFace = interactionVolume.getFacePosition(3);
                crossProductVector1 = centerPos - globalPosFace;
                crossProductVector2 = edgeCoord1 - edgeCoord6;
                faceArea = crossProduct(crossProductVector1, crossProductVector2).two_norm()/2.0;
                interactionVolume.setFaceArea(faceArea, 3);

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

                ElementPointer& element1 = interactionVolume.getSubVolumeElement(0);

                bool hasFaceOne = false;
                bool hasFaceTwo = false;
                IntersectionIterator isIt = problem_.gridView().ibegin(*element1);
                IntersectionIterator isEndIt = problem_.gridView().iend(*element1);
                for (; isIt != isEndIt; ++isIt)
                {
                    if (isIt->indexInInside() == interactionVolume.getIndexOnElement(0, 1))
                    {
                        if (isIt->neighbor())
                        {
                            ElementPointer outside = isIt->outside();
                            if (element1->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 2);
                                interactionVolume.setSubVolumeElement(outside, 3);

                                hasFaceOne = true;
                                if (hasFaceTwo)
                                    break;
                            }
                        }
                    }
                    if (isIt->indexInInside() == interactionVolume.getIndexOnElement(0, 2))
                    {
                        if (isIt->neighbor())
                        {
                            ElementPointer outside = isIt->outside();
                            if (element1->level() > outside->level())
                            {
                                interactionVolume.setSubVolumeElement(outside, 4);
                                interactionVolume.setSubVolumeElement(outside, 5);

                                hasFaceTwo = true;
                                if (hasFaceOne)
                                    break;
                            }
                        }
                    }
                }
            }

            break;
        }
        //hanging-node interaction volume of type 2 or 6
    case 6:
        {
            InteractionVolume hangingNodeVolume;

            std::vector<int> elemIdxOld;
            for (int i = 0; i < InteractionVolume::subVolumeTotalNum; i++)
            {
                if (interactionVolume.hasSubVolumeElement(i))
                    elemIdxOld.push_back(i);
            }

            std::set<int> zeroFaceIdxVec;
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    int fIdx = IndexTranslator::getFaceIndexFromElements(elemIdxOld[i], elemIdxOld[j]);
                    if (fIdx >= 0)
                        zeroFaceIdxVec.insert(fIdx);
                }
            }

            std::vector<int> elemIdxNew(6);
            int zeroFaceIdx = 0;

            bool isFace = true;
            std::set<int>::iterator it = zeroFaceIdxVec.begin();
            for (; it != zeroFaceIdxVec.end(); ++it)
            {
                zeroFaceIdx = *it;
                isFace = true;
                //                    bool isFace = true;
                for (int i = 0; i < 6; i++)
                {
                    elemIdxNew[i] = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elemIdxOld[i]);
                    if (elemIdxNew[i] == 6 || elemIdxNew[i] == 7)
                    {
                        isFace = false;
                        break;
                    }
                }
                if (isFace)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        hangingNodeVolume.setSubVolumeElement(interactionVolume.getSubVolumeElement(elemIdxOld[i]),
                                                              elemIdxNew[i]);
                    }

                    break;
                }
            }

            for (int elem = 0; elem < InteractionVolume::subVolumeTotalNum; elem++)
            {
                for (int i = 0; i < 3; i++)
                {
                    int faceIdxOld = IndexTranslator::getFaceIndexFromSubVolume(elem, i);
                    int faceIdxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, faceIdxOld);

                    for (int j = 0; j < 3; j++)
                    {
                        int elemNew = IndexTranslator::getNewElemIdxFromOldFaceIdxto0(zeroFaceIdx, elem);
                        int faceIdxNewTest = IndexTranslator::getFaceIndexFromSubVolume(elemNew, j);

                        if (faceIdxNew == faceIdxNewTest)
                        {
                            hangingNodeVolume.setNormal(interactionVolume.getNormal(elem, i),
                                                        elemNew, j);
                            hangingNodeVolume.setIndexOnElement(
                                                                interactionVolume.getIndexOnElement(elem, i), elemNew, j);
                        }
                    }
                }
            }

            for (int i = 0; i < InteractionVolume::fluxFacesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewFaceIdxFromOldIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setFaceArea(interactionVolume.getFaceArea(i), idxNew);
                hangingNodeVolume.setFacePosition(interactionVolume.getFacePosition(i), idxNew);
            }

            for (int i = 0; i < InteractionVolume::fluxEdgesTotalNum; i++)
            {
                int idxNew = IndexTranslator::getNewEdgeIdxFromOldFaceIdxto0(zeroFaceIdx, i);
                hangingNodeVolume.setEdgePosition(interactionVolume.getEdgePosition(i), idxNew);
            }

            interactionVolume = hangingNodeVolume;

            interactionVolume.setCenterPosition(centerPos);

            interactionVolume.setHangingNodeType(InteractionVolume::sixSmallCells);

            ElementPointer& element3 = interactionVolume.getSubVolumeElement(2);

            IntersectionIterator isIt = problem_.gridView().ibegin(*element3);
            IntersectionIterator isEndIt = problem_.gridView().iend(*element3);
            for (; isIt != isEndIt; ++isIt)
            {
                if (isIt->indexInInside() == interactionVolume.getIndexOnElement(2, 2))
                {
                    if (isIt->neighbor())
                    {
                        ElementPointer outside = isIt->outside();
                        if (element3->level() > outside->level())
                        {
                            interactionVolume.setSubVolumeElement(outside, 6);
                            interactionVolume.setSubVolumeElement(outside, 7);

                            break;
                        }
                    }
                }
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
            break;
        }
    default:
        {
            interactionVolume.printInteractionVolumeInfo();
            DUNE_THROW(Dune::NotImplemented, "Hanging node shape not implemented");
            break;
        }
    }
}

//! \brief Stores interaction volumes for each grid vertex
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainerAdaptive<TypeTag>::storeInteractionVolumeInfo()
{
    std::vector < std::vector<int> > elemVertMap(problem_.gridView().size(dim), std::vector<int>(8, -1));

    //Add elements to the interaction volumes and store element-vertex map
    ElementIterator eEndIt = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt; ++eIt)
    {
        ElementPointer ePtr = *eIt;
        asImp_().storeSubVolumeElements(*ePtr, elemVertMap);
    }

    for (unsigned int i = 0; i < asImp_().interactionVolumes_.size(); i++)
        if (asImp_().interactionVolumes_[i].getElementNumber() == 0)
            asImp_().interactionVolumes_[i].printInteractionVolumeInfo();

    // Store information related to DUNE intersections for all interaction volumes
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt; ++eIt)
    {
        ElementPointer ePtr = *eIt;
        asImp_().storeIntersectionInfo(*ePtr, elemVertMap);
    }

    faceVertices_.clear();
    faceVertices_.resize(problem_.gridView().size(0));

    // Complete storage of the interaction volumes using the previously stored information
    // about the orientation and relationship of the DUNE elements in the interaction volumes (see doc/docextra/3dmpfa)
    VertexIterator vEndIt = problem_.gridView().template end<dim>();
    for (VertexIterator vIt = problem_.gridView().template begin<dim>(); vIt != vEndIt; ++vIt)
    {
        int vIdxGlobal = problem_.variables().index(*vIt);

        InteractionVolume& interactionVolume = asImp_().interactionVolumes_[vIdxGlobal];

        if (interactionVolume.getElementNumber() == 8)
        {
            asImp_().storeInnerInteractionVolume(interactionVolume, *vIt);
        }
        else if (interactionVolume.isBoundaryInteractionVolume())
        {
            asImp_().storeBoundaryInteractionVolume(interactionVolume, *vIt);
        }
        //hanging node!
        else
        {
            storeHangingNodeInteractionVolume(interactionVolume, *vIt);
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

            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 0), globalIdx1, interactionVolume.getIndexOnElement(0, 0));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 1), globalIdx1, interactionVolume.getIndexOnElement(0, 1));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 2), globalIdx1, interactionVolume.getIndexOnElement(0, 2));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 0), globalIdx2, interactionVolume.getIndexOnElement(1, 0));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 1), globalIdx2, interactionVolume.getIndexOnElement(1, 1));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 2), globalIdx2, interactionVolume.getIndexOnElement(1, 2));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 0), globalIdx3, interactionVolume.getIndexOnElement(2, 0));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 1), globalIdx4, interactionVolume.getIndexOnElement(3, 1));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 0), globalIdx5, interactionVolume.getIndexOnElement(4, 0));
            asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 0), globalIdx6, interactionVolume.getIndexOnElement(5, 0));

            faceVertices_[globalIdx1][interactionVolume.getIndexOnElement(0, 0)].insert(vIdxGlobal);
            faceVertices_[globalIdx1][interactionVolume.getIndexOnElement(0, 1)].insert(vIdxGlobal);
            faceVertices_[globalIdx1][interactionVolume.getIndexOnElement(0, 2)].insert(vIdxGlobal);
            faceVertices_[globalIdx2][interactionVolume.getIndexOnElement(1, 0)].insert(vIdxGlobal);
            faceVertices_[globalIdx2][interactionVolume.getIndexOnElement(1, 1)].insert(vIdxGlobal);
            faceVertices_[globalIdx2][interactionVolume.getIndexOnElement(1, 2)].insert(vIdxGlobal);
            faceVertices_[globalIdx3][interactionVolume.getIndexOnElement(2, 0)].insert(vIdxGlobal);
            faceVertices_[globalIdx4][interactionVolume.getIndexOnElement(3, 1)].insert(vIdxGlobal);
            faceVertices_[globalIdx5][interactionVolume.getIndexOnElement(4, 0)].insert(vIdxGlobal);
            faceVertices_[globalIdx6][interactionVolume.getIndexOnElement(5, 0)].insert(vIdxGlobal);

            if (interactionVolume.getHangingNodeType() != InteractionVolume::twoSmallCells)
            {
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 1), globalIdx3, interactionVolume.getIndexOnElement(2, 1));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 2), globalIdx3, interactionVolume.getIndexOnElement(2, 2));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 0), globalIdx4, interactionVolume.getIndexOnElement(3, 0));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 2), globalIdx4, interactionVolume.getIndexOnElement(3, 2));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 1), globalIdx5, interactionVolume.getIndexOnElement(4, 1));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 2), globalIdx5, interactionVolume.getIndexOnElement(4, 2));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 1), globalIdx6, interactionVolume.getIndexOnElement(5, 1));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 2), globalIdx6, interactionVolume.getIndexOnElement(5, 2));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 0), globalIdx7, interactionVolume.getIndexOnElement(6, 0));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 1), globalIdx7, interactionVolume.getIndexOnElement(6, 1));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 2), globalIdx7, interactionVolume.getIndexOnElement(6, 2));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 0), globalIdx8, interactionVolume.getIndexOnElement(7, 0));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 1), globalIdx8, interactionVolume.getIndexOnElement(7, 1));
                asImp_().addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 2), globalIdx8, interactionVolume.getIndexOnElement(7, 2));

                faceVertices_[globalIdx3][interactionVolume.getIndexOnElement(2, 1)].insert(vIdxGlobal);
                faceVertices_[globalIdx3][interactionVolume.getIndexOnElement(2, 2)].insert(vIdxGlobal);
                faceVertices_[globalIdx4][interactionVolume.getIndexOnElement(3, 0)].insert(vIdxGlobal);
                faceVertices_[globalIdx4][interactionVolume.getIndexOnElement(3, 2)].insert(vIdxGlobal);
                faceVertices_[globalIdx5][interactionVolume.getIndexOnElement(4, 1)].insert(vIdxGlobal);
                faceVertices_[globalIdx5][interactionVolume.getIndexOnElement(4, 2)].insert(vIdxGlobal);
                faceVertices_[globalIdx6][interactionVolume.getIndexOnElement(5, 1)].insert(vIdxGlobal);
                faceVertices_[globalIdx6][interactionVolume.getIndexOnElement(5, 2)].insert(vIdxGlobal);
                faceVertices_[globalIdx7][interactionVolume.getIndexOnElement(6, 0)].insert(vIdxGlobal);
                faceVertices_[globalIdx7][interactionVolume.getIndexOnElement(6, 1)].insert(vIdxGlobal);
                faceVertices_[globalIdx7][interactionVolume.getIndexOnElement(6, 2)].insert(vIdxGlobal);
                faceVertices_[globalIdx8][interactionVolume.getIndexOnElement(7, 0)].insert(vIdxGlobal);
                faceVertices_[globalIdx8][interactionVolume.getIndexOnElement(7, 1)].insert(vIdxGlobal);
                faceVertices_[globalIdx8][interactionVolume.getIndexOnElement(7, 2)].insert(vIdxGlobal);
            }
        }

    }
}
} // end of Dune namespace
#endif
