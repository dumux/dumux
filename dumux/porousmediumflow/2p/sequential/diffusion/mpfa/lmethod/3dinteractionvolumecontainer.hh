// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPModel
 * \brief Interactionvolume container for 3-d MPFA L-method.
 */
#ifndef DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_HH
#define DUMUX_FVMPFAL3D_INTERACTIONVOLUMECONTAINER_HH

// dumux environment
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/linteractionvolume3d.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Interactionvolume container for 3-d MPFA L-method.
 *
 * Container class which stores MPFA-interaction-volume information for each vertex of a DUNE grid.
 * Each <tt>InteractionVolume</tt> object stores the information which is necessary to calculate MPFA transmissibility matrices:
 *
 * - relationship and orientation of the elements around a vertex (see doc/docextra/3dmpfa)
 * - geometric information, such as element/face/edge positions, normals, ...
 */
template<class TypeTag>
class FvMpfaL3dInteractionVolumeContainer
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Implementation = GetPropType<TypeTag, Properties::MPFAInteractionVolumeContainer>;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;
    using ElementGeometry = typename Element::Geometry;

    using Intersection = typename GridView::Intersection;
    using IntersectionGeometry = typename Intersection::Geometry;

    using GlobalPosition = typename ElementGeometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

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
    //! Type for storing an MPFA-interaction-volume. (Usually of type FvMpfaL3dInteractionVolume or FvMpfaL3dInteractionVolumeAdaptive)
    using InteractionVolume = GetPropType<TypeTag, Properties::MPFAInteractionVolume>;

private:
    using GlobalInteractionVolumeVector = std::vector<InteractionVolume>;
    using FaceAreaVector = std::vector<Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2*dim> >;
protected:
    void storeSubVolumeElements(const Element& element, std::vector < std::vector<int> >& elemVertMap);
    void storeIntersectionInfo(const Element& element, std::vector < std::vector<int> >& elemVertMap);
    void storeInnerInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex, bool sameLevel = true);
    void storeBoundaryInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex);
private:
    void storeInteractionVolumeInfo();
public:

    /*!
     * \brief Updates the interaction volume container
     *
     * Rebuilds and stores the interaction volumes for the entire grid
     */
    void update()
    {
        interactionVolumes_.clear();
        realFluxFaceArea_.clear();

        realFluxFaceArea_.resize(problem_.gridView().size(dim),
                                 Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2 * dim>(Dune::FieldVector<Scalar, 2>(0.0)));
        interactionVolumes_.resize(problem_.gridView().size(dim), InteractionVolume(problem_.gridView().grid()));

        asImp_().storeInteractionVolumeInfo();
    }


    /*!
     * \brief Initializes the interaction volume container
     *
     * Builds and stores the interaction volumes for the entire grid
     */
    void initialize(bool solveTwice = true)
    {
        update();

        return;
    }

    /*!
     * \brief Returns an interaction volume
     *
     * \param vertexIdx Global index of a vertex in the DUNE grid
     */
    InteractionVolume& interactionVolume(int vertexIdx)
    {
        return interactionVolumes_[vertexIdx];
    }

    /*!
     * \brief Returns an interaction volume
     *
     * \param vertexIdx Global index of a vertex in the DUNE grid
     */
    InteractionVolume& interactionVolume(int vertexIdx) const
    {
        return interactionVolumes_[vertexIdx];
    }

    //! Returns the interaction volumes container
    GlobalInteractionVolumeVector& interactionVolumesGlobal()
    {
        return interactionVolumes_;
    }

    //! Returns the interaction volumes container
    GlobalInteractionVolumeVector& interactionVolumesGlobal() const
    {
        return interactionVolumes_;
    }

    /*!
     * \brief Returns the area weighting factor for the fluxes
     *
     * \param interactionVolume An interaction volume object
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param elemLocalIdx Local index of an element in the interaction volume
     * \param localFaceIdx  Local index of a flux face with respect to an element of the interaction volume
     *
     * \return Ratio of the element face area and the flux face area through which fluxes are calculated by the MPFA method
     *  (1 if an element does not touches the domain boundary!)
     */
    Scalar faceAreaFactor(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = getRealFaceArea(interactionVolume, elemGlobalIdx, elemLocalIdx, localFaceIdx);
        factor /= getRealFluxFaceArea(interactionVolume, elemGlobalIdx, elemLocalIdx, localFaceIdx);

        return factor;
    }

    /*!
     * \brief Returns the area weighting factor for the fluxes
     *
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param indexInInside Local index of the face in the DUNE reference element
     *
     * \return Ratio of the element face area and the flux face area through which fluxes are calculated by the MPFA method
     *  (1 if an element does not touches the domain boundary!)
     */
    Scalar faceAreaFactor(int elemGlobalIdx, int indexInInside)
    {
        Scalar factor = getRealFaceArea(elemGlobalIdx, indexInInside);
        factor /= getRealFluxFaceArea(elemGlobalIdx, indexInInside);

        return factor;
    }

    /*!
     * \brief Returns the area trough which fluxes are calculated by the MPFA
     *
     * \param interactionVolume An interaction volume object
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param elemLocalIdx Local index of an element in the interaction volume
     * \param localFaceIdx  Local index of a flux face with respect to an element of the interaction volume
     *
     * \return flux face area (equal to the element face area if an element does not touches the domain boundary!)
     */
    Scalar getRealFluxFaceArea(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][interactionVolume.getIndexOnElement(elemLocalIdx, localFaceIdx)][fluxFaceArea];

        return factor;
    }

    /*!
     * \brief  Returns the area trough which fluxes are calculated by the MPFA
     *
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param indexInInside Local index of the face in the DUNE reference element
     *
     * \return flux face area (equal to the element face area if an element does not touches the domain boundary!)
     */
    Scalar getRealFluxFaceArea(int elemGlobalIdx, int indexInInside)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][indexInInside][fluxFaceArea];

        return factor;
    }

    /*!
     * \brief Returns the face area of the element
     *
     * \param interactionVolume An interaction volume object
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param elemLocalIdx Local index of an element in the interaction volume
     * \param localFaceIdx  Local index of a flux face with respect to an element of the interaction volume
     *
     * \return the face area of the element
     */
    Scalar getRealFaceArea(InteractionVolume& interactionVolume, int elemGlobalIdx, int elemLocalIdx, int localFaceIdx)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][interactionVolume.getIndexOnElement(elemLocalIdx, localFaceIdx)][realFaceArea];

        return factor;
    }

    /*!
     * \brief Returns the face area of the element
     *
     * \param elemGlobalIdx Global index of an element in the DUNE grid
     * \param indexInInside Local index of the face in the DUNE reference element
     *
     * \return the face area of the element
     */
    Scalar getRealFaceArea(int elemGlobalIdx, int indexInInside)
    {
        Scalar factor = realFluxFaceArea_[elemGlobalIdx][indexInInside][realFaceArea];

        return factor;
    }

    /*!
     * \brief Constructs a FvMpfaL3dInteractionVolumeContainer object
     *
     * \param problem A problem class object
     */
    FvMpfaL3dInteractionVolumeContainer(Problem& problem) :
        problem_(problem)
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }
    }

protected:
    void addRealFluxFaceArea_(Scalar faceArea, int eIdxGlobal, int fIdx)
    {
        realFluxFaceArea_[eIdxGlobal][fIdx][fluxFaceArea] += faceArea;
    }
    void addRealFaceArea_(Scalar faceArea, int eIdxGlobal, int fIdx)
    {
        realFluxFaceArea_[eIdxGlobal][fIdx][realFaceArea] += faceArea;
    }

    Problem& problem_;

    GlobalInteractionVolumeVector interactionVolumes_;
    FaceAreaVector realFluxFaceArea_;
private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

/*!
 * \brief Function for storing the elements of an interaction volume and
 * constructing a map from a vertex to its surrounding elements
 *
 * Stores an element in all interaction volumes it belongs to. Additionally,
 * the global index of an element is stored at the position of its local index according to a DUNE reference element
 * for every vertex of the element.
 *
 * \param element A level 0 Entity of a DUNE grid
 * \param elemVertMap Vector containing the global vertex-element map
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements(const Element& element,
                                                                          std::vector < std::vector<int> >& elemVertMap)
{
    int eIdxGlobal = problem_.variables().index(element);
    int vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 0, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 7);
    elemVertMap[vIdxGlobal][7] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 1, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 6);
    elemVertMap[vIdxGlobal][6] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 2, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 5);
    elemVertMap[vIdxGlobal][5] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 3, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 4);
    elemVertMap[vIdxGlobal][4] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 4, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 3);
    elemVertMap[vIdxGlobal][3] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 5, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 2);
    elemVertMap[vIdxGlobal][2] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 6, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 1);
    elemVertMap[vIdxGlobal][1] = eIdxGlobal;

    vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, 7, dim);
    interactionVolumes_[vIdxGlobal].setSubVolumeElement(element, 0);
    elemVertMap[vIdxGlobal][0] = eIdxGlobal;

}

/*!
 * \brief Stores information with respect to DUNE intersections in the interaction volumes
 *
 * Stores information with respect to DUNE intersections, such as normals,
 * in the interaction volumes. Assumes a local storage following the DUNE
 * reference element index, which is performed by the function
 * FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements(const Element& element,
 * std::vector < std::vector<int> >& elemVertMap).
 *
 * \param element A level 0 Entity of a DUNE grid
 * \param elemVertMap Vector containing the global vertex-element map
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeIntersectionInfo(const Element& element,
                                                                         std::vector < std::vector<int> >& elemVertMap)
{
    BoundaryTypes bcType;

    int eIdxGlobal = problem_.variables().index(element);

    const ElementGeometry& geometry = element.geometry();
    const auto refElement = referenceElement(geometry);

    int levelI = element.level();

    // run through all intersections
    for (const auto& intersection : intersections(problem_.gridView(), element))
    {
        int indexInInside = intersection.indexInInside();

        DimVector normal = intersection.centerUnitOuterNormal();

        const IntersectionGeometry& isGeometry = intersection.geometry();

        Scalar faceVol = isGeometry.volume();

        const DimVector& globalPosFace = isGeometry.center();

        bool takeIntersection = true;
        if (intersection.neighbor())
        {
            auto outside = intersection.outside();
            int eIdxGlobalJ = problem_.variables().index(outside);

            if (levelI == outside.level() && eIdxGlobal > eIdxGlobalJ)
                takeIntersection = false;
            if (levelI < outside.level())
                takeIntersection = false;
        }

        if (takeIntersection)
        {
            addRealFaceArea_(faceVol, eIdxGlobal, indexInInside);
            if (intersection.neighbor())
            {
                int eIdxGlobalJ = problem_.variables().index(intersection.outside());
                addRealFaceArea_(faceVol, eIdxGlobalJ, intersection.indexInOutside());
            }

            for (int i = 0; i < isGeometry.corners(); i++)
            {
                int localVertIdx = refElement.subEntity(indexInInside, 1, i, dim);

                int vIdxGlobal = problem_.variables().vertexMapper().subIndex(element, localVertIdx, dim);

                InteractionVolume& interactionVolume = interactionVolumes_[vIdxGlobal];

                if (elemVertMap[vIdxGlobal][0] == eIdxGlobal)
                {
                    if (indexInInside == 1)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 0, 0);
                        interactionVolume.setNormal(normal, 0, 0);
                        interactionVolume.setFacePosition(globalPosFace, 0);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 4, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 4, 0);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][1] == eIdxGlobal)
                {
                    if (indexInInside == 3)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 1, 0);
                        interactionVolume.setNormal(normal, 1, 0);
                        interactionVolume.setFacePosition(globalPosFace, 1);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 5, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 5, 0);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][2] == eIdxGlobal)
                {
                    if (indexInInside == 2)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 2, 0);
                        interactionVolume.setNormal(normal, 2, 0);
                        interactionVolume.setFacePosition(globalPosFace, 3);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 6, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 6, 0);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][3] == eIdxGlobal)
                {
                    if (indexInInside == 0)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 3, 0);
                        interactionVolume.setNormal(normal, 3, 0);
                        interactionVolume.setFacePosition(globalPosFace, 2);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 7, 0);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 7, 0);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][4] == eIdxGlobal)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 4, 0);
                        interactionVolume.setNormal(normal, 4, 0);
                        interactionVolume.setFacePosition(globalPosFace, 8);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 6, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 6, 1);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][5] == eIdxGlobal)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 5, 0);
                        interactionVolume.setNormal(normal, 5, 0);
                        interactionVolume.setFacePosition(globalPosFace, 9);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 4, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 4, 1);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][6] == eIdxGlobal)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 6, 0);
                        interactionVolume.setNormal(normal, 6, 0);
                        interactionVolume.setFacePosition(globalPosFace, 11);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 7, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 7, 1);
                        }
                    }
                }
                if (elemVertMap[vIdxGlobal][7] == eIdxGlobal)
                {
                    if (indexInInside == 4)
                    {
                        interactionVolume.setIndexOnElement(indexInInside, 7, 0);
                        interactionVolume.setNormal(normal, 7, 0);
                        interactionVolume.setFacePosition(globalPosFace, 10);

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

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

                        if (intersection.neighbor())
                        {
                            int indexInOutside = intersection.indexInOutside();

                            interactionVolume.setIndexOnElement(indexInOutside, 5, 1);
                            DimVector normalOutside = normal;
                            normalOutside *= -1;

                            interactionVolume.setNormal(normalOutside, 5, 1);
                        }
                    }
                }
                if (intersection.boundary())
                {
                    if (elemVertMap[vIdxGlobal][0] == eIdxGlobal)
                    {
                        if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 0);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 0);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 0);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 3);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 3);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 3);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 8);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 8);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 8);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][1] == eIdxGlobal)
                    {
                        if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 1);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 1);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 1);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 0);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 0);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 0);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 9);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 9);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 9);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][2] == eIdxGlobal)
                    {
                        if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 3);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 3);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 3);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 2);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 2);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 2);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 11);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 11);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 11);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][3] == eIdxGlobal)
                    {
                        if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 2);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 2);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 2);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 1);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 1);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 1);
                            }
                        }
                        else if (indexInInside == 5)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 10);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 10);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 10);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][4] == eIdxGlobal)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 8);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 8);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 8);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 4);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 4);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 4);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 7);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 7);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 7);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][5] == eIdxGlobal)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 9);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 9);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 9);
                            }
                        }
                        else if (indexInInside == 3)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 5);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 5);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 5);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 4);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 4);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 4);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][6] == eIdxGlobal)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 11);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 11);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 11);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 7);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 7);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 7);
                            }
                        }
                        else if (indexInInside == 1)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 6);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 6);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 6);
                            }
                        }
                    }
                    if (elemVertMap[vIdxGlobal][7] == eIdxGlobal)
                    {
                        if (indexInInside == 4)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 10);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 10);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 10);
                            }
                        }
                        else if (indexInInside == 0)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 6);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 6);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 6);
                            }
                        }
                        else if (indexInInside == 2)
                        {
                            problem_.boundaryTypes(bcType, intersection);
                            PrimaryVariables boundValues(0.0);

                            interactionVolume.setBoundary(bcType, 5);
                            if (bcType.isNeumann(pressureEqIdx))
                            {
                                problem_.neumann(boundValues, intersection);
                                boundValues *= faceVol/4.0;
                                interactionVolume.setNeumannCondition(boundValues, 5);
                            }
                            if (bcType.hasDirichlet())
                            {
                                problem_.dirichlet(boundValues, intersection);
                                interactionVolume.setDirichletCondition(boundValues, 5);
                            }
                        }
                    }
                }
            }
        }

    }
}

/*!
 * \brief Stores additional information which can be constructed for interaction volumes of non-boundary vertices.
 *
 * Stores additional information which can be constructed for interaction volumes of non-boundary vertices:
 *
 *  - edge coordinates (coordinates of edge-continuity-points)
 *  - flux face areas
 *
 *  Assumes a local storage following the DUNE reference element index, which is performed by the
 *  function FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements(const Element& element,
 *                                                                                       std::vector < std::vector<int> >& elemVertMap).
 *
 * \param interactionVolume An interaction volume object
 * \param vertex The vertex (level dim entity) for which the interaction volume is stored
 * \param sameLevel Level indicator: true if all elements of an interaction volume are of the same level
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeInnerInteractionVolume(InteractionVolume& interactionVolume,
                                                                               const Vertex& vertex, bool sameLevel)
{
    const DimVector& centerPos = vertex.geometry().center();

    interactionVolume.setCenterPosition(centerPos);

    if (sameLevel)
    {
        auto element1 = interactionVolume.getSubVolumeElement(0);
        auto element8 = interactionVolume.getSubVolumeElement(7);

        const ElementGeometry& geometry1 = element1.geometry();
        const ElementGeometry& geometry8 = element8.geometry();

        const auto refElement = referenceElement(geometry1);

        DimVector edgeCoord(geometry1.global(refElement.position(9, dim - 1)));
        interactionVolume.setEdgePosition(edgeCoord, 2);
        edgeCoord = geometry1.global(refElement.position(3, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 0);
        edgeCoord = geometry1.global(refElement.position(11, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 5);

        edgeCoord = geometry8.global(refElement.position(4, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 4);
        edgeCoord = geometry8.global(refElement.position(6, dim - 1));
        interactionVolume.setEdgePosition(edgeCoord, 3);
        edgeCoord = geometry8.global(refElement.position(0, dim - 1));
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

/*!
 * \brief Stores additional information for interaction volumes of boundary vertices.
 *
 * Stores additional information for interaction volumes of boundary vertices:
 *
 *  - boundary conditions
 *  - information for flux weighting along boundary faces (see  Wolff 2013: http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/, or
 *  M. Wolff, Y. Cao, B. Flemisch, R. Helmig, and B. Wohlmuth (2013a). Multi-point flux
 * approximation L-method in 3D: numerical convergence and application to two-phase
 * flow through porous media. In P. Bastian, J. Kraus, R. Scheichl, and M. Wheeler,
 * editors, Simulation of Flow in Porous Media - Applications in Energy and Environment. De Gruyter.)
 *
 * Assumes a local storage following the DUNE reference element index, which is performed by the
 * function FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeSubVolumeElements(const Element& element,
 *                                                                                      std::vector < std::vector<int> >& elemVertMap).
 *
 * \param interactionVolume An interaction volume object
 * \param vertex The vertex (level dim entity) for which the interaction volume is stored
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeBoundaryInteractionVolume(InteractionVolume& interactionVolume,
                                                                                  const Vertex& vertex)
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(0));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 0, 0)/4.0,0);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 0, 1)/4.0,3);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 0, 2)/4.0,8);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 1, 0)/4.0,1);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 1, 1)/4.0,0);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 1, 2)/4.0,9);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 2, 0)/4.0,3);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 2, 1)/4.0,2);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 2, 2)/4.0,11);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 3, 0)/4.0,2);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 3, 1)/4.0,1);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 3, 2)/4.0,10);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 4, 0)/4.0,8);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 4, 1)/4.0,4);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 4, 2)/4.0,7);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(5));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 5, 0)/4.0,9);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 5, 1)/4.0,5);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 5, 2)/4.0,4);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 6, 0)/4.0,11);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 6, 1)/4.0,7);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 6, 2)/4.0,6);
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

                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(7));
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 7, 0)/4.0,10);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 7, 1)/4.0,6);
                interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal, 7, 2)/4.0,5);
            }
            break;
        }
    case 2:
        {
            // edge
            if (interactionVolume.hasSubVolumeElement(0))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(0));
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(11);

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 1)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 2)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 0)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 2)/4.0,9);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 0)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 2)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 1)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 2)/4.0,11);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 0)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 1)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 1)/4.0,4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 2)/4.0,7);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(7))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(7));
                if (interactionVolume.hasSubVolumeElement(5))
                {
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(5));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 0)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 1)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 5, 0)/4.0,9);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 5, 2)/4.0,4);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 0)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 2)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 6, 0)/4.0,11);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 6, 1)/4.0,7);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 1)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 2)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 3, 0)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 3, 1)/4.0,1);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(5))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(5));
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(10);

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 5, 1)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 5, 2)/4.0,4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 0)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 1)/4.0,0);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 5, 0)/4.0,9);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 5, 1)/4.0,5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 0)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 2)/4.0,7);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(6))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                if (interactionVolume.hasSubVolumeElement(4))
                {
                    interactionVolume.setOutsideFace(1);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(9);
                    interactionVolume.setOutsideFace(10);
                    interactionVolume.setOutsideFace(0);
                    interactionVolume.setOutsideFace(2);
                    interactionVolume.setOutsideFace(3);

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 6, 0)/4.0,11);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 6, 2)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 0)/4.0,8);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 4, 1)/4.0,4);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 6, 1)/4.0,7);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 6, 2)/4.0,6);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 0)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 1)/4.0,2);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(3))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    interactionVolume.setOutsideFace(4);
                    interactionVolume.setOutsideFace(5);
                    interactionVolume.setOutsideFace(6);
                    interactionVolume.setOutsideFace(7);
                    interactionVolume.setOutsideFace(3);
                    interactionVolume.setOutsideFace(8);
                    interactionVolume.setOutsideFace(11);

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 3, 0)/4.0,2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 3, 2)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 1)/4.0,0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 2)/4.0,9);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 3, 1)/4.0,1);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 3, 2)/4.0,10);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 0)/4.0,3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 2)/4.0,11);

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
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(0));
                if (interactionVolume.hasSubVolumeElement(1))
                {
                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                    if (interactionVolume.hasSubVolumeElement(2) && interactionVolume.hasSubVolumeElement(3))
                    {
                        interactionVolume.setOutsideFace(4);
                        interactionVolume.setOutsideFace(5);
                        interactionVolume.setOutsideFace(6);
                        interactionVolume.setOutsideFace(7);

                        int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                        int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 2)/4.0, 8);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 2)/4.0, 9);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 2, 2)/4.0, 10);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 3, 2)/4.0, 11);

                        return;
                    }
                    if (interactionVolume.hasSubVolumeElement(4) && interactionVolume.hasSubVolumeElement(5))
                    {
                        interactionVolume.setOutsideFace(2);
                        interactionVolume.setOutsideFace(6);
                        interactionVolume.setOutsideFace(10);
                        interactionVolume.setOutsideFace(11);

                        int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                        int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(5));
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 1)/4.0, 3);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 1, 0)/4.0, 1);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 4, 2)/4.0, 7);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 5, 1)/4.0, 5);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                    int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                    int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 0, 0)/4.0, 0);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 2, 1)/4.0, 2);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 4, 1)/4.0, 4);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 6, 2)/4.0, 6);

                    return;
                }
            }
            if (interactionVolume.hasSubVolumeElement(7))
            {
                int eIdxGlobal1 = problem_.variables().index(interactionVolume.getSubVolumeElement(7));
                if (interactionVolume.hasSubVolumeElement(5))
                {
                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(5));
                    if (interactionVolume.hasSubVolumeElement(1) && interactionVolume.hasSubVolumeElement(3))
                    {
                        interactionVolume.setOutsideFace(3);
                        interactionVolume.setOutsideFace(7);
                        interactionVolume.setOutsideFace(8);
                        interactionVolume.setOutsideFace(11);

                        int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));
                        int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 1)/4.0, 6);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 5, 2)/4.0, 4);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 1, 1)/4.0, 0);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 3, 0)/4.0, 2);

                        return;
                    }
                    if (interactionVolume.hasSubVolumeElement(4) && interactionVolume.hasSubVolumeElement(6))
                    {
                        interactionVolume.setOutsideFace(0);
                        interactionVolume.setOutsideFace(1);
                        interactionVolume.setOutsideFace(2);
                        interactionVolume.setOutsideFace(3);

                        int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(4));
                        int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 0)/4.0, 10);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 5, 0)/4.0, 9);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 4, 0)/4.0, 8);
                        interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 6, 0)/4.0, 11);

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

                    int eIdxGlobal2 = problem_.variables().index(interactionVolume.getSubVolumeElement(6));
                    int eIdxGlobal3 = problem_.variables().index(interactionVolume.getSubVolumeElement(2));
                    int eIdxGlobal4 = problem_.variables().index(interactionVolume.getSubVolumeElement(3));
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal1, 7, 2)/4.0, 5);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal2, 6, 1)/4.0, 7);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal3, 2, 0)/4.0, 3);
                    interactionVolume.setFaceArea(getRealFaceArea(interactionVolume, eIdxGlobal4, 3, 1)/4.0, 1);

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


/*!
 * \brief Stores interaction volumes for each grid vertex
 */
template<class TypeTag>
void FvMpfaL3dInteractionVolumeContainer<TypeTag>::storeInteractionVolumeInfo()
{
    std::vector < std::vector<int> > elemVertMap(problem_.gridView().size(dim), std::vector<int>(8, -1));

    //Add elements to the interaction volumes and store element-vertex map
    for (const auto& element : elements(problem_.gridView()))
        storeSubVolumeElements(element, elemVertMap);

    for (unsigned int i = 0; i < interactionVolumes_.size(); i++)
        if (interactionVolumes_[i].getElementNumber() == 0)
            interactionVolumes_[i].printInteractionVolumeInfo();

    // Store information related to DUNE intersections for all interaction volumes
    for (const auto& element : elements(problem_.gridView()))
        storeIntersectionInfo(element, elemVertMap);

    // Complete storage of the interaction volumes using the previously stored information
    // about the orientation and relationship of the DUNE elements in the interaction volumes (see doc/docextra/3dmpfa)
    for (const auto& vertex : vertices(problem_.gridView()))
    {
        int vIdxGlobal = problem_.variables().index(vertex);

        InteractionVolume& interactionVolume = interactionVolumes_[vIdxGlobal];

        if (interactionVolume.getElementNumber() == 8)
        {
            storeInnerInteractionVolume(interactionVolume, vertex);
        }
        else if (interactionVolume.isBoundaryInteractionVolume())
        {
            storeBoundaryInteractionVolume(interactionVolume, vertex);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,"Interaction volume is no boundary volume but consists of less than 8 elements");
        }

        // Store information about the MPFA flux face areas for correct flux weighting of
        // fluxes though faces which intersect the domain boundary
        // (see  M. Wolff, Y. Cao, B. Flemisch, R. Helmig, and B. Wohlmuth (2013a). Multi-point flux
        // approximation L-method in 3D: numerical convergence and application to two-phase
        // flow through porous media. In P. Bastian, J. Kraus, R. Scheichl, and M. Wheeler,
        // editors, Simulation of Flow in Porous Media - Applications in Energy and Environment. De Gruyter.)
        if (!interactionVolume.isBoundaryInteractionVolume())
        {
            auto element1 = interactionVolume.getSubVolumeElement(0);
            auto element2 = interactionVolume.getSubVolumeElement(1);
            auto element3 = interactionVolume.getSubVolumeElement(2);
            auto element4 = interactionVolume.getSubVolumeElement(3);
            auto element5 = interactionVolume.getSubVolumeElement(4);
            auto element6 = interactionVolume.getSubVolumeElement(5);
            auto element7 = interactionVolume.getSubVolumeElement(6);
            auto element8 = interactionVolume.getSubVolumeElement(7);

            int eIdxGlobal1 = problem_.variables().index(element1);
            int eIdxGlobal2 = problem_.variables().index(element2);
            int eIdxGlobal3 = problem_.variables().index(element3);
            int eIdxGlobal4 = problem_.variables().index(element4);
            int eIdxGlobal5 = problem_.variables().index(element5);
            int eIdxGlobal6 = problem_.variables().index(element6);
            int eIdxGlobal7 = problem_.variables().index(element7);
            int eIdxGlobal8 = problem_.variables().index(element8);

            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 0), eIdxGlobal1, interactionVolume.getIndexOnElement(0, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 1), eIdxGlobal1, interactionVolume.getIndexOnElement(0, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(0, 2), eIdxGlobal1, interactionVolume.getIndexOnElement(0, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 0), eIdxGlobal2, interactionVolume.getIndexOnElement(1, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 1), eIdxGlobal2, interactionVolume.getIndexOnElement(1, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(1, 2), eIdxGlobal2, interactionVolume.getIndexOnElement(1, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 0), eIdxGlobal3, interactionVolume.getIndexOnElement(2, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 1), eIdxGlobal3, interactionVolume.getIndexOnElement(2, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(2, 2), eIdxGlobal3, interactionVolume.getIndexOnElement(2, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 0), eIdxGlobal4, interactionVolume.getIndexOnElement(3, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 1), eIdxGlobal4, interactionVolume.getIndexOnElement(3, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(3, 2), eIdxGlobal4, interactionVolume.getIndexOnElement(3, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 0), eIdxGlobal5, interactionVolume.getIndexOnElement(4, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 1), eIdxGlobal5, interactionVolume.getIndexOnElement(4, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(4, 2), eIdxGlobal5, interactionVolume.getIndexOnElement(4, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 0), eIdxGlobal6, interactionVolume.getIndexOnElement(5, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 1), eIdxGlobal6, interactionVolume.getIndexOnElement(5, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(5, 2), eIdxGlobal6, interactionVolume.getIndexOnElement(5, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 0), eIdxGlobal7, interactionVolume.getIndexOnElement(6, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 1), eIdxGlobal7, interactionVolume.getIndexOnElement(6, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(6, 2), eIdxGlobal7, interactionVolume.getIndexOnElement(6, 2));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 0), eIdxGlobal8, interactionVolume.getIndexOnElement(7, 0));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 1), eIdxGlobal8, interactionVolume.getIndexOnElement(7, 1));
            addRealFluxFaceArea_(interactionVolume.getFaceArea(7, 2), eIdxGlobal8, interactionVolume.getIndexOnElement(7, 2));
        }

    }
}
} // end namespace Dumux
#endif
