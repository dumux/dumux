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

#ifndef DUMUX_FVMPFAL3DINTERACTIONVOLUME_HH
#define DUMUX_FVMPFAL3DINTERACTIONVOLUME_HH

#include "properties.hh"

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA 3D method that does not change with time
 * @author Markus Wolff
 */

namespace Dumux
{

//! \cond \private
// Mapper for local interaction volume indices (see doc/docextra/3dmpfa).
class IndexTranslator
{
public:
    enum
    {
        subVolumeTotalNum = 8,
        fluxFacesTotalNum = 12,
        fluxFacesNumOnSubVolume = 3,
        fluxEdgesTotalNum = 6,
        edgesNumOnFluxFace = 2
    };

    static int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return faceIndexOnSubVolume_[subVolumeIdx][subVolumeFaceIdx];
    }

    static int getEdgeIndexFromSubVolumeFace(int subVolumeIdx, int subVolumeFaceIdx, int subVolumeEdgeIdx)
    {
        return edgeIndexOnSubVolumeFace_[subVolumeIdx][subVolumeFaceIdx][subVolumeEdgeIdx];
    }

    static int getFaceIndexFromElements(int elem1Idx, int elem2Idx)
    {
        return faceIndexFromElements_[elem1Idx][elem2Idx];
    }


private:
    static const int faceIndexOnSubVolume_[subVolumeTotalNum][fluxFacesNumOnSubVolume];
    static const int edgeIndexOnSubVolumeFace_[subVolumeTotalNum][fluxFacesNumOnSubVolume][edgesNumOnFluxFace];
    static const int faceIndexFromElements_[subVolumeTotalNum][subVolumeTotalNum];
};

const int IndexTranslator::faceIndexOnSubVolume_[subVolumeTotalNum][fluxFacesNumOnSubVolume] =
{
        {0, 3, 8},
        {1, 0, 9},
        {3, 2, 11},
        {2, 1, 10},
        {8, 4, 7},
        {9, 5, 4},
        {11, 7, 6},
        {10, 6, 5}
};

const int IndexTranslator::edgeIndexOnSubVolumeFace_[subVolumeTotalNum][fluxFacesNumOnSubVolume][edgesNumOnFluxFace] =
{
        {{2, 0},{0, 5},{5, 2}},
        {{3, 0},{0, 2},{2, 3}},
        {{5, 0},{0, 4},{4, 5}},
        {{4, 0},{0, 3},{3, 4}},
        {{2, 5},{1, 2},{5, 1}},
        {{3, 2},{1, 3},{2, 1}},
        {{5, 4},{1, 5},{4, 1}},
        {{4, 3},{1, 4},{3, 1}}
};

const int IndexTranslator::faceIndexFromElements_[subVolumeTotalNum][subVolumeTotalNum] =
{
        {-1, 0, 3, -1, 8, -1, -1, -1},
        {0, -1, -1, 1, -1, 9, -1, -1},
        {3, -1, -1, 2, -1, -1, 11, -1},
        {-1, 1, 2, -1, -1, -1, -1, 10},
        {8, -1, -1, -1, -1, 4, 7, -1},
        {-1, 9, -1, -1, 4, -1, -1, 5},
        {-1, -1, 11, -1, 7, -1, -1, 6},
        {-1, -1, -1, 10, -1, 5, 6, -1}
};
//! \endcond

//! \ingroup IMPET mpfa
/*! \brief Class including the information of a 3d interaction volume of a MPFA L-method that does not change with time.
 *
 * Includes information needed to calculate the transmissibility matrices of an L-interaction-volume.
 *
 */
template<class TypeTag>
class FvMpfaL3dInteractionVolume
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename Grid::template Codim<0>::EntitySeed;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using DimVector = Dune::FieldVector<Scalar, dim>;
    using FieldVectorVector = Dune::FieldVector<DimVector, dim>;
    using FieldVectorVector2 = Dune::FieldVector<DimVector, 2>;
    using FieldVectorVectorVector = Dune::FieldVector<FieldVectorVector2, dim>;
    using IndexVector = Dune::FieldVector<int, dim>;
    using BCTypeVector = std::vector<BoundaryTypes>;
    using BCVector = std::vector<PrimaryVariables>;

public:
    //An interaction volume around a vertex includes in general 12 flux faces (see doc/docextra/3dmpfa).
    //If the vertex is on the boundary these flux faces can be inside the domain on the boundary or outside the domain.
    enum FaceTypes
    {
        inside = 1,//!< Flux face is inside the model domain
        boundary = 0,//!< Flux face is a boundary face
        outside = -1,//!< Flux face is outside the model domain
    };

    enum
    {
        subVolumeTotalNum = IndexTranslator::subVolumeTotalNum,//!< Number of sub-volumes in the interaction volume
        fluxFacesTotalNum = IndexTranslator::fluxFacesTotalNum,//!< Number of flux faces in the interaction volume
        fluxEdgesTotalNum = IndexTranslator::fluxEdgesTotalNum//!< Number of edges in the interaction volume
    };

    //! Constructs a FvMpfaL3dInteractionVolume object
    /**
     */
    FvMpfaL3dInteractionVolume(const Grid& grid)
    : grid_(&grid)
    , normal_(FieldVectorVector(DimVector(0.0)))
    , facePos_(DimVector(0.0))
    , edgePos_((DimVector(0.0)))
    , faceArea_(0.0)
    , faceType_(fluxFacesTotalNum, inside)
    , indexOnElement_(IndexVector(0.0))
    , elements_(subVolumeTotalNum)
    , centerVertexPos_(0)
    , elementNum_(0)
    {}

    //! Reset the interaction volume (deletes stored data)
    void reset()
    {
        elements_.clear();
        elements_.resize(subVolumeTotalNum);
        centerVertexPos_ = 0;
        indexOnElement_ = IndexVector(0.0);
        faceType_.clear();
        faceType_.resize(fluxFacesTotalNum, inside);
        faceArea_ = 0.0;
        facePos_ = DimVector(0.0);
        edgePos_ = DimVector(0.0);
        normal_ = FieldVectorVector(DimVector(0.0));
        elementNum_ = 0;
        boundaryTypes_.clear();
        neumannValues_.clear();
        dirichletValues_.clear();
    }

    //! Store position of the central vertex
    /*
     * \param centerVertexPos Position of the central vertex
     */
    void setCenterPosition(const DimVector &centerVertexPos)
    {
        centerVertexPos_ = centerVertexPos;
    }
    //! Store a dune element as a sub volume element
    /*!
     *  \param element The element
     *  \param subVolumeIdx The local element index in the interaction volume
     */
    void setSubVolumeElement(const Element& element, int subVolumeIdx)
    {
        if (!hasSubVolumeElement(subVolumeIdx))
        {
            elements_[subVolumeIdx].push_back(element.seed());
            elementNum_++;
        }
        else
        {
            elements_[subVolumeIdx].insert(elements_[subVolumeIdx].begin(), element.seed());
        }
    }

    //! Store the position of a flux face
    /*!
     *  \param pos Position of face center
     *  \param fIdx The interaction volume face index
     */
    void setFacePosition(const DimVector& pos, int fIdx)
    {
        facePos_[fIdx] = pos;
    }

    //! Store the center of the edges
    /*!
     *  \param pos Position of face center
     *  \param edgeIdx The interaction volume edge index
     */
    void setEdgePosition(const DimVector& pos, int edgeIdx)
    {
        edgePos_[edgeIdx] = pos;
    }

    //! Store the flux face areas
    /*!
     *  \param faceArea The flux face area
     *  \param fIdx The interaction volume face index
     */
    void setFaceArea(Scalar faceArea, int fIdx)
    {
        faceArea_[fIdx] = faceArea;
    }

    //! Store the normals
    /*!
     *  \param normal The normal vector
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdx The local face index in the interaction volume element
     */
    void setNormal(DimVector& normal, int subVolumeIdx, int subVolumeFaceIdx)
    {
        normal_[subVolumeIdx][subVolumeFaceIdx] = normal;
    }

    //! Store the types of boundary conditions
    /*!
     *  \param boundaryTypes BoundaryTypes object
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setBoundary(BoundaryTypes& boundaryTypes, int subVolumeFaceIdx)
    {
        if (boundaryTypes_.size() == 0)
        {
            boundaryTypes_.resize(fluxFacesTotalNum);
        }
        boundaryTypes_[subVolumeFaceIdx] = boundaryTypes;
        faceType_[subVolumeFaceIdx] = boundary;
    }

    //! Define a flux face to be outside the model domain (for boundary vertices)
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setOutsideFace(int subVolumeFaceIdx)
    {
        faceType_[subVolumeFaceIdx] = outside;
    }

    //! Store Dirichlet boundary conditions
    /*!
     *  \param condition Vector of primary variables
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setDirichletCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (dirichletValues_.size() == 0)
        {
            dirichletValues_.resize(fluxFacesTotalNum, PrimaryVariables(0.0));
        }
        dirichletValues_[subVolumeFaceIdx] = condition;
    }

    //! Store Neumann boundary conditions
    /*!
     *  \param condition Vector phase fluxes
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setNeumannCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (neumannValues_.size() == 0)
        {
            neumannValues_.resize(fluxFacesTotalNum, PrimaryVariables(0.0));
        }
        neumannValues_[subVolumeFaceIdx] = condition;
    }

    //! Store the local dune face index dependent on the local interaction volume indices
    /*!
     * \param indexInInside Face index of the Dune reference element
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdx The local face index in the interaction volume element
     */
    void setIndexOnElement(int indexInInside, int subVolumeIdx, int subVolumeFaceIdx)
    {
        indexOnElement_[subVolumeIdx][subVolumeFaceIdx] = indexInInside;
    }

    //! The position of the central vertex
    DimVector& getCenterPosition()
    {
        return centerVertexPos_;
    }

    //! The local interaction volume face index dependent on the local sub volume indices
    /*
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     * \return The interaction volume face index
     */
    int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdx);
    }

    //! The local interaction volume edge index dependent on the local sub volume indices
    /*
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdx The local face index in the interaction volume element
     * \param subVolumeEdgeIdx The local edge index in the interaction volume element
     *
     * \return The interaction volume edge index
     */
    int getEdgeIndexFromSubVolumeFace(int subVolumeIdx, int subVolumeFaceIdx, int subVolumeEdgeIdx)
    {
        return IndexTranslator::getEdgeIndexFromSubVolumeFace(subVolumeIdx, subVolumeFaceIdx, subVolumeEdgeIdx);
    }

    //! Number of stored dune elements
    int getElementNumber()
    {
        return elementNum_;
    }

    //! The local dune face index dependent on the local interaction volume indices
    /*!
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     * \return Face index of the Dune reference element
     */
    int getIndexOnElement(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return indexOnElement_[subVolumeIdx][subVolumeFaceIdx];
    }

    //! The dune element of the interaction volume sub-volume
    /*!
     * \param subVolumeIdx The local element index in the interaction volume
     *
     * \return The interaction volume sub-element.
     */
    Element getSubVolumeElement(int subVolumeIdx)
    {
        if (hasSubVolumeElement(subVolumeIdx))
            return grid_->entity(elements_[subVolumeIdx][0]);
        else
        {
            std::cout<<"Problems when calling getSubVolumeElement("<<subVolumeIdx<<")\n";
            printInteractionVolumeInfo();
            DUNE_THROW(Dune::RangeError, "element not in interaction volume!");
        }
    }

    //! Check if a dune element is stored for an interaction volume sub-volume
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *
     * \return returns <tt>true</tt> if an element pointer is stored at position  <tt>subVolumeIdx</tt>.
     */
    bool hasSubVolumeElement(int subVolumeIdx)
    {
        return (elements_[subVolumeIdx].size() > 0);
    }

    //! The boundary types
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Object containing information about the boundary types.
     */
    BoundaryTypes& getBoundaryType(int subVolumeFaceIdx)
    {
        return boundaryTypes_[subVolumeFaceIdx];
    }

    //! Check if an interaction volume face is outside the model domain (for boundary vertices)
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *   \return returns <tt>true</tt> if the flux face is outside the model domain.
     */
    bool isOutsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == outside);
    }

    //! Check if an interaction volume face is inside the model domain (for boundary vertices)
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *   \return returns <tt>true</tt> if the flux face is inside the model domain.
     */
    bool isInsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == inside);
    }

    //! Check if an interaction volume face is a model domain boundary (for boundary vertices)
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *   \return returns <tt>true</tt> if the flux face is a boundary face.
     */
    bool isBoundaryFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == boundary);
    }

    //! Check if an interaction volume is not located around a boundary vertex
    /*!
     * \return returns <tt>true</tt> if an interaction volume is not located around a boundary vertex
     */
    bool isInnerVolume()
    {
        for (unsigned int i = 0; i < faceType_.size(); i++)
        {
            if (isOutsideFace(i) || isBoundaryFace(i))
                return false;
        }
        return true;
    }

    //! Check if an interaction volume is located around a boundary vertex
    /*!
     * \return returns <tt>true</tt> if an interaction volume located around a boundary vertex
     */
    bool isBoundaryInteractionVolume()
    {
        return (boundaryTypes_.size() > 0);
    }

    //! The Dirichlet boundary values
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Vector of primary variables.
     */
    PrimaryVariables& getDirichletValues(int subVolumeFaceIdx)
    {
        return dirichletValues_[subVolumeFaceIdx];
    }

    //! The Neumann boundary values
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Vector of phase fluxes.
     */
    PrimaryVariables& getNeumannValues(int subVolumeFaceIdx)
    {
        return neumannValues_[subVolumeFaceIdx];
    }

    //! The face normal
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     *  \return Normal vector
     */
    DimVector& getNormal(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return normal_[subVolumeIdx][subVolumeFaceIdx];
    }

    //! The position of the face center
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     *  \return Position of face center
     */
    DimVector& getFacePosition(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return facePos_[IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdx)];
    }

    //! The position of the edge center
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdx The local face index in the interaction volume element
     *  \param subVolumeEdgeIdx The local edge index in the interaction volume element
     *
     *  \return Position of the edge center
     */
    DimVector& getEdgePosition(int subVolumeIdx, int subVolumeFaceIdx, int subVolumeEdgeIdx)
    {
        return edgePos_[IndexTranslator::getEdgeIndexFromSubVolumeFace(subVolumeIdx, subVolumeFaceIdx, subVolumeEdgeIdx)];
    }

    //! The interaction volume flux face area
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     *  \return Area of the flux face
     */
    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return faceArea_[IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdx)];
    }

    //! The position of the face center
    /*!
     *  \param fIdx The local face index in the interaction volume
     *
     *  \return Position of face center
     */
    DimVector& getFacePosition(int fIdx)
    {
        return facePos_[fIdx];
    }

    //! The position of the edge center
    /*!
     *  \param edgeIdx The local edge index in the interaction volume
     *
     *  \return Position of the edge center
     */
    DimVector& getEdgePosition(int edgeIdx)
    {
        return edgePos_[edgeIdx];
    }

    //! The interaction volume flux face area
    /*!
     *  \param fIdx The local face index in the interaction volume
     *
     *  \return Area of the flux face
     */
    Scalar& getFaceArea(int fIdx)
    {
        return faceArea_[fIdx];
    }

    //! Print the stored interaction volume data
    void printInteractionVolumeInfo()
    {
        std::cout<<"\nNumber of stored elements: "<<elementNum_<<"\n";
        std::cout<<" center position: "<<centerVertexPos_<<"\n";
        for (int i = 0; i < subVolumeTotalNum; i++)
        {
            if (elements_[i].size() > 0)
            {
            std::cout<<"element "<<i<<":\n";
            std::cout<<"element level: "<<grid_->entity(elements_[i][0]).level()<<"\n";
            std::cout<<"element position: "<<grid_->entity(elements_[i][0]).geometry().center()<<"\n";
            std::cout<<"element volume: "<<grid_->entity(elements_[i][0]).geometry().volume()<<"\n";
            std::cout<<"face indices on element: "<<indexOnElement_[i]<<"\n";
            std::cout<<"face normals on element: "<<normal_[i]<<"\n";
            std::cout<<"face areas on element: ";
            for (int j = 0; j < 3; j++)
            {
                std::cout<<getFaceArea(i, j)<<" ";
            }
            std::cout<<"\n";
            std::cout<<"face position on element: ";
            for (int j = 0; j < 3; j++)
            {
                std::cout<<getFacePosition(i, j)<<" ";
            }
            std::cout<<"\n";
            std::cout<<"edgePos:";
            for (int j = 0; j < 3; j++)
            {
                std::cout<<"    ";
                for (int k = 0; k < 2; k++)
                {
                    std::cout<<getEdgePosition(i, j, k) <<" ";
                }
            }
            std::cout<<"\n";
            }
        }
    }

//    void printInteractionVolumeInfoToFile(std::ofstream& dataFile)
//    {
//        dataFile<<"\nNumber of stored elements: "<<elementNum_<<"\n";
//        dataFile<<" center position: "<<centerVertexPos_<<"\n";
//        for (int i = 0; i < subVolumeTotalNum; i++)
//        {
//            if (elements_[i].size() > 0)
//            {
//            dataFile<<"element "<<i<<":\n";
//            dataFile<<"element level: "<<elements_[i][0]->level()<<"\n";
//            dataFile<<"element position: "<<elements_[i][0]->geometry().center()<<"\n";
//            dataFile<<"element volume: "<<elements_[i][0]->geometry().volume()<<"\n";
//            dataFile<<"face indices on element: "<<indexOnElement_[i]<<"\n";
//            dataFile<<"face normals on element: "<<normal_[i]<<"\n";
//            dataFile<<"face areas on element: ";
//            for (int j = 0; j < 3; j++)
//            {
//                dataFile<<getFaceArea(i, j)<<" ";
//            }
//            dataFile<<"\n";
//            dataFile<<"face position on element: ";
//            for (int j = 0; j < 3; j++)
//            {
//                dataFile<<getFacePosition(i, j)<<" ";
//            }
//            dataFile<<"\n";
//            dataFile<<"edgePos:";
//            for (int j = 0; j < 3; j++)
//            {
//                dataFile<<"    ";
//                for (int k = 0; k < 2; k++)
//                {
//                    dataFile<<getEdgePosition(i, j, k) <<" ";
//                }
//            }
//            dataFile<<"\n";
//
//        }
//        dataFile<<"face types: ";
//        for (int i = 0; i < fluxFacesTotalNum; i++)
//        {
//            dataFile<<faceType_[i]<<" ";
//        }
//        dataFile<<"\n";
//        if (isBoundaryInteractionVolume())
//        {
//        dataFile<<"Boundaries:\n";
//        dataFile<<"dirichlet: ";
//        for (int i = 0; i < fluxFacesTotalNum; i++)
//        {
//            dataFile<<boundaryTypes_[i].isDirichlet(0)<<" ";
//            dataFile<<boundaryTypes_[i].isDirichlet(1)<<" ";
//        }
//        dataFile<<"\n";
//        dataFile<<"neumann: ";
//        for (int i = 0; i < fluxFacesTotalNum; i++)
//        {
//            dataFile<<boundaryTypes_[i].isNeumann(0)<<" ";
//            dataFile<<boundaryTypes_[i].isNeumann(1)<<" ";
//        }
//        dataFile<<"\n";
//        dataFile<<"values: ";
//        for (int i = 0; i < fluxFacesTotalNum; i++)
//        {
//            if (dirichletValues_.size() > 0)
//            dataFile<<dirichletValues_[i]<<" ";
//            if (neumannValues_.size() > 0)
//            dataFile<<neumannValues_[i]<<" ";
//        }
//        dataFile<<"\n";
//        }
//    }

private:
    const Grid* grid_;
    Dune::FieldVector<FieldVectorVector, subVolumeTotalNum> normal_;
    Dune::FieldVector<DimVector, fluxFacesTotalNum> facePos_;
    Dune::FieldVector<DimVector, fluxEdgesTotalNum> edgePos_;
    Dune::FieldVector<Scalar, fluxFacesTotalNum> faceArea_;
    BCTypeVector boundaryTypes_;
    std::vector<int> faceType_;
    Dune::FieldVector<IndexVector, subVolumeTotalNum> indexOnElement_;
    std::vector<std::vector<ElementSeed> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
    DimVector centerVertexPos_;
    int elementNum_;
};
}
#endif
