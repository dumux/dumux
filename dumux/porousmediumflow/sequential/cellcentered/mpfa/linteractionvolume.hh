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

#ifndef DUMUX_FVMPFALINTERACTIONVOLUME_HH
#define DUMUX_FVMPFALINTERACTIONVOLUME_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA L-method that does not change with time.
 */

#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

namespace Dumux
{

//! \ingroup IMPET
/*! \brief Class including the information of an interaction volume of a MPFA L-method that does not change with time.
 *
 * Includes information needed to calculate the transmissibility matrix of an L-interaction-volume.
 *
 */
template<class TypeTag>
class FVMPFALInteractionVolume
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
        {
            dim = GridView::dimension
        };

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename Grid::template Codim<0>::EntitySeed;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using DimVector = Dune::FieldVector<Scalar, dim>;
    using FieldVectorVector = Dune::FieldVector<DimVector, dim>;
    using IndexVector = Dune::FieldVector<int, dim>;
    using BCTypeVector = std::vector<BoundaryTypes>;
    using BCVector = std::vector<PrimaryVariables>;

public:
    enum FaceTypes
        {
            inside = 1,
            boundary = 0,
            outside = -1
        };

    //! Constructs a FVMPFALInteractionVolume object
    FVMPFALInteractionVolume(const Grid& grid)
    : grid_(&grid),
        stored_(false), normal_(FieldVectorVector(DimVector(0.0))),
        facePos_(FieldVectorVector(DimVector(0.0))),
        faceArea_(DimVector(0.0)),
        faceType_(2*dim, inside),
        indexOnElement_(IndexVector(0.0)),
        elements_(2*dim),
        centerVertexPos_(0),
        elementNum_(0)
    {
        faceIndexOnSubVolume_[0][0] = 0;
        faceIndexOnSubVolume_[0][1] = 3;
        faceIndexOnSubVolume_[1][0] = 1;
        faceIndexOnSubVolume_[1][1] = 0;
        faceIndexOnSubVolume_[2][0] = 2;
        faceIndexOnSubVolume_[2][1] = 1;
        faceIndexOnSubVolume_[3][0] = 3;
        faceIndexOnSubVolume_[3][1] = 2;
    }

    //! Delete stored information
    void reset()
    {
        stored_ = false;
        elements_.clear();
        elements_.resize(2*dim);
        centerVertexPos_ = 0;
        indexOnElement_ = IndexVector(0.0);
        faceType_.clear();
        faceType_.resize(2*dim, inside);
        faceArea_ = DimVector(0.0);
        facePos_ = FieldVectorVector(DimVector(0.0));
        normal_ = FieldVectorVector(DimVector(0.0));
        elementNum_ = 0;
        boundaryTypes_.clear();
        neumannValues_.clear();
        dirichletValues_.clear();
    }

    //! Mark storage as completed
    void setStored()
    {
        stored_ = true;
    }

    //! Returns true if information has already been stored
    bool isStored() const
    {
        return stored_;
    }

    //! Store position of the center vertex of the interaction volume
    void setCenterPosition(DimVector &centerVertexPos)
    {
        centerVertexPos_ = centerVertexPos;
    }

    //! Store an element of the interaction volume
    /*!
     *  \param element The element
     *  \param subVolumeIdx The local element index in the interaction volume
     */
    void setSubVolumeElement(const Element& element, int subVolumeIdx)
    {
        elements_[subVolumeIdx].push_back(element.seed());
        elementNum_++;
    }

    //! Store position of the center of a flux face
    /*!
     *  \param pos Position of face center
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setFacePosition(const DimVector& pos, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        facePos_[subVolumeIdx][subVolumeFaceIdxInInside] = pos;
    }

    //! Store a flux face area
    /*!
     *  \param faceArea The flux face area
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setFaceArea(Scalar& faceArea, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        faceArea_[subVolumeIdx][subVolumeFaceIdxInInside] = faceArea;
    }

    //! Store a flux face normal
    /*!
     *  \param normal The normal vector
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setNormal(DimVector& normal, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        normal_[subVolumeIdx][subVolumeFaceIdxInInside] = normal;
    }

    //! Store boundary condtion types for a flux face
    /*!
     *  \param boundaryTypes BoundaryTypes object
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setBoundary(BoundaryTypes& boundaryTypes, int subVolumeFaceIdx)
    {
        if (boundaryTypes_.size() == 0)
        {
            boundaryTypes_.resize(2 * dim);
        }
        boundaryTypes_[subVolumeFaceIdx] = boundaryTypes;
        faceType_[subVolumeFaceIdx] = boundary;
    }

    //! Mark a flux face to be outside the model domain
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setOutsideFace(int subVolumeFaceIdx)
    {
        faceType_[subVolumeFaceIdx] = outside;
    }

    //! Store Dirichlet boundary condtions for a flux face
    /*!
     *  \param condition Vector of primary variables
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setDirichletCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (dirichletValues_.size() == 0)
        {
            dirichletValues_.resize(2 * dim);
        }
        dirichletValues_[subVolumeFaceIdx] = condition;
    }

    //! Store Neumann boundary condtions for a flux face
    /*!
     *  \param condition Vector phase fluxes
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     */
    void setNeumannCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (neumannValues_.size() == 0)
        {
            neumannValues_.resize(2 * dim);
        }
        neumannValues_[subVolumeFaceIdx] = condition;
    }

    //! Store map from local interaction volume numbering to numbering of the Dune reference element.
    /*!
     * \param indexInInside Face index of the Dune reference element
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setIndexOnElement(int indexInInside, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        indexOnElement_[subVolumeIdx][subVolumeFaceIdxInInside] = indexInInside;
    }

    //! Get position vector of central vertex
    DimVector& getCenterPosition()
    {
        return centerVertexPos_;
    }

    //! Get number of stored elements
    int getElementNumber()
    {
        return elementNum_;
    }

    //! Map from local interaction volume numbering to numbering of the Dune reference element.
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

    //! Map from local interaction volume numbering on element to numbering on interaction volume.
    /*!
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdx The local face index in the interaction volume element
     *
     * \return Local face index int the interaction volume
     */
    int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return faceIndexOnSubVolume_[subVolumeIdx][subVolumeFaceIdx];
    }

    //! Get an element of the interaction volume.
    /*!
     * \param subVolumeIdx The local element index in the interaction volume
     *
     * \return The interaction volume sub-element.
     */
    Element getSubVolumeElement(int subVolumeIdx)
    {
        return grid_->entity(elements_[subVolumeIdx][0]);
    }

    //! Get boundary condtion types for a flux face
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Object containing information about the boundary types.
     */
    BoundaryTypes& getBoundaryType(int subVolumeFaceIdx)
    {
        return boundaryTypes_[subVolumeFaceIdx];
    }

    //! Returns true if an interaction volume flux face is outside the model domain.
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     */
    bool isOutsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == outside);
    }

    //! Returns true if an interaction volume flux face is inside the model domain and no boundary face.
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     */
    bool isInsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == inside);
    }

    //! Returns true if the interaction volume is completely inside the model domain.
    bool isInnerVolume()
    {
        for (unsigned int i = 0; i < faceType_.size(); i++)
        {
            if (isOutsideFace(i) || isBoundaryFace(i))
                return false;
        }
        return true;
    }

    //! Returns true if an interaction volume flux face is a boundary face.
    /*!
     *   \param subVolumeFaceIdx The local face index in the interaction volume
     */
    bool isBoundaryFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == boundary);
    }

    //! Get the Dirichlet boundary condtions for a flux face
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Vector of primary variables.
     */
    PrimaryVariables& getDirichletValues(int subVolumeFaceIdx)
    {
        return dirichletValues_[subVolumeFaceIdx];
    }

    //! Get the Neumann boundary condtions for a flux face
    /*!
     *  \param subVolumeFaceIdx The local face index in the interaction volume
     *
     *  \return Vector of phase fluxes.
     */
    PrimaryVariables& getNeumannValues(int subVolumeFaceIdx)
    {
        return neumannValues_[subVolumeFaceIdx];
    }

    //! Get a flux face normal
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     *
     *  \return Normal vector
     */
    DimVector& getNormal(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return normal_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //! Get position of the center of a flux face
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     *
     *  \return Position of face center
     */
    DimVector& getFacePosition(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return facePos_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //! Get a flux face area
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     *
     *  \return Area of the flux face
     */
    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceArea_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //! Outputs stored information
    void printInteractionVolumeInfo()
    {
        std::cout<<"\nNumber of stored elements: "<<elementNum_<<"\n";
        std::cout<<" center position: "<<centerVertexPos_<<"\n";
        for (int i = 0; i < 2*dim; i++)
        {
            if (elements_[i].size() > 0)
            {
                std::cout<<"element "<<i<<":\n";
                std::cout<<"element position: "<<elements_[i][0]->geometry().center()<<"\n";
                std::cout<<"face indices on element: "<<indexOnElement_[i]<<"\n";
                std::cout<<"face normals on element: "<<normal_[i]<<"\n";
                std::cout<<"face areas on element: "<<faceArea_[i]<<"\n";
                std::cout<<"face position on element: "<<facePos_[i]<<"\n";
                std::cout<<"face type: "<<faceType_[i]<<"\n";
            }
        }
    }

private:
    const Grid* grid_;
    bool stored_;
    Dune::FieldVector<FieldVectorVector, 2*dim> normal_;
    Dune::FieldVector<FieldVectorVector, 2*dim> facePos_;
    Dune::FieldVector<DimVector, 2*dim> faceArea_;
    BCTypeVector boundaryTypes_;
    std::vector<int> faceType_;
    Dune::FieldVector<IndexVector, 2*dim> indexOnElement_;
    Dune::FieldVector<IndexVector, 2*dim> faceIndexOnSubVolume_;
    std::vector<std::vector<ElementSeed> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
    DimVector centerVertexPos_;
    int elementNum_;
};
}
#endif
