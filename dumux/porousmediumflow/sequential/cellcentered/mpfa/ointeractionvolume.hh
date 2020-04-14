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
#ifndef DUMUX_FVMPFAOINTERACTIONVOLUME_HH
#define DUMUX_FVMPFAOINTERACTIONVOLUME_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA O-method that does not change with time.
 */

#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

namespace Dumux
{

//! \ingroup IMPET
/*! \brief Class including the information of an interaction volume of a MPFA O-method that does not change with time.
 *
 * Includes information needed to calculate the transmissibility matrix of an O-interaction-volume.
 *
 */
template<class TypeTag>
class FVMPFAOInteractionVolume
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
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
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

    //! Constructs a FVMPFAOInteractionVolume object
    FVMPFAOInteractionVolume(const Grid& grid)
    : grid_(&grid),
        stored_(false), permTimesNu_(FieldVectorVector(DimVector(0.0))),
        nu_(FieldVectorVector(DimVector(0.0))), normal_(FieldVectorVector(DimVector(0.0))),
        faceArea_(DimVector(0.0)), dF_(0.0),
        faceType_(2*dim, inside),
        indexOnElement_(IndexVector(0.0)),
        elements_(2*dim)
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

    //! Store an element of the interaction volume
    /*!
     *  \param element The element
     *  \param subVolumeIdx The local element index in the interaction volume
     */
    void setSubVolumeElement(const Element& element, int subVolumeIdx)
    {
        elements_[subVolumeIdx].push_back(element.seed());
    }

    //! Store the \f$ dF \f$ for the transmissiblity calculation
    /*!
     *  \param dF Value of  \f$ dF \f$
     *  \param subVolumeIdx The local element index in the interaction volume
     */
    void setDF(Scalar dF, int subVolumeIdx)
    {
        dF_[subVolumeIdx] = dF;
    }

    //! Store a flux face area
    /*!
     *  \param faceArea The flux face area
     *  \param subVolumeIdx The local element index in the interaction volume
     *   \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setFaceArea(Scalar& faceArea, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        faceArea_[subVolumeIdx][subVolumeFaceIdxInInside] = faceArea;
    }

    //! Store \f$ \boldsymbol K \boldsymbol \nu\f$  for the transmissiblity calculation
    /*!
     *  \param nu \f$ \boldsymbol \nu \f$ vector
     *  \param perm Permeability matrix
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     */
    void setPermTimesNu(DimVector& nu, DimMatrix& perm, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInside] = 0;

        perm.mv(nu, permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInside]);

        nu_[subVolumeIdx][subVolumeFaceIdxInInside] = nu;
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
//        std::cout<<"old BCType = "<<boundaryType_[subVolumeFaceIdx]<<"\n";
//        std::cout<<"new BCType = "<<boundaryType<<"\n";
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

    //! Get \f$ \boldsymbol \nu \f$ vector used for the transmissiblity calculation
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInside The local face index in the interaction volume element
     *
     *  \return \f$ \boldsymbol \nu \f$ vector
     */
    DimVector& getNu(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return nu_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //! Get \f$ \boldsymbol n^\text{T} \boldsymbol K \boldsymbol \nu \f$  for the transmissiblity calculation
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInsideN The local face index in the interaction volume element
     *                                   for the normal \f$ \boldsymbol n \f$
     *  \param subVolumeFaceIdxInInsideNu The local face index in the interaction volume element
     *                                   for the vector \f$ \boldsymbol \nu \f$
     *
     *  \return \f$ \boldsymbol n^\text{T} \boldsymbol K \boldsymbol \nu \f$
     */
    Scalar getNtkNu(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return normal_[subVolumeIdx][subVolumeFaceIdxInInsideN] * permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInsideNu];
    }

    //! Get \f$ \boldsymbol n^\text{T} k_{r\alpha} \boldsymbol K \boldsymbol \nu\f$  for the transmissiblity calculation
    /*!
     * \param relPerm relative permeability value (\f$ \boldsymbol n^\text{T} k_{r\alpha} \f$)
     * \param subVolumeIdx The local element index in the interaction volume
     * \param subVolumeFaceIdxInInsideN The local face index in the interaction volume element
     *                                  for the normal \f$ \boldsymbol n \f$
     * \param subVolumeFaceIdxInInsideNu The local face index in the interaction volume element
     *                                  for the vector \f$ \boldsymbol \nu \f$
     *
     *  \return \f$ \boldsymbol n^\text{T} k_{r\alpha} \boldsymbol K \boldsymbol \nu \f$
     */
    Scalar getNtkrkNu(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN,
                      int subVolumeFaceIdxInInsideNu) const
    {
        DimVector krKNu(permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInsideNu]);

        krKNu *= relPerm;

        return normal_[subVolumeIdx][subVolumeFaceIdxInInsideN] * krKNu;
    }

    //! Get \f$ \frac{1}{dF} \left(\boldsymbol n^\text{T} \boldsymbol K \boldsymbol \nu \right) \f$
    //! for the transmissiblity calculation
    /*!
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInsideN The local face index in the interaction volume element
     *                                   for the normal \f$ \boldsymbol n \f$
     *  \param subVolumeFaceIdxInInsideNu The local face index in the interaction volume element
     *                                   for the vector \f$ \boldsymbol \nu \f$
     *
     *  \return \f$ \frac{1}{dF} \left(\boldsymbol n^\text{T} \boldsymbol K \boldsymbol \nu \right) \f$
     */
    Scalar getNtkNu_df(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return  faceArea_[subVolumeIdx][subVolumeFaceIdxInInsideN]
                * getNtkNu(subVolumeIdx, subVolumeFaceIdxInInsideN, subVolumeFaceIdxInInsideNu)
                / dF_[subVolumeIdx];
    }

    //! Get \f$ \frac{1}{dF} \left(\boldsymbol n^\text{T} k_{r\alpha} \boldsymbol K \boldsymbol \nu \right) \f$
    //! for the transmissiblity calculation
    /*!
     *  \param relPerm relative permeability value (\f$ \boldsymbol n^\text{T} k_{r\alpha} \f$)
     *  \param subVolumeIdx The local element index in the interaction volume
     *  \param subVolumeFaceIdxInInsideN The local face index in the interaction volume element
     *                                   for the normal \f$ \boldsymbol n \f$
     *  \param subVolumeFaceIdxInInsideNu The local face index in the interaction volume element
     *                                    for the vector \f$ \boldsymbol \nu \f$
     *
     *  \return \f$ \frac{1}{dF} \left(\boldsymbol n^\text{T} k_{r\alpha} \boldsymbol K \boldsymbol \nu \right) \f$
     */
    Scalar getNtkrkNu_df(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return  faceArea_[subVolumeIdx][subVolumeFaceIdxInInsideN]
                * getNtkrkNu(relPerm, subVolumeIdx, subVolumeFaceIdxInInsideN, subVolumeFaceIdxInInsideNu)
                / dF_[subVolumeIdx];
    }

private:
    const Grid* grid_;
    bool stored_;
    Dune::FieldVector<FieldVectorVector, 2*dim> permTimesNu_;
    Dune::FieldVector<FieldVectorVector, 2*dim> nu_;
    Dune::FieldVector<FieldVectorVector, 2*dim> normal_;
    Dune::FieldVector<DimVector, 2*dim> faceArea_;
    Dune::FieldVector<Scalar, 2*dim> dF_;
    BCTypeVector boundaryTypes_;
    std::vector<int> faceType_;
    Dune::FieldVector<IndexVector, 2*dim> indexOnElement_;
    Dune::FieldVector<IndexVector, 2*dim> faceIndexOnSubVolume_;
    std::vector<std::vector<ElementSeed> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
};
}
#endif
