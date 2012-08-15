/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
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
#ifndef DUMUX_FVMPFAOINTERACTIONVOLUME_HH
#define DUMUX_FVMPFAOINTERACTIONVOLUME_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA method that does not change with time
 * @author Markus Wolff
 */

#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>

namespace Dumux
{

template<class TypeTag>
class FVMPFAOInteractionVolume
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<DimVector, dim> FieldVectorVector;
    typedef Dune::FieldVector<int, dim> IndexVector;
    typedef std::vector<BoundaryTypes> BCTypeVector;
    typedef std::vector<PrimaryVariables> BCVector;

public:
    enum FaceTypes
    {
        inside = 1,
        boundary = 0,
        outside = -1
    };

    //! Constructs a InteractionVolumeInfo object
    /**
     */
    FVMPFAOInteractionVolume() :
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

    void setStored()
    {
        stored_ = true;
    }

    bool isStored() const
    {
        return stored_;
    }

    void setSubVolumeElement(ElementPointer pointer, int subVolumeIdx)
    {
        elements_[subVolumeIdx].push_back(pointer);
    }

    void setDF(Scalar dF, int subVolumeIdx)
    {
        dF_[subVolumeIdx] = dF;
    }

    void setFaceArea(Scalar& faceArea, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        faceArea_[subVolumeIdx][subVolumeFaceIdxInInside] = faceArea;
    }

    void setPermTimesNu(DimVector& nu, DimMatrix& perm, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInside] = 0;

        perm.mv(nu, permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInside]);

        nu_[subVolumeIdx][subVolumeFaceIdxInInside] = nu;
    }

    void setNormal(DimVector& normal, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        normal_[subVolumeIdx][subVolumeFaceIdxInInside] = normal;
    }

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

    void setOutsideFace(int subVolumeFaceIdx)
    {
        faceType_[subVolumeFaceIdx] = outside;
    }

    void setDirichletCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (dirichletValues_.size() == 0)
        {
            dirichletValues_.resize(2 * dim);
        }
        dirichletValues_[subVolumeFaceIdx] = condition;
    }

    void setNeumannCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (neumannValues_.size() == 0)
        {
            neumannValues_.resize(2 * dim);
        }
        neumannValues_[subVolumeFaceIdx] = condition;
    }

    void setIndexOnElement(int indexInInside, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        indexOnElement_[subVolumeIdx][subVolumeFaceIdxInInside] = indexInInside;
    }

    int getIndexOnElement(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return indexOnElement_[subVolumeIdx][subVolumeFaceIdx];
    }

    int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return faceIndexOnSubVolume_[subVolumeIdx][subVolumeFaceIdx];
    }

    ElementPointer& getSubVolumeElement(int subVolumeIdx)
    {
        return elements_[subVolumeIdx][0];
    }

    BoundaryTypes& getBoundaryType(int subVolumeFaceIdx)
    {
        return boundaryTypes_[subVolumeFaceIdx];
    }

    bool isOutsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == outside);
    }

    bool isInsideFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == inside);
    }

    bool isInnerVolume()
    {
        for (int i = 0; i < faceType_.size(); i++)
        {
            if (isOutsideFace(i) || isBoundaryFace(i))
                return false;
        }
        return true;
    }

    bool isBoundaryFace(int subVolumeFaceIdx)
    {
        return (faceType_[subVolumeFaceIdx] == boundary);
    }

    PrimaryVariables& getDirichletValues(int subVolumeFaceIdx)
    {
        return dirichletValues_[subVolumeFaceIdx];
    }

    PrimaryVariables& getNeumannValues(int subVolumeFaceIdx)
    {
        return neumannValues_[subVolumeFaceIdx];
    }

    DimVector& getNormal(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return normal_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceArea_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    DimVector& getNu(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return nu_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //returns n^TKK_rnu
    Scalar getNTKNu(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return normal_[subVolumeIdx][subVolumeFaceIdxInInsideN] * permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInsideNu];
    }

    //returns n^TKK_rnu
    Scalar getNTKrKNu(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        DimVector krKNu(permTimesNu_[subVolumeIdx][subVolumeFaceIdxInInsideNu]);

        krKNu *= relPerm;

        return normal_[subVolumeIdx][subVolumeFaceIdxInInsideN] * krKNu;
    }

    //returns n^TKK_rnu/dF
    Scalar getNTKNu_by_dF(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return  faceArea_[subVolumeIdx][subVolumeFaceIdxInInsideN]*getNTKNu(subVolumeIdx, subVolumeFaceIdxInInsideN, subVolumeFaceIdxInInsideNu) / dF_[subVolumeIdx];
    }

    //returns n^TKK_rnu/dF
    Scalar getNTKrKNu_by_dF(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return  faceArea_[subVolumeIdx][subVolumeFaceIdxInInsideN]*getNTKrKNu(relPerm, subVolumeIdx, subVolumeFaceIdxInInsideN, subVolumeFaceIdxInInsideNu) / dF_[subVolumeIdx];
    }

private:
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
    std::vector<std::vector<ElementPointer> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
};
}
#endif
