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
#ifndef DUMUX_FVMPFALINTERACTIONVOLUME_HH
#define DUMUX_FVMPFALINTERACTIONVOLUME_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA method that does not change with time
 * @author Markus Wolff
 */

namespace Dumux
{

template<class TypeTag>
class FVMPFALInteractionVolume
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

    //! Constructs a FVMPFALInteractionVolumeInfo object
    /**
     */
    FVMPFALInteractionVolume() :
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
    }

    void setStored()
    {
        stored_ = true;
    }

    bool isStored() const
    {
        return stored_;
    }

    void setCenterPosition(DimVector &centerVertexPos)
    {
        centerVertexPos_ = centerVertexPos;
    }

    void setSubVolumeElement(ElementPointer pointer, int subVolumeIdx)
    {
        elements_[subVolumeIdx].push_back(pointer);
        elementNum_++;
    }

    void setFacePosition(const DimVector& pos, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        facePos_[subVolumeIdx][subVolumeFaceIdxInInside] = pos;
    }

    void setFaceArea(Scalar& faceArea, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        faceArea_[subVolumeIdx][subVolumeFaceIdxInInside] = faceArea;
    }

    void setNormal(DimVector& normal, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        normal_[subVolumeIdx][subVolumeFaceIdxInInside] = normal;
    }

    void setBoundary(BoundaryTypes& boundaryTypes, int subVolumeFaceIdx)
    {
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

    DimVector& getCenterPosition()
    {
        return centerVertexPos_;
    }

    void setIndexOnElement(int indexInInside, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        indexOnElement_[subVolumeIdx][subVolumeFaceIdxInInside] = indexInInside;
    }

    int getElementNumber()
    {
        return elementNum_;
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

    DimVector& getFacePosition(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return facePos_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceArea_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

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
    bool stored_;
    Dune::FieldVector<FieldVectorVector, 2*dim> normal_;
    Dune::FieldVector<FieldVectorVector, 2*dim> facePos_;
    Dune::FieldVector<DimVector, 2*dim> faceArea_;
    BCTypeVector boundaryTypes_;
    std::vector<int> faceType_;
    Dune::FieldVector<IndexVector, 2*dim> indexOnElement_;
    Dune::FieldVector<IndexVector, 2*dim> faceIndexOnSubVolume_;
    std::vector<std::vector<ElementPointer> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
    DimVector centerVertexPos_;
    int elementNum_;
};
}
#endif
