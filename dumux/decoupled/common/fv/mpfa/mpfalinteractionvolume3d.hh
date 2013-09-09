/*****************************************************************************
 *   Copyright (C) 2011 by Yufei Cao                                         *
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_FVMPFAL3DINTERACTIONVOLUME_HH
#define DUMUX_FVMPFAL3DINTERACTIONVOLUME_HH

#include "fvmpfaproperties.hh"

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA 3D method that does not change with time
 * @author Markus Wolff
 */

namespace Dumux
{

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

    static int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceIndexOnSubVolume_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    static int getEdgeIndexFromSubVolumeFace(int subVolumeIdx, int subVolumeFaceIdxInInside, int edgeIdx)
    {
        return edgeIndexOnSubVolumeFace_[subVolumeIdx][subVolumeFaceIdxInInside][edgeIdx];
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

template<class TypeTag>
class FvMpfaL3dInteractionVolume
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum
    {
        dim = GridView::dimension, 
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<DimVector, dim> FieldVectorVector;
    typedef Dune::FieldVector<DimVector, 2> FieldVectorVector2;
    typedef Dune::FieldVector<FieldVectorVector2, dim> FieldVectorVectorVector;
    typedef Dune::FieldVector<int, dim> IndexVector;
    typedef std::vector<BoundaryTypes> BCTypeVector;
    typedef std::vector<PrimaryVariables> BCVector;

public:
    enum FaceTypes
    {
        inside = 1,
        boundary = 0,
        outside = -1,
    };

    enum
    {
        subVolumeTotalNum = IndexTranslator::subVolumeTotalNum,
        fluxFacesTotalNum = IndexTranslator::fluxFacesTotalNum,
        fluxEdgesTotalNum = IndexTranslator::fluxEdgesTotalNum
    };

    //! Constructs a FvMpfaL3dInteractionVolumeAdaptiveInfo object
    /**
     */
    FvMpfaL3dInteractionVolume() :
        normal_(FieldVectorVector(DimVector(0.0))),
        facePos_(DimVector(0.0)),
        edgePos_((DimVector(0.0))),
        faceArea_(0.0),
        faceType_(fluxFacesTotalNum, inside),
        indexOnElement_(IndexVector(0.0)),
        elements_(subVolumeTotalNum),
        centerVertexPos_(0),
        elementNum_(0)
    {}

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

    void setCenterPosition(const DimVector &centerVertexPos)
    {
        centerVertexPos_ = centerVertexPos;
    }

    void setSubVolumeElement(ElementPointer pointer, int subVolumeIdx)
    {
        if (!hasSubVolumeElement(subVolumeIdx))
        {
            elements_[subVolumeIdx].push_back(pointer);
            elementNum_++;
        }
        else
        {
            elements_[subVolumeIdx].insert(elements_[subVolumeIdx].begin(), pointer);
        }
    }

    void setFacePosition(const DimVector& pos, int faceIdx)
    {
        facePos_[faceIdx] = pos;
    }

    void setEdgePosition(const DimVector& pos, int edgeIdx)
    {
        edgePos_[edgeIdx] = pos;
    }

    void setFaceArea(Scalar faceArea, int faceIdx)
    {
        faceArea_[faceIdx] = faceArea;
    }

    void setNormal(DimVector& normal, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        normal_[subVolumeIdx][subVolumeFaceIdxInInside] = normal;
    }

    void setBoundary(BoundaryTypes& boundaryTypes, int subVolumeFaceIdx)
    {
        if (boundaryTypes_.size() == 0)
        {
            boundaryTypes_.resize(fluxFacesTotalNum);
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
            dirichletValues_.resize(fluxFacesTotalNum, PrimaryVariables(0.0));
        }
        dirichletValues_[subVolumeFaceIdx] = condition;
    }

    void setNeumannCondition(PrimaryVariables& condition, int subVolumeFaceIdx)
    {
        if (neumannValues_.size() == 0)
        {
            neumannValues_.resize(fluxFacesTotalNum, PrimaryVariables(0.0));
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

    int getFaceIndexFromSubVolume(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdxInInside);
    }

    int getEdgeIndexFromSubVolumeFace(int subVolumeIdx, int subVolumeFaceIdxInInside, int edgeIdx)
    {
        return IndexTranslator::getEdgeIndexFromSubVolumeFace(subVolumeIdx, subVolumeFaceIdxInInside, edgeIdx);
    }

    int getElementNumber()
    {
        return elementNum_;
    }

    int getIndexOnElement(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return indexOnElement_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    ElementPointer& getSubVolumeElement(int subVolumeIdx)
    {
        if (hasSubVolumeElement(subVolumeIdx))
            return elements_[subVolumeIdx][0];
        else
        {
            std::cout<<"Problems when calling getSubVolumeElement("<<subVolumeIdx<<")\n";
            printInteractionVolumeInfo();
            DUNE_THROW(Dune::RangeError, "element not in interaction volume!");
        }
    }

    bool hasSubVolumeElement(int subVolumeIdx)
    {
        return (elements_[subVolumeIdx].size() > 0);
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

    bool isBoundaryInteractionVolume()
    {
        return (boundaryTypes_.size() > 0);
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
        return facePos_[IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdxInInside)];
    }

    DimVector& getEdgePosition(int subVolumeIdx, int subVolumeFaceIdxInInside, int subVolumeEdgeIdxInInside)
    {
        return edgePos_[IndexTranslator::getEdgeIndexFromSubVolumeFace(subVolumeIdx, subVolumeFaceIdxInInside, subVolumeEdgeIdxInInside)];
    }    

    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceArea_[IndexTranslator::getFaceIndexFromSubVolume(subVolumeIdx, subVolumeFaceIdxInInside)];
    }

    DimVector& getFacePosition(int faceIdx)
    {
        return facePos_[faceIdx];
    }

    DimVector& getEdgePosition(int edgeIdx)
    {
        return edgePos_[edgeIdx];
    }

    Scalar& getFaceArea(int faceIdx)
    {
        return faceArea_[faceIdx];
    }

    void printInteractionVolumeInfo()
    {
        std::cout<<"\nNumber of stored elements: "<<elementNum_<<"\n";
        std::cout<<" center position: "<<centerVertexPos_<<"\n";
        for (int i = 0; i < subVolumeTotalNum; i++)
        {
            if (elements_[i].size() > 0)
            {
            std::cout<<"element "<<i<<":\n";
            std::cout<<"element level: "<<elements_[i][0].level()<<"\n";
            std::cout<<"element position: "<<elements_[i][0]->geometry().center()<<"\n";
            std::cout<<"element volume: "<<elements_[i][0]->geometry().volume()<<"\n";
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
//            dataFile<<"element level: "<<elements_[i][0].level()<<"\n";
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
    Dune::FieldVector<FieldVectorVector, subVolumeTotalNum> normal_;
    Dune::FieldVector<DimVector, fluxFacesTotalNum> facePos_;
    Dune::FieldVector<DimVector, fluxEdgesTotalNum> edgePos_;
    Dune::FieldVector<Scalar, fluxFacesTotalNum> faceArea_;
    BCTypeVector boundaryTypes_;
    std::vector<int> faceType_;
    Dune::FieldVector<IndexVector, subVolumeTotalNum> indexOnElement_;
    std::vector<std::vector<ElementPointer> > elements_;
    BCVector neumannValues_;
    BCVector dirichletValues_;
    DimVector centerVertexPos_;
    int elementNum_;
};
}
#endif
