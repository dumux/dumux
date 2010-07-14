// $Id: interactionvolume.hh 3732 2010-06-11 13:27:20Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_INTERACTIONVOLUME_HH
#define DUMUX_INTERACTIONVOLUME_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA method that does not change with time
 * @author Markus Wolff
 */

namespace Dumux
{

template<class DirichletType, class NeumannType, BoundaryConditions::Flags bcType> class ReturnBoundaryCondition
{
};

template<class DirichletType, class NeumannType> class ReturnBoundaryCondition<DirichletType, NeumannType, BoundaryConditions::dirichlet>
{
public:
    typedef DirichletType ReturnType;
    typedef typename ReturnType::value_type BCType;

    static ReturnType& returnBoundaryCondition(DirichletType& dirichletBC, NeumannType& neumannBC)
    {
        return dirichletBC;
    }
};

template<class DirichletType, class NeumannType> class ReturnBoundaryCondition<DirichletType, NeumannType, BoundaryConditions::neumann>
{
public:
    typedef NeumannType ReturnType;
    typedef typename ReturnType::value_type BCType;

    static ReturnType& returnBoundaryCondition(DirichletType& dirichletBC, NeumannType& neumannBC)
    {
        return neumannBC;
    }
};

template<class TypeTag>
class InteractionVolume
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<FieldVector, dim> FieldVectorVector;
    typedef Dune::FieldVector<int, dim> IndexVector;
    typedef std::vector<Scalar> DirichletBCVector;
    typedef std::vector<std::vector<Scalar> > NeumannBCVector;

public:
    enum
    {
        inside = 999, outside = -999
    };

    //! Constructs a InteractionVolumeInfo object
    /**
     */
    InteractionVolume() :
        stored_(false), normalTimesPerm_(FieldVectorVector(FieldVector(0.0))), nu_(FieldVectorVector(FieldVector(0.0))), normalByFace_(FieldVectorVector(FieldVector(0.0))), faceArea_(FieldVector(0.0)), dF_(0.0), boundaryType_(inside),  isDirichletSat_(false), indexOnElement_(IndexVector(0.0)), elements_(2*dim), dirichletBC_(0), dirichletSatBC_(0), neumannBC_(0, std::vector<Scalar>(dim,0.0))
    {
    }

    void setStored()
    {
        stored_ = true;
    }

    bool isStored() const
    {
        return stored_;
    }

    void setSubVolumeElement(ElementPointer& pointer, int subVolumeIdx)
    {
        elements_[subVolumeIdx].push_back(pointer);
    }

    void setDF(Scalar dF, int subVolumeIdx)
    {
        dF_[subVolumeIdx] = dF;
    }

    void setNormalTimesPerm(FieldVector& normal, Scalar& faceArea, FieldMatrix& perm, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInside] = 0;

//        std::cout<<"perm = "<<perm<<
//                "\n normal "<<normal<<"\n";

        for (int col = 0; col < dim; col++)
        {
            for (int row = 0; row < dim; row++)
            {
                normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInside][col] += normal[row] * perm[row][col] * faceArea;
            }
        }
        normalByFace_[subVolumeIdx][subVolumeFaceIdxInInside] = normal;
        normalByFace_[subVolumeIdx][subVolumeFaceIdxInInside] /= faceArea;
        faceArea_[subVolumeIdx][subVolumeFaceIdxInInside] = faceArea;
    }

    void setNu(FieldVector& nu, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        nu_[subVolumeIdx][subVolumeFaceIdxInInside] = nu;
    }

    void setBoundary(int boundaryType, int subVolumeFaceIdx)
    {
//        std::cout<<"old BCType = "<<boundaryType_[subVolumeFaceIdx]<<"\n";
//        std::cout<<"new BCType = "<<boundaryType<<"\n";
        boundaryType_[subVolumeFaceIdx] = boundaryType;
    }

    void setSatBoundDirichlet(int boundaryType, int subVolumeFaceIdx)
    {
        if (boundaryType == BoundaryConditions::dirichlet)
        isDirichletSat_[subVolumeFaceIdx] = true;
    }

    void setBoundaryCondition(Scalar condition, int subVolumeFaceIdx)
    {
        if (!dirichletBC_.size())
        {
            dirichletBC_.resize(2 * dim, 0.0);
        }
        dirichletBC_[subVolumeFaceIdx] = condition;
    }

    void setBoundaryCondition(std::vector<Scalar>& condition, int subVolumeFaceIdx)
    {
        if (!neumannBC_.size())
        {
            neumannBC_.resize(2 * dim, std::vector<Scalar>(2, 0.0));
        }
        neumannBC_[subVolumeFaceIdx].swap(condition);
    }

    void setDirichletSat(Scalar condition, int subVolumeFaceIdx)
    {
        if (!dirichletSatBC_.size())
        {
            dirichletSatBC_.resize(2 * dim, 0.0);
        }
        dirichletSatBC_[subVolumeFaceIdx] = condition;
    }

    void setIndexOnElement(int indexInInside, int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        indexOnElement_[subVolumeIdx][subVolumeFaceIdxInInside] = indexInInside;
    }

    int getIndexOnElement(int subVolumeIdx, int subVolumeFaceIdx)
    {
        return indexOnElement_[subVolumeIdx][subVolumeFaceIdx];
    }

    ElementPointer& getSubVolumeElement(int subVolumeIdx)
    {
        return elements_[subVolumeIdx][0];
    }

    int getBoundaryType(int subVolumeFaceIdx)
    {
        return boundaryType_[subVolumeFaceIdx];
    }

    bool isDirichletSatBound(int subVolumeFaceIdx)
    {
        return isDirichletSat_[subVolumeFaceIdx];
    }

    template<BoundaryConditions::Flags bcType> typename ReturnBoundaryCondition<DirichletBCVector, NeumannBCVector, bcType>::BCType& getBoundaryCondition(
            int subVolumeFaceIdx)
    {
        return ReturnBoundaryCondition<DirichletBCVector, NeumannBCVector, bcType>::returnBoundaryCondition(dirichletBC_, neumannBC_)[subVolumeFaceIdx];
    }

    Scalar& getDirichletSat(int subVolumeFaceIdx)
    {
        return dirichletSatBC_[subVolumeFaceIdx];
    }

    FieldVector& getUnitOuterNormalByFace(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return normalByFace_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    Scalar& getFaceArea(int subVolumeIdx, int subVolumeFaceIdxInInside)
    {
        return faceArea_[subVolumeIdx][subVolumeFaceIdxInInside];
    }

    //returns n^TKK_rnu
    Scalar getNTKNu(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInsideN] * nu_[subVolumeIdx][subVolumeFaceIdxInInsideNu];
    }

    //returns n^TKK_rnu
    Scalar getNTKKrNu(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        FieldVector nTKKr(normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInsideN]);
        nTKKr *= relPerm;

        return nTKKr * nu_[subVolumeIdx][subVolumeFaceIdxInInsideNu];
    }

    //returns n^TKK_rnu/dF
    Scalar getNTKNu_by_dF(int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        return normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInsideN] * nu_[subVolumeIdx][subVolumeFaceIdxInInsideNu]
                / dF_[subVolumeIdx];
    }

    //returns n^TKK_rnu/dF
    Scalar getNTKKrNu_by_dF(Scalar& relPerm, int subVolumeIdx, int subVolumeFaceIdxInInsideN, int subVolumeFaceIdxInInsideNu) const
    {
        FieldVector nTKKr(normalTimesPerm_[subVolumeIdx][subVolumeFaceIdxInInsideN]);
        nTKKr *= relPerm;

        return nTKKr * nu_[subVolumeIdx][subVolumeFaceIdxInInsideNu] / dF_[subVolumeIdx];
    }

private:
    bool stored_;
    Dune::FieldVector<FieldVectorVector, 2*dim> normalTimesPerm_;
    Dune::FieldVector<FieldVectorVector, 2*dim> nu_;
    Dune::FieldVector<FieldVectorVector, 2*dim> normalByFace_;
    Dune::FieldVector<Dune::FieldVector<Scalar, dim>, 2*dim> faceArea_;
    Dune::FieldVector<Scalar, 2*dim> dF_;
    Dune::FieldVector<int, 2*dim> boundaryType_;
    Dune::FieldVector<bool, 2*dim> isDirichletSat_;
    Dune::FieldVector<IndexVector, 2*dim> indexOnElement_;
    std::vector<std::vector<ElementPointer> > elements_;
    DirichletBCVector dirichletBC_;
    DirichletBCVector dirichletSatBC_;
    NeumannBCVector neumannBC_;
};
}
#endif
