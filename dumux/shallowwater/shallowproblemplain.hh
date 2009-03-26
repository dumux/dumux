// $Id$

#ifndef DUNE_SHALLOWPROBLEMPLAIN_HH
#define DUNE_SHALLOWPROBLEMPLAIN_HH

#include "dumux/shallowwater/shallowproblembase.hh"
#include"dumux/shallowwater/shallowvariableclass.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem in shallow water
template<class Grid, class Scalar, class VC> class ShallowProblemPlain :
    public ShallowProblemBase<Grid, Scalar,VC>
{
    enum
    {   dim=Grid::dimension, m=1, blocksize=2*Grid::dimension};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar,dim+1> SystemType;
    typedef Dune::SolidSurfaceBase<Grid,Scalar> Surface;

private:
    FieldVector<Scalar,dim> lowerLeft_;
    FieldVector<Scalar,dim> upperRight_;
    Scalar eps_;

public:
    BoundaryConditions::Flags bctype(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {

        int boundarySetting;
        
        boundarySetting =2;

        switch (boundarySetting)
        {
        case 1:

            if (faceGlobalPos[0] < eps_)
            {
                return Dune::BoundaryConditions::dirichlet;
            }

            if (faceGlobalPos[0] >= upperRight_[0] - eps_)
            {
                return Dune::BoundaryConditions::neumann;
            }
            break;

        case 2:

            if (faceGlobalPos[0]< eps_ || faceGlobalPos[0] >= upperRight_[0]
                    - eps_)
            {
                return Dune::BoundaryConditions::neumann;
            }
            break;
        }
    }

    Scalar dirichletWaterDepth(const GlobalPosition& faceGlobalPos) const
    {
        return 0;
    }

    Scalar dirichletVelocity(const GlobalPosition& faceGlobalPos) const
    {
        return 0.0;
    }

    SystemType neumannWaterDepth(const GlobalPosition& faceGlobalPos) const
    {
        SystemType boundaryFlux(0);
        return boundaryFlux;
    }

    SystemType neumannVelocity(const GlobalPosition& faceGlobalPos) const
    {
        SystemType boundaryFlux(0);
        return boundaryFlux;
    }

    SystemType neumannFlux(const GlobalPosition& faceGlobalPos,
            SystemType& helpFlux) const
    {
        if (faceGlobalPos[0] <= lowerLeft_+eps_)
        {
            SystemType boundaryFlux(0);
            return helpFlux;
        }

        if (faceGlobalPos[0]>= upperRight_[0] - eps_)
        {
            SystemType boundaryFlux(0);
            return helpFlux;
        }

    }

    Scalar setInitialWaterDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        int waterDepthDistribution; // 1 = constant, 2 = high at left side, 3 = high in the middle

        waterDepthDistribution = 2;

        switch (waterDepthDistribution)
        {
        case 1:
            return 0.0001;
            break;
        case 2:
            if (globalPos[0]< 40)
            {
                return 0.5;
            }
            else
            {
                return 0;
            }
            break;
        case 3:
            if (globalPos[0] > 45 && globalPos[0]<55 )
            {
                return 0.5;
            }
            else
            {
                return 0.2;
            }
            break;
        }

    }

    VelType setInitialVelocity(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {

        VelType initialVelocity(0);

        return initialVelocity;

    }

    Scalar setSource(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0.0000; //der Umrechnungsfaktor von l/(s ha) auf m/s ist 1E-7, Standardwert für Stuttgart 125 l/(s ha)
    }

    ShallowProblemPlain(VC& variableobject, Surface& surfaceobject,
            FieldVector<Scalar,dim>& L, FieldVector<Scalar,dim>& H) :
        ShallowProblemBase<Grid, Scalar, VC>(variableobject, surfaceobject),
                lowerLeft_(L), upperRight_(H), eps_(1e-8)
    {
    }

};

}
#endif
