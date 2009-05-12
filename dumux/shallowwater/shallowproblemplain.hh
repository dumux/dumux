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
    Scalar gravity_;

public:

    BoundaryConditions::Flags bctypeConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {

        if (faceGlobalPos[0] < eps_ )
        {
            return Dune::BoundaryConditions::neumann;
        }
        else if (faceGlobalPos[0] > upperRight_[0]- eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        else
        {
            // DUNE_THROW(NotImplemented,"Problem with BC!");
        }

    }

    BoundaryConditions::Flags bctypeMomentum(
            const GlobalPosition& faceGlobalPos, const Element& element,
            const LocalPosition& localPos) const
    {

        if (faceGlobalPos[0] < eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        else if (faceGlobalPos[0] > upperRight_[0]- eps_)
        {
            return Dune::BoundaryConditions::neumann;
        }

        else
        {
            //DUNE_THROW(NotImplemented, "Problem with BC!");
        }
    }

    Scalar dirichletConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return 0.33;

    }

    VelType neumannConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos, VelType flux = 0) const
    {
        VelType hu(0.18);
        
        return hu;
    }

    VelType dirichletMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        // return momentum (neumann boundary condition for conti part = hu)
         VelType hu(0.18);  

        return hu;
    }

    VelType neumannMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos, VelType flux = 0) const
    {
        VelType boundaryFlux_(0);
        
        //free flow condition
            boundaryFlux_ = flux;
            return boundaryFlux_;      
    }
          
        
    Scalar setInitialWaterDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        int waterDepthDistribution; // type of waterdepth distribution, options see below

        waterDepthDistribution = 4;
        Scalar initWaterDepth = 0;

        switch (waterDepthDistribution)
        {
        case 1: //costant

            initWaterDepth = 0.0;
            return initWaterDepth;
            break;

        case 2: // dam 
            if (globalPos[0]<= 100 /*&& globalPos[dim-1]< 20*/)
            {
                initWaterDepth = 1;
            }
            else
            {
                initWaterDepth = 0;
            }
            break;

        case 3: // pulse
            if (globalPos[0]> 11 && globalPos[dim-1]> 11 && globalPos[0]< 14
                    && globalPos[dim-1]< 14)
            {
                initWaterDepth = 0.5;
            }
            else
            {
                initWaterDepth = 0.2;
            }
            break;
        case 4: //constant water level over irregular topography          

            initWaterDepth = 0.33-this->surface.evalBottomElevation(globalPos);
            break;
        }

        return initWaterDepth;
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
    
    Scalar defineGravity() const
    {
        return gravity_;
    }
    
    ShallowProblemPlain(VC& variableobject, Surface& surfaceobject,
            FieldVector<Scalar,dim>& L, FieldVector<Scalar,dim>& H) :
        ShallowProblemBase<Grid, Scalar, VC>(variableobject, surfaceobject),
                lowerLeft_(L), upperRight_(H), eps_(1e-8), gravity_(9.81)
    {
    }

};

}
#endif
