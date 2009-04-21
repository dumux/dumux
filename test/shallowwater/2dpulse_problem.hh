// $Id: shallowproblemplain.hh 1501 2009-03-26 15:52:22Z anneb $

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

    BoundaryConditions::Flags bctypeConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return Dune::BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeMomentum(
            const GlobalPosition& faceGlobalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return 0;

    }

    VelType neumannConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            VelType flux = VelType(0.0)) const
    {
        VelType hv(0);

        //free flow condition: return what you get from fvshallowwater

      //  hv = flux;

        // no flow condition


        return hv;
    }

    VelType dirichletMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        VelType hv(0);

        //free flow condition: return what you get from fvshallowwater

        //hv = flux;

        // no flow condition


        return hv;
    }

    VelType neumannMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            VelType flux = VelType(0.0)) const
    {
        VelType boundaryFlux_(0);

        //free flow condition
        boundaryFlux_ = flux;
        return boundaryFlux_;
    }

    Scalar setInitialWaterDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        Scalar initWaterDepth = 0;

        if (globalPos[0]> 23 && globalPos[dim-1]> 23 && globalPos[0]< 27
                && globalPos[dim-1]< 27)
        {
            initWaterDepth = 0.5;
        }
        else
        {
            initWaterDepth = 0.2;
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

    ShallowProblemPlain(VC& variableobject, Surface& surfaceobject,
            FieldVector<Scalar,dim>& L, FieldVector<Scalar,dim>& H) :
        ShallowProblemBase<Grid, Scalar, VC>(variableobject, surfaceobject),
                lowerLeft_(L), upperRight_(H), eps_(1e-8)
    {
    }

};

}
#endif
