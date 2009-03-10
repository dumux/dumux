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
    BoundaryConditions::Flags bctype(const GlobalPosition& globalPos,
                                     const Element& element, const LocalPosition& localPos) const
    {

        if (globalPos[0] < eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        else
        {
            return Dune::BoundaryConditions::neumann;
        }
    }

    Scalar dirichlet(const GlobalPosition& globalPos, const Element& element,
                     const LocalPosition& localPos) const
    {
        return 0;
    }

    SystemType neumann(const GlobalPosition& globalPos, const Element& element,
                       const LocalPosition& localPos) const
    {
        SystemType boundaryFlux;

        if (globalPos[0]== upperRight_[0] - eps_)
        {
            boundaryFlux[0]=0;
            boundaryFlux[1]=0;
            boundaryFlux[2]=0;

        }
        return boundaryFlux;

    }

    Scalar setInitWDepth(const GlobalPosition& globalPos,
                         const Element& element, const LocalPosition& localPos) const
    {
        return 0.00;
    }

    VelType setInitVel(const GlobalPosition& globalPos, const Element& element,
                       const LocalPosition& localPos) const
    {
        VelType initVel_(0);
        return initVel_;

    }

    Scalar setSource(const GlobalPosition& globalPos, const Element& element,
                     const LocalPosition& localPos) const
    {
        return 125E-7; //der Umrechnungsfaktor von l/(s ha) auf m/s ist 1E-7, Standardwert für Stuttgart 125 l/(s ha)
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
