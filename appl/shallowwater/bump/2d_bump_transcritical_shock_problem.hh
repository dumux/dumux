// $Id: shallowproblemplain.hh 1501 2009-03-26 15:52:22Z anneb $

#ifndef DUNE_2D_BUMP_TRANSRITICAL_PROBLEM_HH
#define DUNE_2D_BUMP_TRANSRITICAL_PROBLEM_HH

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
            const Element& element, const LocalPosition& localPos,
            const Scalar& froudeNumber) const
    {
        if (faceGlobalPos[0] > upperRight_[0]-eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        else
        {
            return Dune::BoundaryConditions::neumann;
        }
    }

    BoundaryConditions::Flags bctypeMomentum(
            const GlobalPosition& faceGlobalPos, const Element& element,
            const LocalPosition& localPos, const Scalar& froudeNumber) const
    {
        if (faceGlobalPos[0]< eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        else
        {
            return Dune::BoundaryConditions::neumann;
        }
    }

    Scalar dirichletConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return 0.33;

    }

    Scalar neumannConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            Scalar flux = 0) const
    {
        Scalar contiFlux = 0;

        if (faceGlobalPos[0]< eps_)
        {
            contiFlux = 0.18;
        }
        else if (faceGlobalPos[0]> upperRight_[0]-eps_)
        {
            contiFlux = flux;
        }
        else if (faceGlobalPos[1]< eps_ || faceGlobalPos[1] < upperRight_[1]
                -eps_)
        {
            contiFlux = 0;
        }
        return contiFlux;
    }

    VelType dirichletMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        VelType momentum_(0);

        momentum_[0] = 0.18;
        momentum_[1] = 0;

        return momentum_;
    }

    VelType neumannMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            const Scalar& waterDepth, VelType flux = VelType(0.0)) const
    {
        VelType momentumFlux_(0);

        if (faceGlobalPos[1]< eps_)
        {
            momentumFlux_[0] = 0;
            momentumFlux_[1] = 0.5*9.81*waterDepth*waterDepth*(-1);
        }
        if (faceGlobalPos[1] > upperRight_[1]-eps_)
        {
            momentumFlux_[0] = 0;
            momentumFlux_[1] = 0.5*9.81*waterDepth*waterDepth;
        }
        else
        {
            momentumFlux_ = flux;
        }

        return momentumFlux_;
    }

    Scalar setInitialWaterDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        Scalar initWaterDepth = 0;
        initWaterDepth = 0.33 - this->surface.evalBottomElevation(globalPos);

        return initWaterDepth;
    }

    VelType setInitialVelocity(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        VelType initialVelocity(0);
        initialVelocity[0]= 0;
        initialVelocity[1]= 0;

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
                lowerLeft_(L), upperRight_(H), eps_(1e-12)
    {
    }

};

}
#endif
