// $Id: noflow_pulseproblem.hh 1816 2009-05-15 13:55:10Z anneb $

#ifndef DUNE_NOFLOW_PULSE_PROBLEM_HH
#define DUNE_NOFLOW_PULSE_PROBLEM_HH

#include "dumux/shallowwater/shallowproblem.hh"
#include"dumux/shallowwater/shallowvariableclass.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem in shallow water
template<class Grid, class Scalar, class VC> class NoflowPulseProblem :
    public ShallowProblem<Grid, Scalar,VC>
{
    enum
    {   dim=Grid::dimension, m=1, blocksize=2*Grid::dimension};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar,dim+1> SystemType;
    typedef Dune::BottomProperties<Grid,Scalar> Surface;

private:
    FieldVector<Scalar,dim> lowerLeft_;
    FieldVector<Scalar,dim> upperRight_;
    Scalar eps_;

public:

    BoundaryConditions::Flags bctypeConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos, const Scalar& froudeNumbe) const
    {
        return Dune::BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeMomentum(
            const GlobalPosition& faceGlobalPos, const Element& element,
            const LocalPosition& localPos, const Scalar& froudeNumbe) const
    {
        return Dune::BoundaryConditions::dirichlet;
    }

    Scalar dirichletConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return 0;

    }

    Scalar neumannConti(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            Scalar flux = 0) const
    {
        Scalar contiFlux_ = 0;
        return contiFlux_;
    }

    VelType dirichletMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        VelType momentum(0);

        return momentum;
    }

    VelType neumannMomentum(const GlobalPosition& faceGlobalPos,
            const Element& element, const LocalPosition& localPos,
            const Scalar& waterDepth, VelType flux = VelType(0.0)) const
    {
        VelType momentumFlux_(0);
       
        return momentumFlux_;
    }

    Scalar setInitialWaterDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        Scalar initWaterDepth = 0;

        Scalar help = -(globalPos[0]-20)*(globalPos[0]-20)-(globalPos[1]-20)
                *(globalPos[1]-20);
        help /=10;

        initWaterDepth =0.2 + 1.0 * 0.5 * exp(help);
        
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

    NoflowPulseProblem(VC& variableobject, Surface& surfaceobject,
            FieldVector<Scalar,dim>& L, FieldVector<Scalar,dim>& H) :
        ShallowProblem<Grid, Scalar, VC>(variableobject, surfaceobject),
                lowerLeft_(L), upperRight_(H), eps_(1e-8)
    {
    }

};

}
#endif
