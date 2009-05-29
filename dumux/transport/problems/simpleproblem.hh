// $Id$

#ifndef DUNE_SIMPLEPROBLEM_HH
#define DUNE_SIMPLEPROBLEM_HH

#include <dumux/material/matrixproperties.hh>
#include "dumux/transport/transportproblem.hh"
#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class HomogeneousLinearSoil: public HomogeneousSoil<Grid,Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;

    virtual double porosity(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos) const
    {
        return 1.0;
    }

    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos, const double T) const
    {
        std::vector<double> param(2);
        //linear law parameters
        param[0] = 0; // minimal capillary pressure
        param[1] = 0; // maximal capillary pressure

        return param;
    }

    virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid,Scalar>::linear;
    }

    HomogeneousLinearSoil()
        :HomogeneousSoil<Grid,Scalar>()
    {}
};

//! \ingroup transportProblems
//! @brief example class for a transport problem
template<class GridView, class Scalar, class VC>
class SimpleProblem : public TransportProblem<GridView, Scalar, VC> {
    enum {dim=GridView::dimension, numEq=1};
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

private:
    Scalar left;
    Scalar right;

public:
    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                         const LocalPosition& localPos) const
    {
        if (globalPos[0] > right-1E-8 || globalPos[0] < left+1e-8)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
                         const LocalPosition& localPos) const
    {
        if (globalPos[0] < left+1e-8)
            return 1;
        else
            return 0;
    }

    virtual Scalar neumannSat (const GlobalPosition& globalPos, const Element& element,
                               const LocalPosition& localPos, Scalar helpFactor) const
    {
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                    const LocalPosition& localPos) const
    {
        return 0;
    }

    virtual Scalar neumann(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
                               const LocalPosition& localPos) const
    {
        return neumannSat(globalPos, element, localPos, 0);
    }

    BoundaryConditions::Flags bctype (const GlobalPosition& globalPos, const Element& element,
                const IntersectionIterator& intersectionIt,
                                   const LocalPosition& localPos) const
        {
            return bctypeSat(globalPos, element, localPos);
        }

     SimpleProblem(VC& variableobj, Fluid& wettingPhase, Fluid& nonwettingPhase, Matrix2p<Grid,Scalar>& soil, TwoPhaseRelations<Grid,Scalar>& materialLaw  = *(new TwoPhaseRelations<Grid,Scalar>), GlobalPosition& Left = 0, GlobalPosition& Right = 1)
        : TransportProblem<GridView, Scalar, VC>(variableobj, wettingPhase, nonwettingPhase, soil, materialLaw), left(Left[0]), right(Right[0])
    {    }

    SimpleProblem(VC& variableobj,Fluid& wettingPhase, Fluid& nonwettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid,Scalar>& materialLaw  = *(new TwoPhaseRelations<Grid,Scalar>))
        : TransportProblem<GridView, Scalar, VC>(variableobj, wettingPhase, nonwettingPhase, soil,materialLaw), left(0), right(1)
    {    }
};

}
#endif
