// $Id: simpleproblem.hh 1511 2009-03-31 09:23:02Z markus $

#ifndef DUNE_DYNAMICPCPROBLEM_HH
#define DUNE_DYNAMICPCPROBLEM_HH

#include <dumux/material/matrixproperties.hh>
#include "dumux/transport/transportproblem.hh"
#include <dumux/material/property_baseclasses.hh>

namespace Dune
{
template<class Grid, class Scalar, class VC>
class DynamicPcProblem : public TransportProblem<Grid, Scalar, VC> {
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;

public:
    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                         const LocalPosition& localPos) const
    {
        if (globalPos[0] > upperRight_[0]-1E-8 || globalPos[0] < lowerLeft_[0]+1e-8)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
                         const LocalPosition& localPos) const
    {
        if (globalPos[0] < lowerLeft_[0]+1e-8)
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

    DynamicPcProblem(VC& variableobj,Fluid& wettingPhase, Fluid& nonwettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid,Scalar>& materialLaw  = *(new TwoPhaseRelations<Grid,Scalar>), GlobalPosition& left = 0, GlobalPosition& right = 1)
        : TransportProblem<Grid, Scalar, VC>(variableobj, wettingPhase, nonwettingPhase, soil,materialLaw), lowerLeft_(left), upperRight_(right)
    {    }
};

}
#endif
