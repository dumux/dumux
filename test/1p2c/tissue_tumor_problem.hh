// $Id$
#ifndef DUNE_TISSUE_TUMOR_PROBLEM_HH
#define DUNE_TISSUE_TUMOR_PROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/phaseproperties/phaseproperties1p.hh>
#include<dumux/1p2c/1p2cproblem.hh>
#include"tissue_soilproperties.hh"

/**
 * @file
 * @brief  Base class for defining an instance of the OnePhase problem
 * @author Karin Erbertseder
 */

namespace Dune
{
//! base class that defines the parameters for the boundary conditions and the source/sink term
/*! Template parameters are:
 *
 *    - Grid   a DUNE grid type
 *    - Scalar type used for return values
 */
template<class Grid, class Scalar>
class TissueTumorProblem : public OnePTwoCProblem<Grid, Scalar> {
    typedef typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator IntersectionIterator;
    enum {dim=Grid::dimension, numEq=2};
    enum {konti = 0, transport = 1};    // Solution vector index

public:

    virtual FieldVector<Scalar,numEq> q (const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                         const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        if(golbalPos[0]>10 && golbalPos[0]<12 && golbalPos[1]>10 && golbalPos[1]<12)
            values[0]=1.5e-6;
        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::dirichlet);

        //        if( golbalPos[0]<1e-1 || golbalPos[0]> 22-1e-1 )
        //            values = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual void dirichletIndex(const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<CoordScalar,dim>& localPos, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual FieldVector<Scalar,numEq>dirichlet(const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        if(golbalPos[0]<1e-1)
        {
            values[0] = -931;
            values[1] = 1.1249e-8;
        }

        else if(golbalPos[0]>22-1e-1)
        {
            values[0] = -1067;      //Dirichlet RB für Druck
            values[1] = 1.1249e-8;       //Dirichlet RB für mol fraction x
        }

        else
        {
            values[0] =-931-((136/22)*golbalPos[0]); //AB für Druck p
            values[1] = 1.1249e-8; //AB für mole fraction x
        }
        return values;
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                         const IntersectionIterator& intersectionIt, const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        //        if(golbalPos[0] < 1.e-1)
        //        {
        //            values[0] = -3.8676e-8;
        //            values[1] = -4.35064e-16;
        //        }
        return values;
    }

    // Initial Conditions for global vector x, element e and local vector xi
    virtual FieldVector<Scalar,numEq> initial (const FieldVector<CoordScalar,dim>& golbalPos, const Element& element,
                                               const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        values[0] =-931-((136/22)*golbalPos[0]); //AB für Druck p

        //        values[0] = -1067;
        //        values[1] = 1.1249e-8;
        if(golbalPos[0]>10 && golbalPos[0]<12 && golbalPos[1]>10 && golbalPos[1]<12)
            values[1]=0;
        else
            values[1]=1.1249e-8;

        return values;
    }



    TissueTumorProblem(Liquid_GL& phase, Matrix2p<Grid, Scalar>& soil)
        : OnePTwoCProblem<Grid, Scalar>(phase, soil)
    {

    }
};

}
#endif
