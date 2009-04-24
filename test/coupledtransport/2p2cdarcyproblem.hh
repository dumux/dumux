#ifndef DUNE_TWOPTWOCDARCYPROBLEM_HH
#define DUNE_TWOPTWOCDARCYPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/phaseproperties/phaseproperties_waterair.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/2p2c/2p2cproblem.hh>

/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
//! class that defines the parameters of an air injection under a low permeable layer
/*! Problem definition of an air injection under a low permeable layer. Air enters the domain
 * at the right boundary and migrates upwards.
 * Problem was set up using the rect2d.dgf grid.
 *
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar>
class TwoPTwoCDarcyProblem : public TwoPTwoCProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,numEq> SolutionVector;

public:
    enum {pWIdx = 0, switchIdx = 1}; // phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // phase state


    /////////////////////////////
    // TYPE of the boundaries
    /////////////////////////////
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const GlobalPosition& globalPos, const Element& element,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const LocalPosition& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (globalPos[1] < eps_)
            values = BoundaryConditions::dirichlet;

        return values;
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    virtual SolutionVector dirichlet(const GlobalPosition& globalPos, const Element& element,
                              const IntersectionIterator& intersectionIt,
                              const LocalPosition& localPos) const
    {
        return initial(globalPos, element, localPos);
    }

    /////////////////////////////
    // NEUMANN boundaries
    /////////////////////////////
    virtual SolutionVector neumann(const GlobalPosition& globalPos, const Element& element,
                              const IntersectionIterator& intersectionIt,
                              const LocalPosition& localPos) const
    {
        SolutionVector values(0);

        return values;
    }

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    virtual SolutionVector initial (const GlobalPosition& globalPos, const Element& element,
                                    const LocalPosition& localPos) const
    {

        SolutionVector values;
        Scalar densityW_ = 1000.0;

        values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - globalPos[1]);
        values[switchIdx] = 1e-6;

        return values;
    }


    int initialPhaseState (const GlobalPosition& globalPos, const Element& element,
                           const LocalPosition& localPos) const
    {

        enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states

        return gasPhase;
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    virtual SolutionVector q (const GlobalPosition& globalPos, const Element& element,
                              const LocalPosition& localPos) const
    {
        SolutionVector values(0);

        return values;
    }

    //////////////////////////////

    virtual FieldVector<Scalar,dim> gravity () const
    {
        return gravity_;
    }

    double depthBOR () const
    {
        return depthBOR_;
    }

    TwoPTwoCDarcyProblem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, Scalar>& soil, Scalar depthBOR,
                         TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>),
                         MultiComp& multicomp = *(new CWaterAir))
        : TwoPTwoCProblem<Grid,Scalar>(liq, gas, soil, multicomp, law),
          depthBOR_(depthBOR)
    {
        gravity_[0] = 0;
        gravity_[1] = 0;//-9.81;

        eps_ = 1e-10;
    }

private:
    //      Scalar densityW_, densityN_;
    FieldVector<Scalar,dim> gravity_;
    Scalar depthBOR_;
    Scalar eps_;
};

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/** \todo Please doc me! */

template<class Grid, class Scalar>
class TwoPTwoCDarcySoil: public Matrix2p<Grid,Scalar>
{
    enum {dim=Grid::dimension, numEq=1};

public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return permeability_;
    }
    virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.3;
    }

    virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
    {
        return 0.2;
    }

    virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
    {
        return 0.0;
    }

    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
    {
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 2.; // lambda
        param[1] = 1e4; // entry-pressures

        return param;
    }

    virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid,Scalar>::brooks_corey;
    }

    TwoPTwoCDarcySoil():Matrix2p<Grid,Scalar>()
    {
        permeability_ = 0.;
        for(int i = 0; i < dim; i++)
            permeability_[i][i] = 1e-2;
    }

    ~TwoPTwoCDarcySoil()
    {}

private:
    FieldMatrix permeability_;
};


} //end namespace
#endif
