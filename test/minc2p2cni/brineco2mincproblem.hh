// $Id: brineco2problem.hh 609 2008-09-18 17:01:25Z melanie $
#ifndef DUNE_BRINECO2PROBLEM_HH
#define DUNE_BRINECO2PROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

//#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/minc2p2cni/minc2p2cniproblem.hh>
#include "mincco2_soilproperties.hh"


/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch
 */

namespace Dune
{
//! base class that defines the parameters of a diffusion equation
/*! An interface for defining parameters for the stationary diffusion equation
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 *
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - RT    type used for return values
 */
template<class G, class RT, int m>
class BrineCO2Problem : public TwoPTwoCNIProblem<G, RT, m> {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
    enum {dim=G::dimension};
    enum {nCont = m/3};
    //    enum {wPhase = 0, nPhase = 1};
    //    enum {pWIdx = 0, switchIdx = 1, teIdx = 2};


public:

    virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
                                 const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
                                                              const IntersectionIterator& intersectionIt,
                                                              const FieldVector<DT,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::neumann);

        if(x[0] >= 10-1e-2 )
            values = Dune::BoundaryConditions::dirichlet;
        if(x[0] < 1e-2 && x[1] < 5)
            values[2] = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual void dirichletIndex(const FieldVector<DT,dim>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,dim>& xi, FieldVector<int,m>& dirichletIndex) const
    {
        for (int i = 0; i < m; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0);
        for (int nC=0; nC < m; nC+=(m/nCont)){

            if (nC<3){
                values[0]=1.013e5 + (depthBOR_ - x[1])*1045*9.81;
                values[1]=0.0;
                values[2]=283.15 + (depthBOR_ - x[1])*0.03;
            }
            else {
                values[nC]=values[0]*1.0;
                values[nC+1]=0.0;
                values[nC+2]=values[2];
            }
        }

        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt, const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0.0);
        if(x[0] <= 1.e-2 && x[1] <= 5.)
        {
            values[0] = 0.0;//-4.046e-5;
            values[1] = -0.02;
        }
        return values;
    }

    // Initial Conditions for global vector x, element e and local vector xi
    virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
                                       const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0.0);
        for (int nEq=0; nEq < m; nEq+=(m/nCont)){

            if (nEq<3){
                values[0]=1.013e5 + (depthBOR_ - x[1])*1045*9.81;
                values[1]=0.0;
                values[2]=283.15 + (depthBOR_ - x[1])*0.03;
            }
            else {
                values[nEq]=values[0]*1.0;
                values[nEq+1]=0.0;
                values[nEq+2]=values[2];
            }
        }

        return values;
    }

    int initialPhaseState (const FieldVector<DT,dim>& x, const Entity& e,
                           const FieldVector<DT,dim>& xi) const
    {
        enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
        int state;

        state = waterPhase;

        return state;
    }

    FieldVector<RT,dim> gravity () const
    {
        FieldVector<RT,dim> values(0.0);

        values[1] = -9.81;

        return values;
    }


    BrineCO2Problem(Liquid_GL& liq, Gas_GL& gas, MincCO2Soil<G, RT>& soil,
                    TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>),
                    MultiComp& multicomp = *(new CWaterAir), RT depthBOR = 0.0)
        : TwoPTwoCNIProblem<G, RT, m>(liq, gas, soil, multicomp, law)
    {
        depthBOR_ = depthBOR;
    }

private:
    RT depthBOR_;
    //        RT soilDens_, soilHeatCp_, soilLDry_, soilLSw_;
};

}
#endif
