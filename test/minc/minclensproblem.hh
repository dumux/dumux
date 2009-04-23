
#ifndef DUNE_MINCLENSPROBLEM_HH
#define DUNE_MINCLENSPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/minc/mincproblem.hh>
#include "minc_soilproperties.hh"
/**
 * @file
 * @brief  Base class for defining an instance of the Minc problem
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
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar, int numEq>
class MincLensProblem : public MincProblem<Grid, Scalar, numEq> {
    typedef typename Grid::ctype DT;
    enum {dim=Grid::dimension};
    enum {nCont = numEq/2};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

public:
    // Constructor
    MincLensProblem(Fluid& liq1, Fluid& liq2, MincLensSoil<Grid, Scalar>& soil,
                    const FieldVector<DT,dim>& outerLowerLeft = 0., const FieldVector<DT,dim>& outerUpperRight = 0,
                    const FieldVector<DT,dim>& innerLowerLeft = 0., const FieldVector<DT,dim>& innerUpperRight = 0,
                    TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>))
        : MincProblem<Grid,Scalar, numEq>(liq1, liq2, soil, law),
          outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
          innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
          eps_(1e-8*outerUpperRight[0]),
          densityW_(liq1.density()), densityN_(liq2.density())
    {
        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }


    enum {F = 0, M = 1};
    enum {pWFIdx = 0, sNFIdx = 1, pWMIdx = 2, sNMIdx = 3};
    enum {swrIdx = 0, snrIdx = 1, alphaIdx = 2, nIdx = 3};


    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<DT,dim>& x, const Entity& e,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<DT,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::neumann);

        if (x[0] > outerUpperRight_[0] - eps_) {
            //std::cout << "Dirichlet: " << x << std::endl;
            values = BoundaryConditions::dirichlet;
        }
        //        if (values[0] == BoundaryConditions::dirichlet)
        //            std::cout << "Dirichlet: " << x[0] << ", " << x[1] << std::endl;
        //        else
        //            std::cout << "Neumann: " << x[0] << ", " << x[1] << std::endl;

        return values;
    }

    virtual FieldVector<Scalar,numEq> g (const FieldVector<DT,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        for (int nC=0; nC < numEq-1; nC+=2){

            if (x[0] > outerUpperRight_[0] - eps_) {
                Scalar a = -1;
                Scalar b = outerUpperRight_[1];
                // Fracture domain
                if (nC<2){
                    values[pWFIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
                    values[sNFIdx] = outerSnrFracture_+0.0;
                }
                // Matrix domain
                else {
                    values[nC] = values[pWFIdx]*1.0;            //pWM = pWF
                    values[nC+1] = outerSnrMatrix_;
                }
            }
        }
        return values;
    }

    virtual FieldVector<Scalar,numEq> J (const FieldVector<DT,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        if (x[0] < outerLowerLeft_[0] + eps_) {
            values[sNFIdx] = -0.5;
            values[pWFIdx] = -0.0;
        }

        return values;
    }

    virtual FieldVector<Scalar,numEq> initial (const FieldVector<DT,dim>& x, const Entity& e,
                                               const FieldVector<DT,dim>& xi) const
    {

        FieldVector<Scalar,numEq> values;
        values = 0;

        for (int nC=0; nC < numEq-1; nC+=2){
            if (nC==0){
                values[nC] = -densityW_*gravity_[1]*(outerUpperRight_[1] - x[1])+0;
                values[nC+1] = outerSnrFracture_+0.0;
            }
            else
            {
                values[nC] = (-densityW_*gravity_[1]*(outerUpperRight_[1] - x[1])+0)*0.9;
                values[nC+1] = outerSnrMatrix_+0.0;
            }
        }
        return values;
    }

    // function returning SOURCE/SINK terms
    virtual FieldVector<Scalar,numEq> q (const FieldVector<DT,dim>& x, const Entity& element,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<Scalar,dim> gravity () const
    {
        return gravity_;
    }


private:
    FieldVector<DT,dim> outerLowerLeft_;
    FieldVector<DT,dim> outerUpperRight_;
    FieldVector<DT,dim> innerLowerLeft_;
    FieldVector<DT,dim> innerUpperRight_;
    DT eps_;
    Scalar densityW_, densityN_;
    FieldVector<DT,dim> gravity_;
    Scalar volumefraction;

    //        RT outerSwr_, outerSnr_, innerSwr_, innerSnr_;
    Scalar outerSwrFracture_, outerSnrFracture_, innerSwrFracture_, innerSnrFracture_;
    Scalar outerSwrMatrix_, outerSnrMatrix_, innerSwrMatrix_, innerSnrMatrix_;
};

}
#endif
