// $Id: lensproblem.hh 566 2008-09-11 11:38:31Z bernd $

#ifndef DUNE_DMTPROBLEM_HH
#define DUNE_DMTPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/property_baseclasses.hh>
#include<dumux/twophase/twophaseproblem.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch, Klaus Mosthaf
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
 *    - Grid    a DUNE grid type
 *    - Scalar  type used for return values
 */
template<class Grid, class Scalar>
class DMTProblem : public TwoPhaseProblem<Grid, Scalar> {
    typedef typename Grid::ctype DT;
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

public:
    // Constructor
    DMTProblem(Fluid& liq1, Fluid& liq2, Matrix2p<Grid, Scalar>& soil,
               TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>))
        : TwoPhaseProblem<Grid,Scalar>(liq1, liq2, soil, law)
    {
        gravity_[0] = 0;
        gravity_[1] = -9.81;

        pAtm_ = 1.013e5;
    }

    // indices for the two phases
    enum {pWIdx = 0, sNIdx = 1};

    // function returning the BOUNDARY CONDITION TYPE depending on the position
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<DT,dim>& x, const Entity& element,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<DT,dim>& xi) const
    {
        // boundary condition type is set to neumann
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        switch (intersectionIt->boundaryId()) {
        case 5:
            values = BoundaryConditions::dirichlet;
            break;
        case 3:
        case 4:
            values = BoundaryConditions::neumann;
            break;
        case 1:
            values = BoundaryConditions::neumann;
            break;
        case 2:
            values = BoundaryConditions::neumann;
            break;
        }

        return values;
    }

    // definition of DIRICHLET boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> g (const FieldVector<DT,dim>& x, const Entity& element,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        switch (intersectionIt->boundaryId()) {
        case 5:
            values[pWIdx] = pAtm_;
            values[sNIdx] = 0.1;
            break;
        }

        return values;
    }

    // definition of NEUMANN boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> J (const FieldVector<DT,dim>& x, const Entity& element,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    // definition of INITIAL VALUES for pressure and saturation
    virtual FieldVector<Scalar,numEq> initial (const FieldVector<DT,dim>& x, const Entity& element,
                                               const FieldVector<DT,dim>& xi) const
    {

        FieldVector<Scalar,numEq> values;

        double rhoG = 0.68;
        values[pWIdx] = rhoG*gravity_[1]*x[1] + pAtm_;
        values[sNIdx] = 0.1;

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
    FieldVector<DT,dim> gravity_;
    DT pAtm_;
};

} // end namespace
#endif
