// $Id$

#ifndef DUNE_UNIFORMTWOPHASEPROBLEM_HH
#define DUNE_UNIFORMTWOPHASEPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/linearlaw.hh>
#include<dumux/twophase/twophaseproblem_deprecated.hh>

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
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar>
class UniformTwoPhaseProblem : public TwoPhaseProblem<Grid, Scalar> {
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

public:
    virtual const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                                  const FieldVector<Scalar,dim>& localPos)
    {
        return permloc;
    }

    virtual FieldVector<Scalar,numEq> q (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                         const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::neumann);

        if (globalPos[0] > 600-1E-6)// || globalPos[0] < 1e-6)
            //if (globalPos[0] < 1E-6)
            values = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual FieldVector<Scalar,numEq> g (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        if (globalPos[0] < 1e-6) {
            values[0] = 1e6;
            values[1] = 1;
        }
        else {
            values[0] = 0;
            values[1] = 0;
        }

        return values;
    }

    virtual FieldVector<Scalar,numEq> J (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        if (globalPos[0] < 1e-6)
            values[0] = 1;

        return values;
    }

    virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                               const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        values[0] = 1e6 - 1.0/600.0*1e6*globalPos[0];
        //values[1] = 1 - 1.0/600.0*globalPos[0];
        if (globalPos[0] < 1e-6) {
            values[1] = 1e3;
        }

        return values;
    }

    double porosity (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                     const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> gravity () const
    {
        FieldVector<Scalar,dim> values(0);

        return values;
    }

    virtual FieldVector<Scalar,4> materialLawParameters (const FieldVector<Scalar,dim>& globalPos, const Element& e,
                                                         const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,4> values(0);

        return values;
    }

    UniformTwoPhaseProblem(TwoPhaseRelations& law = *(new DeprecatedLinearLaw))
        : TwoPhaseProblem<Grid, Scalar>(law)
    {
        permloc = 0;

        for (int i = 0; i < dim; i++)
            permloc[i][i] = 1.0;
    }

private:
    Dune::FieldMatrix<Scalar,dim,dim> permloc;
};

}
#endif
