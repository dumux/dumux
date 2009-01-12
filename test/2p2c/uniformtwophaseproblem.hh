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
   *    - RT    type used for return values
   */
  template<class Grid, class RT>
  class UniformTwoPhaseProblem : public TwoPhaseProblem<Grid, RT> {
    typedef typename Grid::ctype Scalar;
    enum {n=Grid::dimension, m=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
    virtual const FieldMatrix<Scalar,n,n>& K (const FieldVector<Scalar,n>& x, const Element& e,
                    const FieldVector<Scalar,n>& xi)
    {
        return permloc;
    }

    virtual FieldVector<RT,m> q (const FieldVector<Scalar,n>& x, const Element& e,
                    const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<Scalar,n>& x, const Element& e,
                    const IntersectionIterator& intersectionIt,
                       const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::neumann);

        if (x[0] > 600-1E-6)// || x[0] < 1e-6)
        //if (x[0] < 1E-6)
            values = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual FieldVector<RT,m> g (const FieldVector<Scalar,n>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        if (x[0] < 1e-6) {
            values[0] = 1e6;
            values[1] = 1;
        }
        else {
            values[0] = 0;
            values[1] = 0;
        }

        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<Scalar,n>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        if (x[0] < 1e-6)
            values[0] = 1;

        return values;
    }

    virtual FieldVector<RT,m> initial (const FieldVector<Scalar,n>& x, const Element& e,
                  const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        values[0] = 1e6 - 1.0/600.0*1e6*x[0];
        //values[1] = 1 - 1.0/600.0*x[0];
        if (x[0] < 1e-6) {
            values[1] = 1e3;
        }

        return values;
    }

    double porosity (const FieldVector<Scalar,n>& x, const Element& e,
              const FieldVector<Scalar,n>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<RT,n> gravity () const
    {
        FieldVector<RT,n> values(0);

        return values;
    }

    virtual FieldVector<RT,4> materialLawParameters (const FieldVector<Scalar,n>& x, const Element& e,
              const FieldVector<Scalar,n>& xi) const
    {
        FieldVector<RT,4> values(0);

        return values;
    }

    UniformTwoPhaseProblem(TwoPhaseRelations& law = *(new LinearLaw))
    : TwoPhaseProblem<Grid, RT>(law)
    {
        permloc = 0;

        for (int i = 0; i < n; i++)
            permloc[i][i] = 1.0;
    }

    private:
        Dune::FieldMatrix<Scalar,n,n> permloc;
  };

}
#endif
