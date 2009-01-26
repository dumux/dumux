
#ifndef DUNE_LENSPROBLEM_PNSW_HH
#define DUNE_LENSPROBLEM_PNSW_HH

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
  class LensProblem : public TwoPhaseProblem<Grid, Scalar> {
    typedef typename Grid::ctype DT;
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
    // Constructor
    LensProblem(Fluid& liq1, Fluid& liq2, Matrix2p<Grid, Scalar>& soil,
            const FieldVector<DT,dim>& outerLowerLeft = 0., const FieldVector<DT,dim>& outerUpperRight = 0,
            const FieldVector<DT,dim>& innerLowerLeft = 0., const FieldVector<DT,dim>& innerUpperRight = 0,
            TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>))
    : TwoPhaseProblem<Grid,Scalar>(liq1, liq2, soil, law),
        outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
        innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
        eps_(1e-8*outerUpperRight[0]),
        densityW_(liq1.density()), densityN_(liq2.density())
    {
        height_ = outerUpperRight[1] - outerLowerLeft[1];
        width_ = outerUpperRight[0] - outerLowerLeft[0];

        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }

    // indices for the two phases
    enum {pNIdx = 0, sWIdx = 1};

    // function returning the BOUNDARY CONDITION TYPE depending on the position
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<DT,dim>& x, const Entity& element,
                    const IntersectionIterator& intersectionIt,
                       const FieldVector<DT,dim>& xi) const
    {
        // boundary condition type is set to Neumann
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (x[0] < outerLowerLeft_[0] + eps_ || x[0] > outerUpperRight_[0] - eps_)
        {
            values = BoundaryConditions::dirichlet;
        }

        return values;
    }

    // definition of DIRICHLET boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> g (const FieldVector<DT,dim>& x, const Entity& element,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        if (x[0] < outerLowerLeft_[0] + eps_)
        {
            Scalar a = -(1 + 0.5/height_);
            Scalar b = -a*outerUpperRight_[1];
            values[sWIdx] = 1;
            values[pNIdx] = -densityW_*gravity_[1]*(a*x[1] + b);    // pN = pW for sW = 1 and pC = 0
        }
        else {
            Scalar a = -1;
            Scalar b = outerUpperRight_[1];
            values[sWIdx] = 1;
            values[pNIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
        }

        return values;
    }

    // definition of NEUMANN boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> J (const FieldVector<DT,dim>& x, const Entity& element,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        Scalar lambda = (outerUpperRight_[0] - x[0])/width_;
        if (lambda > 0.5 && lambda < 2.0/3.0 && x[1] > outerUpperRight_[1] - eps_) {
            values[pNIdx] = -0.04;
        }

        return values;
    }

    // definition of INITIAL VALUES for pressure and saturation
    virtual FieldVector<Scalar,numEq> initial (const FieldVector<DT,dim>& x, const Entity& element,
                  const FieldVector<DT,dim>& xi) const
    {

        FieldVector<Scalar,numEq> values;

        values[pNIdx] = -densityW_*gravity_[1]*(height_ - x[1]);

        if (x[0] < outerLowerLeft_[0] + eps_) {
            Scalar a = -(1 + 0.5/height_);
            Scalar b = -a*outerUpperRight_[1];
            values[pNIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
        }

        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            values[sWIdx] = 1;//innerSnr_;
        else
            values[sWIdx] = 1;//outerSnr_;

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
        DT width_, height_;
        DT eps_;
        Scalar densityW_, densityN_;
        FieldVector<DT,dim> gravity_;
  };

} // end namespace
#endif
