// $Id$

#ifndef DUNE_LENSPROBLEM_HH
#define DUNE_LENSPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>
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
class LensProblem : public TwoPhaseProblem<Grid, Scalar> {
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

public:
    enum {pWIdx = 0, sNIdx = 1};
    enum {swrIdx = 0, snrIdx = 1, alphaIdx = 2, nIdx = 3};

    virtual const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& x, const Entity& e,
                                                  const FieldVector<Scalar,dim>& xi)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerK_;
        else
            return outerK_;
    }

    virtual FieldVector<Scalar,numEq> q (const FieldVector<Scalar,dim>& x, const Entity& e,
                                         const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,dim>& x, const Entity& e,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (x[0] < outerLowerLeft_[0] + eps_ || x[0] > outerUpperRight_[0] - eps_) {
            //std::cout << "Dirichlet: " << x << std::endl;
            values = BoundaryConditions::dirichlet;
        }
        //        if (values[0] == BoundaryConditions::dirichlet)
        //            std::cout << "Dirichlet: " << x[0] << ", " << x[1] << std::endl;
        //        else
        //            std::cout << "Neumann: " << x[0] << ", " << x[1] << std::endl;

        return values;
    }

    virtual FieldVector<Scalar,numEq>dirichlet(const FieldVector<Scalar,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        if (x[0] < outerLowerLeft_[0] + eps_) {
            Scalar a = -(1 + 0.5/height_);
            Scalar b = -a*outerUpperRight_[1];
            values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
            values[sNIdx] = outerSnr_;
        }
        else {
            Scalar a = -1;
            Scalar b = outerUpperRight_[1];
            values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
            values[sNIdx] = outerSnr_;
        }

        return values;
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        Scalar lambda = (outerUpperRight_[0] - x[0])/width_;
        if (lambda > 0.5 && lambda < 2.0/3.0 && x[1] > outerUpperRight_[1] - eps_) {
            values[sNIdx] = -0.04;
        }

        return values;
    }

    virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Entity& e,
                                               const FieldVector<Scalar,dim>& xi) const
    {

        FieldVector<Scalar,numEq> values;

        values[pWIdx] = -densityW_*gravity_[1]*(height_ - x[1]);

        if (x[0] < outerLowerLeft_[0] + eps_) {
            Scalar a = -(1 + 0.5/height_);
            Scalar b = -a*outerUpperRight_[1];
            values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
        }

        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            values[sNIdx] = innerSnr_;
        else
            values[sNIdx] = outerSnr_;

        return values;
    }

    double porosity (const FieldVector<Scalar,dim>& x, const Entity& e,
                     const FieldVector<Scalar,dim>& xi) const
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerPorosity_;
        else
            return outerPorosity_;
    }

    virtual FieldVector<Scalar,dim> gravity () const
    {
        return gravity_;
    }

    virtual FieldVector<Scalar,4> materialLawParameters (const FieldVector<Scalar,dim>& x, const Entity& e,
                                                         const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,4> values;

        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]) {
            values[swrIdx] = innerSwr_;
            values[snrIdx] = innerSnr_;
            values[alphaIdx] = innerAlpha_;
            values[nIdx] = innerN_;
        }
        else {
            values[swrIdx] = outerSwr_;
            values[snrIdx] = outerSnr_;
            values[alphaIdx] = outerAlpha_;
            values[nIdx] = outerN_;
        }

        return values;
    }

    LensProblem(TwoPhaseRelations& law = *(new DeprecatedLinearLaw),
                const FieldVector<Scalar,dim> outerLowerLeft = 0, const FieldVector<Scalar,dim> outerUpperRight = 0,
                const FieldVector<Scalar,dim> innerLowerLeft = 0, const FieldVector<Scalar,dim> innerUpperRight = 0,
                Scalar outerK = 4.6e-10, Scalar innerK = 9.05e-13,
                Scalar outerSwr = 0.05, Scalar outerSnr = 0, Scalar innerSwr = 0.18, Scalar innerSnr = 0,
                Scalar outerPorosity = 0.4, Scalar innerPorosity = 0.4,
                Scalar outerAlpha = 0.0037, Scalar innerAlpha = 0.00045,
                Scalar outerN = 4.7, Scalar innerN = 7.3)
        : TwoPhaseProblem<Grid, Scalar>(law),
          outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
          innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
          eps_(1e-8*outerUpperRight[0]),
          densityW_(law.wettingPhase.density()), densityN_(law.nonwettingPhase.density()),
          outerSwr_(outerSwr), outerSnr_(outerSnr), innerSwr_(innerSwr), innerSnr_(innerSnr),
          outerPorosity_(outerPorosity), innerPorosity_(innerPorosity),
          outerAlpha_(outerAlpha), innerAlpha_(innerAlpha),
          outerN_(outerN), innerN_(innerN)
    {
        outerK_[0][0] = outerK_[1][1] = outerK;
        outerK_[0][1] = outerK_[1][0] = 0;

        innerK_[0][0] = innerK_[1][1] = innerK;
        innerK_[0][1] = innerK_[1][0] = 0;

        height_ = outerUpperRight[1] - outerLowerLeft[1];
        width_ = outerUpperRight[0] - outerLowerLeft[0];

        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }

private:
    FieldMatrix<Scalar,dim,dim> outerK_;
    FieldMatrix<Scalar,dim,dim> innerK_;
    FieldVector<Scalar,dim> outerLowerLeft_;
    FieldVector<Scalar,dim> outerUpperRight_;
    FieldVector<Scalar,dim> innerLowerLeft_;
    FieldVector<Scalar,dim> innerUpperRight_;
    Scalar width_, height_;
    Scalar eps_;
    Scalar densityW_, densityN_;
    FieldVector<Scalar,dim> gravity_;
    Scalar outerSwr_, outerSnr_, innerSwr_, innerSnr_;
    Scalar outerPorosity_, innerPorosity_;
    Scalar outerAlpha_, innerAlpha_;
    Scalar outerN_, innerN_;
};

}
#endif
