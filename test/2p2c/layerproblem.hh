// $Id$

#ifndef DUNE_LAYERPROBLEM_HH
#define DUNE_LAYERPROBLEM_HH

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
#include<dumux/2p2c/2p2cproblem.hh>

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
  class LayerProblem : public TwoPTwoCProblem<Grid, Scalar> {
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
    enum {pWIdx = 0, satNIdx = 1};
    enum {swrIdx = 0, snrIdx = 1, alphaIdx = 2, nIdx = 3};
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};

    // permeabilities
    virtual const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos)
    //, const Element& element, const FieldVector<Scalar,dim>& localPos)
    {
        if (globalPos[0] >= innerLowerLeft_[0] && globalPos[0] <= innerUpperRight_[0]
            && globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1])
            return innerK_;
        else
            return outerK_;
    }

    // sources and sinks
    virtual FieldVector<Scalar,numEq> q (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                    const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

/////////////////////////////
// TYPE of the boundaries
/////////////////////////////
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                    const IntersectionIterator& intersectionIt,
                       const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (globalPos[0] < outerLowerLeft_[0] + eps_)
            values = BoundaryConditions::dirichlet;
//        if (globalPos[1] < eps_)
//            values = BoundaryConditions::dirichlet;

        return values;
    }

/////////////////////////////
// INITIAL values
/////////////////////////////
        virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                      const FieldVector<Scalar,dim>& localPos) const
        {

            FieldVector<Scalar,numEq> values;

            values[pWIdx] = -densityW_*gravity_[1]*(depthBOR_ - globalPos[1]);

//            if (globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]
//             && globalPos[0] >= innerLowerLeft_[0])
                values[satNIdx] = 0.2;
//            else
//                values[satNIdx] = 1e-6;

            return values;
        }


        int initialPhaseState (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                      const FieldVector<Scalar,dim>& localPos) const
        {

            enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
            int state;

//            if (globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]
//                  && globalPos[0] >= innerLowerLeft_[0])
//                state = 2;
//            else
                state = 2;

            return state;
        }


/////////////////////////////
// DIRICHLET boundaries
/////////////////////////////
    virtual FieldVector<Scalar,numEq> g (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        values[pWIdx] = -densityW_*gravity_[1]*(depthBOR_ - globalPos[1]);

//        if (globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]
//         && globalPos[0] >= innerLowerLeft_[0])
            values[satNIdx] = 0.2;
//        else
//            values[satNIdx] = 1e-6;

        return values;
    }

/////////////////////////////
// NEUMANN boundaries
/////////////////////////////
    virtual FieldVector<Scalar,numEq> J (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        //Scalar lambda = (globalPos[1])/height_;
//        if (globalPos[1] < 2.0 && globalPos[1] > 1.5)
            values[satNIdx] = -5e-6;

        return values;
    }
//////////////////////////////

    double porosity (const FieldVector<Scalar,dim>& globalPos, const Element& element,
              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] > innerLowerLeft_[0] && globalPos[0] < innerUpperRight_[0]
            && globalPos[1] > innerLowerLeft_[1] && globalPos[1] < innerUpperRight_[1])
            return innerPorosity_;
        else
            return outerPorosity_;
    }

    virtual FieldVector<Scalar,dim> gravity () const
    {
        return gravity_;
    }

    double depthBOR () const
    {
        return depthBOR_;
    }

    virtual FieldVector<Scalar,4> materialLawParameters (const FieldVector<Scalar,dim>& globalPos, const Element& element,
              const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,4> values;

        if (globalPos[0] > innerLowerLeft_[0] && globalPos[0] < innerUpperRight_[0]
            && globalPos[1] > innerLowerLeft_[1] && globalPos[1] < innerUpperRight_[1]) {
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

    LayerProblem(TwoPhaseRelations& law = *(new DeprecatedLinearLaw), MultiComp& multicomp = *(new CWaterAir),
            const FieldVector<Scalar,dim> outerLowerLeft = 0., const FieldVector<Scalar,dim> outerUpperRight = 0.,
            const FieldVector<Scalar,dim> innerLowerLeft = 0., const FieldVector<Scalar,dim> innerUpperRight = 0.,
            const Scalar depthBOR = 0., Scalar outerK = 1.2e-12, Scalar innerK = 1.2e-12,
            Scalar outerSwr = 0.05, Scalar outerSnr = 0.1, Scalar innerSwr = 0.05, Scalar innerSnr = 0.1,
            Scalar outerPorosity = 0.4, Scalar innerPorosity = 0.4,
            Scalar outerAlpha = 0.0037, Scalar innerAlpha = 0.0037,  //0.00045
            Scalar outerN = 4.7, Scalar innerN = 4.7)    //7.3
    : TwoPTwoCProblem<Grid, Scalar>(law, multicomp),
      outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
      innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
      depthBOR_(depthBOR), eps_(1e-8*outerUpperRight[0]),
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
        Scalar depthBOR_, eps_;
        Scalar densityW_, densityN_;
        FieldVector<Scalar,dim> gravity_;
        Scalar outerSwr_, outerSnr_, innerSwr_, innerSnr_;
        Scalar outerPorosity_, innerPorosity_;
        Scalar outerAlpha_, innerAlpha_;
        Scalar outerN_, innerN_;
  };

}
#endif
