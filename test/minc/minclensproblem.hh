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
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>
#include<dumux/minc/mincproblem.hh>

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
   *    - RT    type used for return values
   */
  template<class G, class RT, int m>
  class MincLensProblem : public MincProblem<G, RT, m> {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
    enum {F = 0, M = 1};
    enum {pWFIdx = 0, sNFIdx = 1, pWMIdx = 2, sNMIdx = 3};
    enum {swrIdx = 0, snrIdx = 1, alphaIdx = 2, nIdx = 3};
    enum {nCont=10};

    virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerKFracture_;
        else
            return outerKFracture_;
    }
    virtual const FieldMatrix<DT,n,n>& K1Fracture (const FieldVector<DT,n>& x)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerKFracture_;
        else
            return outerKFracture_;
    }
    virtual const FieldMatrix<DT,n,n>& K1Matrix (const FieldVector<DT,n>& x)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerKMatrix_;
        else
            return outerKMatrix_;
    }


    virtual const FieldMatrix<DT,n,n>& KFracture (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerKFracture_;
        else
            return outerKFracture_;
    }

    virtual const FieldMatrix<DT,n,n>& KMatrix (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi)
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerKMatrix_;
        else
            return outerKMatrix_;
    }

    virtual FieldVector<RT,m> q (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e,
                    const IntersectionIterator& intersectionIt,
                       const FieldVector<DT,n>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

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

    virtual FieldVector<RT,m> g (const FieldVector<DT,n>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);

        for (int nC=0; nC < m; nC+=2){

            if (x[0] > outerUpperRight_[0] - eps_) {
            RT a = -1;
            RT b = outerUpperRight_[1];
            // Fracture domain
            if (nC<2){
            values[pWFIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
            values[sNFIdx] = outerSnrFracture_+0.0;}
            // Matrix domain
            else {
            values[nC] = values[pWFIdx]*1.0;            //pWM = pWF
            values[nC+1] = outerSnrMatrix_+0.0;
            }
            }
        }
        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<DT,n>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                  const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);

//        RT lambda = (outerUpperRight_[0] - x[0])/width_;
            if (x[0] < outerLowerLeft_[0] + eps_) {
            values[sNFIdx] = -0.5;
            values[pWFIdx] = -0.0;
            }

        return values;
    }

    virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e,
                  const FieldVector<DT,n>& xi) const
    {

        FieldVector<RT,m> values;
        for (int nC=0; nC<m; nC+=4){
        values[nC] = -densityW_*gravity_[1]*(outerUpperRight_[1] - x[1])+0;
        values[nC+1] = outerSnrFracture_+0.0;
        values[nC+2] = (-densityW_*gravity_[1]*(outerUpperRight_[1] - x[1])+0)*0.9;
        values[nC+3] = outerSnrMatrix_+0.0;
        }
        return values;
    }

//    double porosity (const FieldVector<DT,n>& x, const Entity& e,
//              const FieldVector<DT,n>& xi) const
//    {
//        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
//            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
//            return innerPorosity_;
//        else
//            return outerPorosity_;
//    }

    double porosityFracture (const FieldVector<DT,n>& x, const Entity& e,
              const FieldVector<DT,n>& xi) const
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerPorosityFracture_;
        else
            return outerPorosityFracture_;
    }

    double porosityMatrix (const FieldVector<DT,n>& x, const Entity& e,
              const FieldVector<DT,n>& xi) const
    {
        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
            return innerPorosityMatrix_;
        else
            return outerPorosityMatrix_;
    }


    virtual FieldVector<RT,n> gravity () const
    {
        return gravity_;
    }

    virtual FieldVector<RT,4> materialLawParametersFracture (const FieldVector<DT,n>& x, const Entity& e,
              const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,4> values;

        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]) {
            values[swrIdx] = innerSwrFracture_;
            values[snrIdx] = innerSnrFracture_;
            values[alphaIdx] = innerAlphaFracture_;
            values[nIdx] = innerNFracture_;
        }
        else {
            values[swrIdx] = outerSwrFracture_;
            values[snrIdx] = outerSnrFracture_;
            values[alphaIdx] = outerAlphaFracture_;
            values[nIdx] = outerNFracture_;
        }

        return values;
    }
//
    virtual FieldVector<RT,4> materialLawParametersMatrix (const FieldVector<DT,n>& x, const Entity& e,
              const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,4> values;

        if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
            && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]) {
            values[swrIdx] = innerSwrMatrix_;
            values[snrIdx] = innerSnrMatrix_;
            values[alphaIdx] = innerAlphaMatrix_;
            values[nIdx] = innerNMatrix_;
        }
        else {
            values[swrIdx] = outerSwrMatrix_;
            values[snrIdx] = outerSnrMatrix_;
            values[alphaIdx] = outerAlphaMatrix_;
            values[nIdx] = outerNMatrix_;
        }

        return values;
    }

    MincLensProblem(DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw),
            const FieldVector<DT,n> outerLowerLeft = 0, const FieldVector<DT,n> outerUpperRight = 0,
            const FieldVector<DT,n> innerLowerLeft = 0, const FieldVector<DT,n> innerUpperRight = 0,
            RT outerKFracture = 4.0e-10, RT innerKFracture = 9.0e-13,
            RT outerKMatrix = 4.0e-10, RT innerKMatrix = 9.0e-13,
            RT outerSwrFracture = 0.05, RT outerSnrFracture = 0,
            RT innerSwrFracture = 0.18, RT innerSnrFracture = 0,
            RT outerSwrMatrix = 0.05, RT outerSnrMatrix = 0,
            RT innerSwrMatrix = 0.18, RT innerSnrMatrix = 0,

            RT outerPorosityFracture = 0.2, RT innerPorosityFracture = 0.2,
            RT outerPorosityMatrix = 0.2, RT innerPorosityMatrix = 0.2,
            RT outerAlphaFracture = 0.0037, RT innerAlphaFracture = 0.00045,
            RT outerAlphaMatrix = 0.0037, RT innerAlphaMatrix = 0.00045,
            RT outerNFracture = 4.7, RT innerNFracture = 7.3,
            RT outerNMatrix = 4.7, RT innerNMatrix = 7.3)

    : MincProblem<G, RT, m>(law),
      outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
      innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
      eps_(1e-8*outerUpperRight[0]),
      densityW_(law.wettingPhase.density()), densityN_(law.nonwettingPhase.density()),
      outerSwrFracture_(outerSwrFracture), outerSnrFracture_(outerSnrFracture),
      innerSwrFracture_(innerSwrFracture), innerSnrFracture_(innerSnrFracture),
      outerSwrMatrix_(outerSwrMatrix), outerSnrMatrix_(outerSnrMatrix),
      innerSwrMatrix_(innerSwrMatrix), innerSnrMatrix_(innerSnrMatrix),
      outerPorosityFracture_(outerPorosityFracture), innerPorosityFracture_(innerPorosityFracture),
      outerPorosityMatrix_(outerPorosityMatrix), innerPorosityMatrix_(innerPorosityMatrix),
      outerAlphaFracture_(outerAlphaFracture), innerAlphaFracture_(innerAlphaFracture),
      outerAlphaMatrix_(outerAlphaMatrix), innerAlphaMatrix_(innerAlphaMatrix),
      outerNFracture_(outerNFracture), innerNFracture_(innerNFracture),
      outerNMatrix_(outerNMatrix), innerNMatrix_(innerNMatrix)
    {
        outerKFracture_[0][0] = outerKFracture_[1][1] = outerKFracture;
        outerKFracture_[0][1] = outerKFracture_[1][0] = 0;

        outerKMatrix_[0][0] = outerKMatrix_[1][1] = outerKMatrix;
        outerKMatrix_[0][1] = outerKMatrix_[1][0] = 0;

        innerKFracture_[0][0] = innerKFracture_[1][1] = innerKFracture;
        innerKFracture_[0][1] = innerKFracture_[1][0] = 0;

        innerKMatrix_[0][0] = innerKMatrix_[1][1] = innerKMatrix;
        innerKMatrix_[0][1] = innerKMatrix_[1][0] = 0;


        height_ = outerUpperRight[1] - outerLowerLeft[1];
        width_ = outerUpperRight[0] - outerLowerLeft[0];

        gravity_[0] = 0;
        gravity_[1] = -9.81;
        volumefraction= 0.7; //The interacting continua volume is 70% of the controlvolume.
    }

    private:
        FieldMatrix<DT,n,n> outerK_;
        FieldMatrix<DT,n,n> innerK_;
        FieldMatrix<DT,n,n> outerKFracture_;
        FieldMatrix<DT,n,n> innerKFracture_;
        FieldMatrix<DT,n,n> outerKMatrix_;
        FieldMatrix<DT,n,n> innerKMatrix_;
        FieldVector<DT,n> outerLowerLeft_;
        FieldVector<DT,n> outerUpperRight_;
        FieldVector<DT,n> innerLowerLeft_;
        FieldVector<DT,n> innerUpperRight_;
        DT width_, height_;
        DT eps_;
        RT densityW_, densityN_;
        FieldVector<DT,n> gravity_;
        RT volumefraction;

//        RT outerSwr_, outerSnr_, innerSwr_, innerSnr_;
        RT outerSwrFracture_, outerSnrFracture_, innerSwrFracture_, innerSnrFracture_;
        RT outerSwrMatrix_, outerSnrMatrix_, innerSwrMatrix_, innerSnrMatrix_;
//        RT outerPorosity_, innerPorosity_;
        RT outerPorosityFracture_, innerPorosityFracture_;
        RT outerPorosityMatrix_, innerPorosityMatrix_;
        RT outerAlphaFracture_, innerAlphaFracture_;
        RT outerAlphaMatrix_, innerAlphaMatrix_;
        RT outerNFracture_, innerNFracture_;
        RT outerNMatrix_, innerNMatrix_;
  };

}
#endif
