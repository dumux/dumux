#ifndef DUNE_FIVESPOTCASE1PROBLEM_HH
#define DUNE_FIVESPOTCASE1PROBLEM_HH

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
 * @brief  Classes defining two cases of the Five-Spot problem (2-D)
 * @author Markus Wolff
 */

namespace Dune {

template<class G, class RT> class FiveSpotCase1Problem :
        public DeprecatedTwoPhaseProblem<G, RT> {
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=2};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::LeafGridView::IntersectionIterator
    IntersectionIterator;

public:
    enum {wPhaseIdx = 0, nwPhaseIdx = 1,pWIdx = 0, sNIdx = 1};
    enum {swrIdx = 0, snrIdx = 1};
    enum {VanGenuchten = 1, alphaIdx = 2, nIdx = 3};
    enum {BrooksCorey = 0,lambdaIdx = 2, p0Idx = 3};

    virtual const FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x,
                                         const Entity& e, const FieldVector<DT,n>& xi) {
        return K_;
    }

    virtual FieldVector<RT,m> q(const FieldVector<DT,n>& x, const Entity& e,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype(
                                                             const FieldVector<DT,n>& x, const Entity& e,
                                                             const IntersectionIterator& intersectionIt,
                                                             const FieldVector<DT,n>& xi) const {
        FieldVector<BoundaryConditions::Flags, m> values(
                                                         BoundaryConditions::neumann);

        if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_)||(x[1]
                                                                           < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_)) {
            values[wPhaseIdx] = BoundaryConditions::dirichlet;
            values[nwPhaseIdx] = BoundaryConditions::dirichlet;
        }
        //        if ((x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_)
        //                ||(x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_)) {
        //            values[wPhaseIdx] = BoundaryConditions::dirichlet;
        //            values[nwPhaseIdx] = BoundaryConditions::dirichlet;
        //        }
        return values;
    }

    virtual FieldVector<RT,m>dirichlet(const FieldVector<DT,n>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);

        if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_)||(x[1]
                                                                           < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_)) {
            values[pWIdx] = pwlowerleftbc_;
            values[sNIdx] = Snr_;
        }

        //        if ((x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_)
        //                ||(x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_)) {
        //            values[pWIdx] = pwupperrightbc_;
        //            values[sNIdx] = 1-Swr_;
        //
        //        }

        return values;
    }

    virtual FieldVector<RT,m> neumann(const FieldVector<DT,n>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);

        //      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) ||
        //         (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_)) {
        //         values[wPhaseIdx] = -0.12;
        //      }
        if ((x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_)
            ||(x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_)) {
            values[nwPhaseIdx] = 0.001;//15x15 cells: bcf 21;
            //    values[nwPhaseIdx] = 0.002;// 30x30 cells: bcf 11
            //    values[nwPhaseIdx] = 0.004;// 60x60 cells: bcf 6
        }
        return values;
    }

    virtual FieldVector<RT,m> initial(const FieldVector<DT,n>& x,
                                      const Entity& e, const FieldVector<DT,n>& xi) const {

        FieldVector<RT,m> values;

        values[pWIdx] = pwlowerleftbc_;
        values[sNIdx] = 1-Swr_;

        return values;
    }

    double porosity(const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const {
        return Porosity_;
    }

    virtual FieldVector<RT,n> gravity() const {
        return gravity_;
    }

    virtual FieldVector<RT,4> materialLawParameters(const FieldVector<DT,n>& x,
                                                    const Entity& e, const FieldVector<DT,n>& xi) const {
        FieldVector<RT,4> values;

        if (chooselaw_) {
            values[swrIdx] = Swr_;
            values[snrIdx] = Snr_;
            values[alphaIdx] = Alpha_;
            values[nIdx] = N_;
        } else {
            values[swrIdx] = Swr_;
            values[snrIdx] = Snr_;
            values[lambdaIdx] = Lambda_;
            values[p0Idx] = p0_;
        }
        return values;
    }

    FiveSpotCase1Problem(DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw), const FieldVector<DT,n> LowerLeft = 0,
                         const FieldVector<DT,n> UpperRight = 0, RT bcf = 11,
                         int chooselaw = BrooksCorey, RT K = 1e-7, RT Swr = 0.2,
                         RT Snr = 0.2, RT Porosity = 0.2, RT Lambda = 2.0, RT p0 = 0,
                         RT Alpha = 0.0037, RT N = 4.7, RT pwlowerleftbc=2e5,
                         RT pwupperrightbc=1.999986e5) :
        DeprecatedTwoPhaseProblem<G, RT>(law),
        LowerLeft_(LowerLeft),
        UpperRight_(UpperRight),
        eps_(1e-8*UpperRight[0]),
        bcf_(bcf),
        densityW_(law.wettingPhase.density()),
        densityN_(law.nonwettingPhase.density()),
        Swr_(Swr),
        Snr_(Snr),
        Porosity_(Porosity),
        Lambda_(Lambda),
        p0_(p0),
        Alpha_(Alpha),
        N_(N),
        pwlowerleftbc_(pwlowerleftbc),
        pwupperrightbc_(pwupperrightbc),
        chooselaw_(chooselaw)
    {
        K_[0][0]=K_[1][1]=K;
        K_[1][0]=K_[0][1]=0;

        height_ = UpperRight[1] - LowerLeft[1];
        width_ = UpperRight[0] - LowerLeft[0];

        gravity_[0] = 0;
        gravity_[1] = 0;
    }

private:
    FieldMatrix<DT,n,n> K_;

    FieldVector<DT,n> LowerLeft_;
    FieldVector<DT,n> UpperRight_;

    DT width_, height_;
    DT eps_;
    RT bcf_;
    RT densityW_, densityN_;
    FieldVector<DT,n> gravity_;
    RT Swr_, Snr_;
    RT Porosity_;
    RT Lambda_;
    RT p0_;
    RT Alpha_;
    RT N_;
    RT pwlowerleftbc_, pwupperrightbc_;
    int chooselaw_;
};
template<class G, class RT> class FiveSpotCase2Problem :
        public DeprecatedTwoPhaseProblem<G, RT> {
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=2};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::LeafGridView::IntersectionIterator
    IntersectionIterator;

public:
    enum {wPhaseIdx = 0, nwPhaseIdx = 1,pWIdx = 0, sNIdx = 1};
    enum {swrIdx = 0, snrIdx = 1};
    enum {VanGenuchten = 1, alphaIdx = 2, nIdx = 3};
    enum {BrooksCorey = 0,lambdaIdx = 2, p0Idx = 3};

    virtual const FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x,
                                         const Entity& e, const FieldVector<DT,n>& xi) {
        return K_;
    }

    virtual FieldVector<RT,m> q(const FieldVector<DT,n>& x, const Entity& e,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype(
                                                             const FieldVector<DT,n>& x, const Entity& e,
                                                             const IntersectionIterator& intersectionIt,
                                                             const FieldVector<DT,n>& xi) const {
        FieldVector<BoundaryConditions::Flags, m> values(
                                                         BoundaryConditions::neumann);

        //        if ((x[0] < LowerLeft_[0] + eps_ && x[1] > UpperRight_[1] - bcf_) || //upper left
        //                (x[1] > UpperRight_[1] - eps_ && x[0] < LowerLeft_[0] + bcf_)|| //upper left
        //                (x[0] > UpperRight_[0] - eps_ && x[1] < LowerLeft_[1] + bcf_)|| //lower right
        //                (x[1] < LowerLeft_[1] + eps_ && x[0] > UpperRight_[0] - bcf_)) { //lower right
        //            values[wPhaseIdx] = BoundaryConditions::dirichlet;
        //            values[nwPhaseIdx] = BoundaryConditions::dirichlet;
        //        }
        if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
            (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_)|| //lower left
            (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_)
            || //upper right
            (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_)) { //upper right
            values[wPhaseIdx] = BoundaryConditions::dirichlet;
            values[nwPhaseIdx] = BoundaryConditions::dirichlet;
        }
        return values;
    }

    virtual FieldVector<RT,m>dirichlet(const FieldVector<DT,n>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);

        //        if ((x[0] < LowerLeft_[0] + eps_ && x[1] > UpperRight_[1] - bcf_)
        //                ||(x[1] > UpperRight_[1] - eps_ && x[0] < LowerLeft_[0] + bcf_)
        //                ||(x[0] > UpperRight_[0] - eps_ && x[1] < LowerLeft_[1] + bcf_)
        //                ||(x[1] < LowerLeft_[1] + eps_ && x[0] > UpperRight_[0] - bcf_)) {
        //            values[pWIdx] = pwoutbc_;
        //            values[sNIdx] = 1-Swr_;
        //        }
        if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
            (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_)|| //lower left
            (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_)|| //upper right
            (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_)) { //upper right
            values[pWIdx] = pwinbc_;
            values[sNIdx] = Snr_;
        }
        return values;
    }

    virtual FieldVector<RT,m> neumann(const FieldVector<DT,n>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,n>& xi) const {
        FieldVector<RT,m> values(0);
        if ((x[0] < LowerLeft_[0] + eps_ && x[1] > UpperRight_[1] - bcf_)
            ||(x[1] > UpperRight_[1] - eps_ && x[0] < LowerLeft_[0] + bcf_)
            ||(x[0] > UpperRight_[0] - eps_ && x[1] < LowerLeft_[1] + bcf_)
            ||(x[1] < LowerLeft_[1] + eps_ && x[0] > UpperRight_[0] - bcf_)) {
            values[nwPhaseIdx] = 0.001;//15x15 cells: bcf 21;
            //            values[nwPhaseIdx] = 0.002;// 30x30 cells: bcf 11
            //            values[nwPhaseIdx] = 0.004;// 60x60 cells: bcf 6
        }
        return values;
    }

    virtual FieldVector<RT,m> initial(const FieldVector<DT,n>& x,
                                      const Entity& e, const FieldVector<DT,n>& xi) const {

        FieldVector<RT,m> values;

        values[pWIdx] = pwinbc_;
        values[sNIdx] = 1-Swr_;

        return values;
    }

    double porosity(const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const {
        return Porosity_;
    }

    virtual FieldVector<RT,n> gravity() const {
        return gravity_;
    }

    virtual FieldVector<RT,4> materialLawParameters(const FieldVector<DT,n>& x,
                                                    const Entity& e, const FieldVector<DT,n>& xi) const {
        FieldVector<RT,4> values;

        if (chooselaw_) {
            values[swrIdx] = Swr_;
            values[snrIdx] = Snr_;
            values[alphaIdx] = Alpha_;
            values[nIdx] = N_;
        } else {
            values[swrIdx] = Swr_;
            values[snrIdx] = Snr_;
            values[lambdaIdx] = Lambda_;
            values[p0Idx] = p0_;
        }
        return values;
    }

    FiveSpotCase2Problem(DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw), const FieldVector<DT,n> LowerLeft = 0,
                         const FieldVector<DT,n> UpperRight = 0, RT bcf = 11, int chooselaw = BrooksCorey, RT K = 1e-7,
                         RT Swr = 0.2, RT Snr = 0.2, RT Porosity = 0.2, RT Lambda = 2.0,
                         RT p0 = 0, RT Alpha = 0.0037, RT N = 4.7, RT pwinbc=2e5,
                         RT pwoutbc=2e5) :
        DeprecatedTwoPhaseProblem<G, RT>(law), LowerLeft_(LowerLeft),
        UpperRight_(UpperRight), eps_(1e-8*UpperRight[0]), bcf_(bcf),
        densityW_(law.wettingPhase.density()),
        densityN_(law.nonwettingPhase.density()),
        Swr_(Swr), Snr_(Snr),
        Porosity_(Porosity), Lambda_(Lambda), p0_(p0), Alpha_(Alpha),
        N_(N),
        pwinbc_(pwinbc),
        pwoutbc_(pwoutbc),
        chooselaw_(chooselaw)
    {
        K_[0][0]=K_[1][1]=K;
        K_[1][0]=K_[0][1]=0;

        height_ = UpperRight[1] - LowerLeft[1];
        width_ = UpperRight[0] - LowerLeft[0];

        gravity_[0] = 0;
        gravity_[1] = 0;
    }

private:
    FieldMatrix<DT,n,n> K_;

    FieldVector<DT,n> LowerLeft_;
    FieldVector<DT,n> UpperRight_;

    DT width_, height_;
    DT eps_;
    RT bcf_;
    RT densityW_, densityN_;
    FieldVector<DT,n> gravity_;
    RT Swr_, Snr_;
    RT Porosity_;
    RT Lambda_;
    RT p0_;
    RT Alpha_;
    RT N_;
    RT pwinbc_, pwoutbc_;
    int chooselaw_;
};
}
#endif
