#ifndef INJECTIONDIFFPROBLEM_HH
#define INJECTIONDIFFPROBLEM_HH

#include "dumux/diffusion/problems/homogeneousproblem.hh"

namespace Dune
{
    //! \ingroup diffusionProblems
    //! example class for diffusion problems
    template<class G, class RT, class VC> class InjectionDiffProblem :
        public HomogeneousProblem<G,RT,VC>
    {
        typedef typename G::ctype DT;
        enum
        {    n=G::dimension};
        typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
        InjectionDiffProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), const FieldVector<DT,n> outerLowerLeft = 0.,
                const FieldVector<DT,n> outerUpperRight = 0.,
                const FieldVector<DT,n> innerLowerLeft = 0.,
                const FieldVector<DT,n> innerUpperRight = 0., const RT depthBOR = 0.,
                const bool cap = false, RT outerK = 7.2e-13, RT innerK = 7.2e-13) :
            HomogeneousProblem<G, RT, VC>(variableobj, law, cap), outerLowerLeft_(outerLowerLeft),
                    outerUpperRight_(outerUpperRight), innerLowerLeft_(innerLowerLeft),
                    innerUpperRight_(innerUpperRight), depthBOR_(depthBOR),
                    densityW_(law.wettingPhase.density()),
                    densityN_(law.nonwettingPhase.density()), eps_(1e-8)
        {
            outerK_[0][0] = outerK_[1][1] = outerK;
            outerK_[0][1] = outerK_[1][0] = 0;

            innerK_[0][0] = innerK_[1][1] = innerK;
            innerK_[0][1] = innerK_[1][0] = 0;

            gravity_[0] = 0;
            gravity_[1] = -9.81;
        }
        // permeabilities
                virtual FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x, const Entity& e,
                        const FieldVector<DT,n>& xi)
                {
                    if (x[0] >= innerLowerLeft_[0]&& x[0] <= innerUpperRight_[0]&& x[1]
                            >= innerLowerLeft_[1]&& x[1] <= innerUpperRight_[1])
                        return innerK_;
                    else
                        return outerK_;
                }

                virtual RT source(const FieldVector<DT,n>& x, const Entity& e,
                        const FieldVector<DT,n>& xi)
                {
                        return 0;
                }

        typename BoundaryConditions::Flags bctype(const FieldVector<DT,n>& x, const Entity& e,
                const FieldVector<DT,n>& xi) const
        {
            if (x[0] < outerLowerLeft_[0] + eps_)
                return BoundaryConditions::dirichlet;
            // all other boundaries
            return BoundaryConditions::neumann;
        }

        RT dirichletPress(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
        {
            if (x[0] < outerLowerLeft_[0] + eps_)
                return (1e5 - densityW_*gravity_[1]*(depthBOR_ - x[1]));
            return 0;
        }

        RT dirchletSat(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
        {
               if (x[0] < outerLowerLeft_[0] + eps_)
                    return 1;
            return 0;
        }

        RT neumannPress(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
        {
            if (x[1] < 25 && x[1] > 15)
                return -1.54e-5;
            return 0;
        }

        const FieldVector<DT,n>& gravity() const
        {
            return gravity_;
        }

private:
        FieldMatrix<DT,n,n> outerK_;
        FieldMatrix<DT,n,n> innerK_;
        FieldVector<DT,n> outerLowerLeft_;
        FieldVector<DT,n> outerUpperRight_;
        FieldVector<DT,n> innerLowerLeft_;
        FieldVector<DT,n> innerUpperRight_;
        DT depthBOR_;
        RT densityW_, densityN_;
        FieldVector<DT,n> gravity_;
        RT eps_;
    };
}

#endif
