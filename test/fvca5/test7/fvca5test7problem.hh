#ifndef FVCA5TEST7PROBLEM_HH
#define FVCA5TEST7PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC>
class FVCA5Test7Problem : public DeprecatedDiffusionProblem<G,RT,VC>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    FVCA5Test7Problem(VC& variables, double delta = 0.2)
        : DeprecatedDiffusionProblem<G,RT,VC>(variables)
    {
        delta_ = 0.2;
    }

    FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
                            const FieldVector<DT,n>& xi)
    {
        double phi1 = x[1] - delta_*(x[0] - 0.5) - 0.475;
        double phi2 = phi1 - 0.05;
        if (phi1 < 0) {
            permloc_[0][0] = 1.0;
            permloc_[0][1] = permloc_[1][0] = 0.0;
            permloc_[1][1] = 1.0;
        }
        else if (phi2 < 0) {
            permloc_[0][0] = 0.01;
            permloc_[0][1] = permloc_[1][0] = 0.0;
            permloc_[1][1] = 0.01;
        }
        else {
            permloc_[0][0] = 1.0;
            permloc_[0][1] = permloc_[1][0] = 0.0;
            permloc_[1][1] = 1.0;
        }

        return permloc_;
    }

    RT source   (const FieldVector<DT,n>& x, const Entity& e,
                 const FieldVector<DT,n>& xi)
    {
        return (0.0);
    }

    typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                                               const FieldVector<DT,n>& xi) const
    {
        return BoundaryConditions::dirichlet;
    }

    RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
                       const FieldVector<DT,n>& xi) const
    {
        return (exact(x));
    }

    RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const
    {
        return 0;
    }

    RT exact (const FieldVector<DT,n>& x) const
    {
        double phi1 = x[1] - delta_*(x[0] - 0.5) - 0.475;
        double phi2 = phi1 - 0.05;
        if (phi1 < 0) {
            return (-phi1);
        }
        else if (phi2 < 0) {
            return (- phi1/0.01);
        }
        else {
            return (- phi2 - 5.0);
        }
    }

    FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
    {
        FieldVector<RT,n> grad(0);
        double phi1 = x[1] - delta_*(x[0] - 0.5) - 0.475;
        double phi2 = phi1 - 0.05;
        if (phi1 < 0) {
            grad[0] = delta_;
            grad[1] = -1.0;
        }
        else if (phi2 < 0) {
            grad[0] = delta_/0.01;
            grad[1] = -1.0/0.01;
        }
        else {
            grad[0] = delta_;
            grad[1] = -1.0;
        }

        return grad;
    }

private:
    FieldMatrix<DT,n,n> permloc_;
    double delta_;
};
}

#endif
