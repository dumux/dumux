#ifndef MCWHORTERDIFFPROBLEM_HH
#define MCWHORTERDIFFPROBLEM_HH

#include "dumux/diffusion/problems/homogeneousproblem.hh"

namespace Dune
{
template<class G,class RT, class VC>
class McWhorterDiffProblem : public HomogeneousProblem<G,RT,VC>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    McWhorterDiffProblem(VC& variableobj, DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw),
                         const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0,
                         RT pleftbc=1.95e5, RT prightbc=0, const bool cap = false)
        : HomogeneousProblem<G,RT,VC>(variableobj, law, cap,1e-10),
          Left_(Left[0]), Right_(Right[0]),
          eps_(1e-8),
          pleftbc_(pleftbc),prightbc_(prightbc)
    {}

    RT source (const FieldVector<DT,n>& x, const Entity& e,
               const FieldVector<DT,n>& xi)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                                               const FieldVector<DT,n>& xi) const
    {
        if (x[0] < eps_)
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;

    }

    RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
                       const FieldVector<DT,n>& xi) const
    {
        if (x[0] < eps_)
            return pleftbc_;
        //all other boundaries
        return prightbc_;
    }

    RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const
    {
        if (x[0] < eps_)
            return 1;
        else
            return 0;
    }

    RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const
    {
        return 0;
    }
private:
    DT Left_;
    DT Right_;

    RT eps_;
    RT pleftbc_, prightbc_;
};
}

#endif
