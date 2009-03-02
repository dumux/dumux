// $Id$

#ifndef PARABOLICPROBLEM_HH
#define PARABOLICPROBLEM_HH

#include "dumux/parabolic/parabolicproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class G, class RT>
class UniformParabolicProblem : public ParabolicProblem<G,RT>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    UniformParabolicProblem(G& g, DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw), const bool cap = false)
        : ParabolicProblem<G,RT>(law, cap)
    { }

    UniformParabolicProblem()
        : ParabolicProblem<G,RT>()
    { }

    const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                        const Dune::FieldVector<DT,n>& xi)
    {
        permloc[0][0] = permloc[1][1] = 1;
        permloc[0][1] = permloc[1][0] = 0;

        return permloc;
    }

    RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
            const Dune::FieldVector<DT,n>& xi) const
    {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                                     const Dune::FieldVector<DT,n>& xi) const
    {
        if (x[0] > 600-1E-6 || x[0] < 1e-6)
            //if (x[0] > 600-1E-6)
            return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
          const Dune::FieldVector<DT,n>& xi) const
    {
        return (x[0] < 1e-6) ? 1 : 0;
    }


    RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
          const Dune::FieldVector<DT,n>& xi) const
    {
        return 0;//(x[0] < 1e-6) ? -1.0/600.0 : 0;
    }

    virtual RT initial (const FieldVector<DT,n>& x, const Entity& e,
                        const FieldVector<DT,n>& xi) const
    {
        return (x[0] < 1e-6) ? 1 : 0;
    }

private:
    Dune::FieldMatrix<DT,n,n> permloc;
};
}

#endif
