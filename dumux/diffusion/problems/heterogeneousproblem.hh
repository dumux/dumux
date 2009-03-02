// $Id$

#ifndef HETEROGENPROBLEM_HH
#define HETEROGENPROBLEM_HH

#include "dumux/diffusion/diffusionproblem_deprecated.hh"
#include "dumux/material/randompermeability.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC>
class HeterogeneousProblem : public DeprecatedDiffusionProblem<G,RT,VC>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    HeterogeneousProblem(VC& variableobj, const G& g, const char* name = "permeab.dat", const bool create = true,
                         DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw), const bool cap = false)
        : DeprecatedDiffusionProblem<G,RT,VC>(variableobj,law, cap), permeability(g, name, create)
    { }

    Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                  const Dune::FieldVector<DT,n>& xi)
    {
        return permeability.K(e);
    }

    RT source   (const Dune::FieldVector<DT,n>& x, const Entity& e,
                 const Dune::FieldVector<DT,n>& xi)
    {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                                     const Dune::FieldVector<DT,n>& xi) const
    {
        if (x[0] > 300-1E-6 || x[0] < 1e-6)
            return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    RT dirichletPress (const Dune::FieldVector<DT,n>& x, const Entity& e,
                       const Dune::FieldVector<DT,n>& xi) const
    {
        return (x[0] < 1e-6) ? 1e6 : 0;
    }

    RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const
    {
        if (x[0] < 1e-6)
            return 0.8;
        else
            return 0.2;
    }

    RT neumannPress (const Dune::FieldVector<DT,n>& x, const Entity& e,
                     const Dune::FieldVector<DT,n>& xi) const
    {
        return 0;
    }

    RandomPermeability<G> permeability;
};
}

#endif
