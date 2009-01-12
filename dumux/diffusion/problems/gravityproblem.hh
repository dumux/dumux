// $Id$

#ifndef GRAVITYPROBLEM_HH
#define GRAVITYPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
    template<class G, class RT, class VC>
    class GravityProblem : public DiffusionProblem<G,RT,VC>
    {
      typedef typename G::ctype DT;
      enum {n=G::dimension};
      typedef typename G::Traits::template Codim<0>::Entity Entity;

    public:
      GravityProblem(VC& variables, G* g, TwoPhaseRelations& law = *(new LinearLaw), FieldVector<DT,n> gravity = *(new FieldVector<DT,n>(0)))
          : DiffusionProblem<G,RT,VC>(variables, law, false, gravity), grid(g)
      { }

      GravityProblem()
        : DiffusionProblem<G,RT,VC>()
      { }

      FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
                      const FieldVector<DT,n>& xi)
      {
          permloc[0][0] = permloc[1][1] = 1e-10;
          permloc[0][1] = permloc[1][0] = 0;

          return permloc;
      }

      RT sat (const Dune::FieldVector<DT,n>& x, const Entity& e,
                      const Dune::FieldVector<DT,n>& xi)
      {
              return (*saturation)[grid->levelIndexSet(e.level()).index(e)];
      }

      RT source   (const FieldVector<DT,n>& x, const Entity& e,
                          const FieldVector<DT,n>& xi)
      {
              return 0;
      }

      typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                           const FieldVector<DT,n>& xi) const
      {
        if (x[1] > 10-1E-6)// || x[1] < 1e-6)
          return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
      }

      RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const
      {
          //if (x[1] > 10-1e-6)
              return (2e5);
      }


      RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const
      {
//          if (x[1] < 1e-8)
//              return (((this->gravity_)[n-1]));
//          else
              return 0;
      }

    private:
        FieldMatrix<DT,n,n> permloc;
        G* grid;
    public:
        Dune::BlockVector<Dune::FieldVector<RT,1> >* saturation;
    };
}

#endif
