#ifndef FVCA5TEST5PROBLEM_HH
#define FVCA5TEST5PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
    template<class G, class RT, class VC>
    class FVCA5Test5Problem : public DiffusionProblem<G,RT,VC>
    {
      typedef typename G::ctype DT;
      enum {n=G::dimension};
      typedef typename G::Traits::template Codim<0>::Entity Entity;

    public:
      FVCA5Test5Problem(VC& variables)
        : DiffusionProblem<G,RT,VC>(variables)
      {  }

      FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
                      const FieldVector<DT,n>& xi)
      {
          double delta = 1.0e-3;
//          double pi = M_PI;
          double rt = x[0]*x[0]+x[1]*x[1];

          permloc_[0][0] = (delta*x[0]*x[0] + x[1]*x[1])/rt;
          permloc_[0][1] = permloc_[1][0] = -(1.0 - delta)*x[0]*x[1]/rt;
          permloc_[1][1] = (x[0]*x[0] + delta*x[1]*x[1])/rt;

          return permloc_;
      }

      RT source   (const FieldVector<DT,n>& x, const Entity& e,
                      const FieldVector<DT,n>& xi)
      {
          double delta = 1.0e-3;
          double pi = 4.0*atan(1.0);
          double rt = x[0]*x[0]+x[1]*x[1];
          double ux = pi*cos(pi*x[0])*sin(pi*x[1]);
          double uy = pi*cos(pi*x[1])*sin(pi*x[0]);
          double kxx = (delta*x[0]*x[0] + x[1]*x[1])/rt;
          double kxy = -(1.0 - delta)*x[0]*x[1]/rt;
          double kyy = (x[0]*x[0] + delta*x[1]*x[1])/rt;
          double f0 = sin(pi*x[0])*sin(pi*x[1])*pi*pi*(1.0 + delta)*(x[0]*x[0] + x[1]*x[1])
                  + cos(pi*x[0])*sin(pi*x[1])*pi*(1.0 - 3.0*delta)*x[0]
                  + cos(pi*x[1])*sin(pi*x[0])*pi*(1.0 - 3.0*delta)*x[1]
                  + cos(pi*x[1])*cos(pi*x[0])*2.0*pi*pi*(1.0 - delta)*x[0]*x[1];

          return ((f0 + 2.0*(x[0]*(kxx*ux + kxy*uy) + x[1]*(kxy*ux + kyy*uy)))/rt);
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
          double pi = 4.0*atan(1.0);

          return (sin(pi*x[0])*sin(pi*x[1]));
      }

      FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
      {
          FieldVector<RT,n> grad(0);
          double pi = 4.0*atan(1.0);
          grad[0] = pi*cos(pi*x[0])*sin(pi*x[1]);
          grad[1] = pi*cos(pi*x[1])*sin(pi*x[0]);

          return grad;
      }

    private:
        FieldMatrix<DT,n,n> permloc_;
    };
}

#endif
