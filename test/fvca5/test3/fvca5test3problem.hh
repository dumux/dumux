#ifndef FVCA5TEST3PROBLEM_HH
#define FVCA5TEST3PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
    template<class G, class RT, class VC>
    class FVCA5Test3Problem : public DeprecatedDiffusionProblem<G,RT,VC>
    {
      typedef typename G::ctype DT;
      enum {n=G::dimension};
      typedef typename G::Traits::template Codim<0>::Entity Entity;

    public:
      FVCA5Test3Problem(VC& variables, double delta = 1.0e-3, double theta = 0.6981317007977316802)
        : DeprecatedDiffusionProblem<G,RT,VC>(variables)
      {
          double cost = cos(theta);
          double sint = sqrt(1.0 - cost*cost);

          permloc_[0][0] = cost*cost + delta*sint*sint;
          permloc_[1][1] = sint*sint + delta*cost*cost;
          permloc_[0][1] = permloc_[1][0] = cost*sint*(1.0 - delta);
      }

      FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
                      const FieldVector<DT,n>& xi)
      {
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
          double u;

          if (fabs(x[1]) < 1e-6) {
              if (x[0] < 0.2)
                  u = 1.0;
              else if (x[0] < 0.3)
                  u = -5.0*x[0] + 2.0;
              else
                  u = 0.5;
          }
          else if (fabs(1.0 - x[1]) < 1e-6) {
              if (x[0] < 0.7)
                  u = 0.5;
              else if (x[0] < 0.8)
                  u = -5.0*x[0] + 4.0;
              else
                  u = 0.0;
          }
          else if (fabs(x[0]) < 1e-6) {
              if (x[1] < 0.2)
                  u = 1.0;
              else if (x[1] < 0.3)
                  u = -5.0*x[1] + 2.0;
              else
                  u = 0.5;
          }
          else {
              if (x[1] < 0.7)
                  u = 0.5;
              else if (x[1] < 0.8)
                  u = -5.0*x[1] + 4.0;
              else
                  u = 0.0;
          }

          return (u);
      }

      RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const
      {
        return 0;
      }

      RT exact (const FieldVector<DT,n>& x) const
      {
          return 0;
      }

      FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
      {
          FieldVector<RT,n> grad(0);
          return grad;
      }

    private:
        FieldMatrix<DT,n,n> permloc_;
    };
}

#endif
