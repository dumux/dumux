#ifndef TESTPROBLEM_HH
#define TESTPROBLEM_HH


/** \todo Please doc me! */

template<class G, class RT>
class TestModelProblem : public Dune::GroundwaterEquationParameters<G,RT>
{
  typedef typename G::ctype DT;
  enum {n=G::dimension};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
public:
  TestModelProblem ()
  {
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        if (i==j)
          large[i][j] = 1;
        else
          large[i][j] = 0;
    large[0][0] = 1;
  }

  virtual const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
                  const Dune::FieldVector<DT,n>& xi) const
  {
    return large;
  }

  virtual RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
                  const Dune::FieldVector<DT,n>& xi) const
  {
    return 0;
  }

  virtual typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
                       const Dune::FieldVector<DT,n>& xi) const
  {
      if (x[0] < 1e-6 || x[0] > 599.9999)
          return Dune::BoundaryConditions::dirichlet;
      else
          return Dune::BoundaryConditions::neumann;
  }

  virtual RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    return (1e6*(1 - 1.0/600.0*x[0] ));
  }

  virtual RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    return 0;
  }

private:
  Dune::FieldMatrix<DT,n,n> large;
};

/** \todo Please doc me! */

template<class G, typename ct, int dim>
class LinearFunction : virtual public Dune::C0GridFunction<G,double,1>,
                  private Dune::FunctionDefault<typename G::ctype,double,G::dimension,1>,
                  private Dune::DifferentiableFunctionDefault<typename G::ctype,double,G::dimension,1>,
                  private Dune::GridFunctionDefault<G,double,1>,
                  private Dune::DifferentiableGridFunctionDefault<G,double,1>
{
public:
  LinearFunction () {}
  double operator() (const Dune::FieldVector<ct,dim>& x) const
  {
    if (dim == 3)
      return (x[0]+x[1]+x[2])/3.0;
    else
      return 0.5*(x[0]+x[1]);
  }
  Dune::FieldVector<ct,dim> grad (const Dune::FieldVector<ct,dim>& x) const
  {
    Dune::FieldVector<ct,dim> g(1.0/dim);
    return g;
  }
  double laplace (const Dune::FieldVector<ct,dim>& x) const
  {
    return 0;
  }
  double eval (int comp, const Dune::FieldVector<ct,dim>& x) const
  {
    return (*this)(x);
  }
  int order () const
  {
    return Dune::InfinitelyDifferentiable;
  }
};

#endif
