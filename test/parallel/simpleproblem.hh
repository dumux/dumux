#ifndef SIMPLEPROBLEM_HH
#define SIMPLEPROBLEM_HH

template<class G, class RT>
class ParallelModelProblem : public Dune::GroundwaterEquationParameters<G,RT>
{
  typedef typename G::ctype DT;
  enum {n=G::dimension};
  typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
  ParallelModelProblem ()
  {
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        if (i==j)
          large[i][j] = 1;
        else
          large[i][j] = 0;
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

  virtual Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
                       const Dune::FieldVector<DT,n>& xi) const
  {
    if (x[0]<0.5+1E-6)
      return Dune::BoundaryConditions::dirichlet;
    return Dune::BoundaryConditions::neumann;
  }

  virtual RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    return 1.0;
  }

  virtual RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
                const Dune::FieldVector<DT,n>& xi) const
  {
    return 1.0;
  }

private:
  Dune::FieldMatrix<DT,n,n> large;
};

#endif
