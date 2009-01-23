#ifndef DARCYPARAMETERS_HH
#define DARCYPARAMETERS_HH

template<class G, class RT>
class DarcyParameters
{
  typedef typename G::ctype DT;
  enum {n=G::dimension, m=1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename Dune::IntersectionIteratorGetter<G,Dune::LeafTag>::IntersectionIterator IntersectionIterator;

public:
  DarcyParameters ()
  {
      // CHANGE also in the Stokes problem!
	for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
		if (i==j)
		  permeability_[i][j] = 1.0e-2;//5.88e-5;
		else
			permeability_[i][j] = 0;
  }

  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return permeability_;
  }

  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return 0;
  }

  Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
					   const Dune::FieldVector<DT,n>& xi) const
  {
 	if (x[0] > 4 - 1e-6)
	  return Dune::BoundaryConditions::dirichlet;

	return Dune::BoundaryConditions::neumann;
  }

	virtual void dirichletIndex(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<DT,n>& xi, Dune::FieldVector<int,m>& dirichletIndex) const
	{
		for (int i = 0; i < m; i++)
			dirichletIndex[i]=i;
		return;
	}

  RT exact(const Dune::FieldVector<DT,n>& x) const
  {
		return (x[0]*x[1]);
  }

  // Dirichlet boundary conditions
  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
				const Dune::FieldVector<DT,n>& xi) const
  {
	return 0;
  }

  Dune::FieldVector<RT,n> exactGrad(const Dune::FieldVector<DT,n>& x) const
  {
	  Dune::FieldVector<DT,n> grad;

	  grad[0] = x[1];
	  grad[1] = x[0];
	  grad[2] = 0.0;

	  return grad;
  }

  // Neumann b.c.
  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
				const Dune::FieldVector<DT,n>& xi) const
  {
	  return 0;
//	  if (x[0] < 0.5 + 1e-6 && x[1] > 0.4 - 1e-6)
//		  return 0;
//	  else
//	  {
//		  Dune::FieldVector<RT,n> KGradU(0);
//		  large.umv(exactGrad(x), KGradU);
//
//		  // ASSUMING face-wise constant normal
//		  Dune::FieldVector<DT, n-1> localDimM1(0);
//		  return -(KGradU*intersectionIt->unitOuterNormal(localDimM1));
//	  }
  }

private:
  Dune::FieldMatrix<DT,n,n> permeability_;
};

#endif
