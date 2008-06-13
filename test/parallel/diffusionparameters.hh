#ifndef DIFFUSIONPARAMETERS_HH
#define DIFFUSIONPARAMETERS_HH

template<class G, class RT>
class DiffusionParameters 
{
  typedef typename G::ctype DT;
  enum {n=G::dimension, m=1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename Dune::IntersectionIteratorGetter<G,Dune::LeafTag>::IntersectionIterator IntersectionIterator;

public:
  DiffusionParameters ()
  {
	for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
		if (i==j)
		  large[i][j] = 1;
		else
		  large[i][j] = 0;
  }

  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return large;
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
	if (x[0] < 1e-6) 
	  return Dune::BoundaryConditions::dirichlet;
	
	return Dune::BoundaryConditions::neumann;
  }

  virtual void dirichletIndex(const Dune::FieldVector<DT,n>& x, const Entity& e,
			      const IntersectionIterator& intersectionIt,
			      const Dune::FieldVector<DT,n>& xi, Dune::FieldVector<int,m>& dirichletIdx) const 
  {
    for (int i = 0; i < m; i++)
      dirichletIdx[i]=i;
    return;
  }

  RT exact(const Dune::FieldVector<DT,n>& x) const
  {
		return (x[0]);//1.0 - x[0] + 2.0*x[1] - 3.0*x[2]);	  
  }
  
  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
			const IntersectionIterator& intersectionIt, 
				const Dune::FieldVector<DT,n>& xi) const
  {
	return (exact(x));
  }
	  
  Dune::FieldVector<RT,n> exactGrad(const Dune::FieldVector<DT,n>& x) const 
  {
	  Dune::FieldVector<DT,n> grad;

	  grad[0] = -1.0; 
	  grad[1] = 0.0; 
	  grad[2] = 0.0;
	  
	  return grad;
  }
  
  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
			const IntersectionIterator& intersectionIt, 
				const Dune::FieldVector<DT,n>& xi) const
  {
	  Dune::FieldVector<RT,n> KGradU(0); 
	  large.umv(exactGrad(x), KGradU);
	  
	  // ASSUMING face-wise constant normal 
	  Dune::FieldVector<DT, n-1> localDimM1(0);
	  return (KGradU*intersectionIt->unitOuterNormal(localDimM1));
  }
	  
private:
  Dune::FieldMatrix<DT,n,n> large;
};

#endif
