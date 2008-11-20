#ifndef DIFFUSIONPARAMETERS_HH
#define DIFFUSIONPARAMETERS_HH
#include <dumux/material/properties.hh>

template<class G, class RT>
class DiffusionParameters
{
  typedef typename G::ctype DT;
  enum {n=G::dimension, m=1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename Dune::IntersectionIteratorGetter<G,Dune::LeafTag>::IntersectionIterator IntersectionIterator;

public:
  DiffusionParameters (Medium& law = *(new Uniform)): materialLaw_(law)
  {
	for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
		if (i==j)
		  large[i][j] = 1e-10;
		else
		  large[i][j] = 0;

    gravity_[0] = 0;
    gravity_[1] = 0;
    gravity_[n-1] = -9.81;

	Temperature = 288.15;
	density_= materialLaw_.density(Temperature);
  }

  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return large;
  }

  const Dune::FieldVector<RT,n>& gravity () const
  	{
  		return gravity_;
  	}

  const RT porosity (const Dune::FieldVector<DT,n>& x, const Entity& e,
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return 0.5;
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
	    Dune::FieldVector<Dune::BoundaryConditions::Flags, 1> values; //(Dune::BoundaryConditions::neumann);

//		std::cout << "global coordinate "<< x << "bctype: boundaryId = " << intersectionIt.boundaryId() << std::endl;

		switch (intersectionIt->boundaryId()) {
                    case 1:
                    case 2:
                    case 5:
                    case 6:
			values = Dune::BoundaryConditions::neumann;
			break;
                    
                    case 3:
                    case 4:
                    default:
			values = Dune::BoundaryConditions::dirichlet;
			break;
		}

		return values;
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
	  return (x[0]);//1.0 - x[0] + 2.0*x[1] - 3.0*x[2]);
  }

  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
				const Dune::FieldVector<DT,n>& xi) const
  {
		RT values(0);

		switch (intersectionIt->boundaryId()) {
		case 3:
			values= 1.001e+5;
			break;
		case 4:
			values= 1.0e+5;
			break;
		}

		return values;
  }

  RT initial   (const Dune::FieldVector<DT,n>& x, const Entity& e,
				  const Dune::FieldVector<DT,n>& xi) const
  {
	return 1.0e+5;
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
	  RT values(0);

	  switch (intersectionIt->boundaryId()) {
	  case 1: case 2: case 5: case 6:
		  values = 0;
		  break;

	  }

	  return values;
		  }

  Medium& materialLaw() {
	  return materialLaw_;
  }
  const RT& Temp(){
	  return Temperature;
  }
private:
  Dune::FieldMatrix<DT,n,n> large;
  RT density_;
  Dune::FieldVector<RT,n> gravity_;
  RT Temperature;
  Medium& materialLaw_;

};

#endif
