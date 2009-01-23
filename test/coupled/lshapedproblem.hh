// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_LSHAPEDPROBLEM_HH
#define DUNE_LSHAPEDPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class G, class RT>
class LShapedProblem : public StokesProblem<G, RT>
{
  typedef typename G::ctype DT;
  enum {dim=G::dimension, m=G::dimension+1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

public:
  virtual FieldVector<RT,m> q(const FieldVector<DT,dim>& x, const Entity& e,
				const FieldVector<DT,dim>& xi) const
  {
    FieldVector<RT,m> result(0);

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
					    const IntersectionIterator& intersectionIt,
					    const FieldVector<DT,dim>& xi) const
  {
    if (x[0] > 1.5 - 1e-6 && x[1] < 0.5 + 1e-6)
      return BoundaryConditions::process;
   	if (x[0] > 4 - 1e-6)
      return BoundaryConditions::neumann;
    else
      return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				const FieldVector<DT,dim>& xi) const
  {
    return velocity(x);
  }

  // function returns the square root of the permeability (scalar) divided by alpha
  virtual RT beaversJosephC(const FieldVector<DT,dim>& x, const Entity& e,
		const IntersectionIterator& intersectionIt,
		const FieldVector<DT,dim>& xi) const
  {
      // CHANGE also in the porous medium problem!
	  double permeability = 1.0e-2;//5.88e-5;
	  double alpha = 1.0;

	  //TODO: divide by viscosity - check?
	  return sqrt(permeability)/alpha;
  }

  //TODO: this is not called yet
  virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  {
    return 0.01;
  }

  virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
  {
	  FieldVector<RT,dim> result(0);

	  if (x[0] < 1e-6)
		  result[0] = 1;//-x[1];
	  if (x[1] < 1e-6)
		  result[0] = 0;//-x[1];
	  if (x[1] > 1-1e-6)
		  result[0] = 0;//-x[1];

    return result;
  }

  virtual RT pressure(const FieldVector<DT,dim>& x) const
  {
    return (x[0]);
  }

  virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
  {
    FieldMatrix<DT, dim, dim> result(0);
//    result[0][1] = -1.0;
//    result[1][0] = -1.0;
    result[0][0] = -1.0;

    return result;
  }

  LShapedProblem()
  {}

};

}
#endif
