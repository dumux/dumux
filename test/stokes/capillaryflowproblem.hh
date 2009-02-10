// $Id: capillaryflowproblem.hh 733 2009-02-09 08:45:27Z kathinka $

#ifndef DUNE_CAPILLARYFLOWPROBLEM_HH
#define DUNE_CAPILLARYFLOWPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class G, class RT>
class CapillaryFlowProblem : public StokesProblem<G, RT>
{
  typedef typename G::ctype DT;
  enum {dim=G::dimension, numEq=G::dimension+1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

public:
  virtual FieldVector<RT,numEq> q(const FieldVector<DT,dim>& x, const Entity& e,
                const FieldVector<DT,dim>& xi) const
  {
    FieldVector<RT,numEq> result(0);

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<DT,dim>& xi) const
  {
	  if (x[0] > 0.0009999)
		  return BoundaryConditions::neumann;

	  return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi) const
  {
	  if (x[0] < 1.0e-8)
		  return velocity(x);
	  else
	  {
		  FieldVector<RT,dim> result(0);
		  return result;
	  }
  }

  virtual FieldVector<RT,dim> J(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi)
  {
      FieldVector<RT,dim> result(0);
      return result;
  }

  virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  {
    return 0.016625;
  }

  virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
  {
    FieldVector<RT,dim> result(0);
    result[0] = -60000000.0 * x[1] * x[1] + 600.0 * x[1]; //entspricht v_m = 1 mm/s
    result[1] = 0.0;


    return result;
  }

  virtual RT pressure(const FieldVector<DT,dim>& x) const
  {
    return (-1995000.0*x[0]);
  }

  virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
  {
    FieldMatrix<DT, dim, dim> result(0);
    result[0][1] = -120000000.0 * x[1] + 600.0;

    return result;
  }

  CapillaryFlowProblem()
  {}
};

}
#endif
