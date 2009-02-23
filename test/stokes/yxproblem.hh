// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_YXPROBLEM_HH
#define DUNE_YXPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class G, class RT>
class YXProblem : public StokesProblem<G, RT>
{
  typedef typename G::ctype DT;
  enum {dim=G::dimension, numEq=G::dimension+1};
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
    {
      FieldVector<RT,dim> result(0);
      result[0] = -x[1];
      result[1] = -x[0];

      return result;
    }

    virtual RT pressure(const FieldVector<DT,dim>& x) const
    {
      return (x[0]*x[1]);
    }

    virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
    {
      FieldMatrix<DT, dim, dim> result(0);
      result[0][1] = -1.0;
      result[1][0] = -1.0;

      return result;
    }

    virtual FieldVector<RT,numEq> q(const FieldVector<DT,dim>& x, const Entity& e,
                const FieldVector<DT,dim>& xi) const
  {
    FieldVector<RT,numEq> result(0);
    result[0] = x[1];
    result[1] = x[0];

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<DT,dim>& xi) const
  {
     if (x[0]< 1e-6)
       return BoundaryConditions::dirichlet;

    return BoundaryConditions::neumann;
  }

  virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi) const
  {
    return velocity(x);
  }

  virtual FieldVector<RT,dim> J(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi)
  {
      FieldVector<RT,dim> result(0);

      // ASSUMING face-wise constant normal
      FieldVector<RT, dim-1> localDimM1(0);
      FieldVector<RT,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

      FieldVector<RT,dim> pN = normal;
      pN *= pressure(x);

      FieldVector<RT,dim> muGradVN(0);
      velocityGradient(x).umv(normal, muGradVN);
      muGradVN *= mu(x, e, xi);

      result = muGradVN;
      result -= pN;

      return result;
  }

  virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  {
    return 1.0;
  }

  YXProblem()
  {}

};

}
#endif
