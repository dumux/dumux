// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_SINPROBLEM_HH
#define DUNE_SINPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class G, class RT>
class SinProblem : public StokesProblem<G, RT>
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
    result[0] = 8.0*pi*pi*cos(2.0*pi*x[0])*sin(2.0*pi*x[1]) + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
    result[1] = -8.0*pi*pi*sin(2.0*pi*x[0])*cos(2.0*pi*x[1]) + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<DT,dim>& xi) const
  {
    return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi) const
  {
    return velocity(x);
  }

  virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  {
    return 1.0;
  }

  virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
  {
    FieldVector<RT,dim> result(0);
    result[0] = cos(2.0*pi*x[0])*sin(2.0*pi*x[1]);
    result[1] = -sin(2.0*pi*x[0])*cos(2.0*pi*x[1]);

    return result;
  }

  virtual RT pressure(const FieldVector<DT,dim>& x) const
  {
    return (sin(4.0*pi*(x[0] + x[1])));
  }

  virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
  {
    FieldMatrix<DT, dim, dim> result(0);
    result[0][0] = -2.0*pi*sin(2.0*pi*x[0])*sin(2.0*pi*x[1]);
    result[0][1] = 2.0*pi*cos(2.0*pi*x[0])*cos(2.0*pi*x[1]);
    result[1][0] = -2.0*pi*cos(2.0*pi*x[0])*cos(2.0*pi*x[1]);
    result[1][1] = 2.0*pi*sin(2.0*pi*x[0])*sin(2.0*pi*x[1]);

    return result;
  }

  SinProblem()
  {
      pi = 4.0*atan(1.0);
  }

  double pi;
};

template<class G, class RT>
class SinProblem2 : public StokesProblem<G, RT>
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
    result[0] = 5.0*pi*pi*cos(2.0*pi*x[0])*sin(pi*x[1]) + 4.0*pi*cos(4.0*pi*(x[0])*x[1]);
    result[1] = -10.0*pi*pi*sin(2.0*pi*x[0])*cos(pi*x[1]) + sin(4.0*pi*(x[0]));

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<DT,dim>& xi) const
  {
    return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<DT,dim>& xi) const
  {
    return velocity(x);
  }

  virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  {
    return 1.0;
  }

  virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
  {
    FieldVector<RT,dim> result(0);
    result[0] = cos(2.0*pi*x[0])*sin(pi*x[1]);
    result[1] = -2.0*sin(2.0*pi*x[0])*cos(pi*x[1]);

    return result;
  }

  virtual RT pressure(const FieldVector<DT,dim>& x) const
  {
    return (sin(4.0*pi*(x[0]))*x[1]);
  }

  virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
  {
    FieldMatrix<DT, dim, dim> result(0);

    return result;
  }

  SinProblem2()
  {
      pi = 4.0*atan(1.0);
  }

  double pi;
};

}
#endif
