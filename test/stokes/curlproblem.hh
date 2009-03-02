// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_CURLPROBLEM_HH
#define DUNE_CURLPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class G, class RT>
class CurlProblem : public StokesProblem<G, RT>
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
        result[0] = -200.0*(2.0*(1.0 - 6.0*x[0] + 6.0*x[0]*x[0])*(x[1] - 3.0*x[1]*x[1] + 2.0*x[1]*x[1]*x[1])
                            + (12.0*x[1] - 6.0)*(x[0]*x[0] - 2.0*x[0]*x[0]*x[0] + x[0]*x[0]*x[0]*x[0]))
            + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
        result[1] = 200.0*(2.0*(1.0 - 6.0*x[1] + 6.0*x[1]*x[1])*(x[0] - 3.0*x[0]*x[0] + 2.0*x[0]*x[0]*x[0])
                           + (12.0*x[0] - 6.0)*(x[1]*x[1] - 2.0*x[1]*x[1]*x[1] + x[1]*x[1]*x[1]*x[1]))
            + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
        result[2] = -32.0*pi*pi*sin(4.0*pi*(x[0] + x[1]));

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

    virtual FieldVector<RT,dim> J(const FieldVector<DT,dim>& x, const Entity& e,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<DT,dim>& xi)
    {
        FieldVector<RT,dim> result(0);
        return result;
    }

    virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
    {
        FieldVector<RT,dim> result(0);
        result[0] = 200.0*((x[1]*(1.0 - x[1])*(1.0 - x[1]) - x[1]*x[1]*(1.0 - x[1]))*x[0]*x[0]*(1.0 - x[0])*(1.0 - x[0]));
        result[1] = -200.0*((x[0]*(1.0 - x[0])*(1.0 - x[0]) - x[0]*x[0]*(1.0 - x[0]))*x[1]*x[1]*(1.0 - x[1])*(1.0 - x[1]));

        return result;
    }

    virtual RT pressure(const FieldVector<DT,dim>& x) const
    {
        return (sin(4.0*pi*(x[0] + x[1])));
    }

    virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
    {
        FieldMatrix<DT, dim, dim> result(0);

        return result;
    }

    CurlProblem()
    {
        pi = 4.0*atan(1.0);
    }

    double pi;
};


}
#endif
