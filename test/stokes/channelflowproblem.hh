// $Id: channelflowproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_CHANNELFLOWPROBLEM_HH
#define DUNE_CHANNELFLOWPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class G, class RT>
class ChannelFlowProblem : public StokesProblem<G, RT>
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
        if (x[0] > 8 - 1e-6)
            return BoundaryConditions::neumann;

        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<DT,dim>& xi) const
    {
        if (x[0] < 1e-6)
            return velocity(x);
        else
        {
            FieldVector<RT,dim> result(0);
            return result;
        }
    }

    virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
    {
        FieldVector<RT,dim> result(0);
        result[0] = 4.0*(x[1] - x[1]*x[1]);
        result[1] = 0.0;

        return result;
    }

    virtual RT pressure(const FieldVector<DT,dim>& x) const
    {
        return (-8.0*x[0]);
    }

    virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
    {
        FieldMatrix<DT, dim, dim> result(0);
        result[0][1] = 4.0*(1.0 - 2.0*x[1]);

        return result;
    }

    ChannelFlowProblem()
    {}

};

}
#endif
