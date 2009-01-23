// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_YXPROBLEM_HH
#define DUNE_YXPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class G, class RT>
class YXProblem : public StokesProblem<G, RT>
{
    typedef typename G::ctype DT;
    enum {dim=G::dimension, m=G::dimension+1};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<RT,m> q(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> result(0);
        result[0] = x[1];
        result[1] = x[0];

        return result;
    }

    virtual BoundaryConditions::Flags bctype(const FieldVector<DT,dim>& x, const Entity& e,
            const IntersectionIterator& intersectionIt,
            const FieldVector<DT,dim>& xi) const
            {
        if (x[0] > 1 - 1e-6)
            return BoundaryConditions::process;

        return BoundaryConditions::dirichlet;
            }

    virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
            const IntersectionIterator& intersectionIt,
            const FieldVector<DT,dim>& xi) const
            {
        return velocity(x);
            }

    virtual RT Jn(const FieldVector<DT,dim>& x, const Entity& e,
            const IntersectionIterator& intersectionIt,
            const FieldVector<DT,dim>& xi) const
            {
        FieldVector<RT, dim-1> localDimM1(0);
        FieldVector<RT,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
        FieldMatrix<DT, dim, dim> gradU = velocityGradient(x);
        FieldVector<RT,dim> gradUn(0);
        gradU.umv(normal, gradUn);
        RT p = pressure(x);
        RT visc = mu(x, e, xi);

        return (p - visc*(normal*gradUn));
            }

    virtual RT beaversJosephC(const FieldVector<DT,dim>& x, const Entity& e,
            const IntersectionIterator& intersectionIt,
            const FieldVector<DT,dim>& xi) const
            {
        return -1.0;
            }

    virtual FieldVector<RT,dim> Jt(const FieldVector<DT,dim>& x, const Entity& e,
            const IntersectionIterator& intersectionIt,
            const FieldVector<DT,dim>& xi) const
            {
        FieldVector<RT, dim-1> localDimM1(0);
        FieldVector<RT,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
        FieldMatrix<DT, dim, dim> gradU = velocityGradient(x);
        FieldVector<RT,dim> gradUn(0);
        gradU.umv(normal, gradUn);
        RT nGradUn = normal*gradUn;
        FieldVector<RT,dim> nGradUnN(normal);
        nGradUnN *= nGradUn;

        FieldVector<RT,dim> result(gradUn);
        result -= nGradUnN;
        result *= mu(x, e, xi);

        return result;
            }

    virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return 1.0;
    }

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

    YXProblem()
    {}

};

}
#endif
