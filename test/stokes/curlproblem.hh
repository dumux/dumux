// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_CURLPROBLEM_HH
#define DUNE_CURLPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class CurlProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator
        IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Entity& e,
                                    const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = -200.0*(2.0*(1.0 - 6.0*x[0] + 6.0*x[0]*x[0])*(x[1] - 3.0*x[1]*x[1] + 2.0*x[1]*x[1]*x[1])
                            + (12.0*x[1] - 6.0)*(x[0]*x[0] - 2.0*x[0]*x[0]*x[0] + x[0]*x[0]*x[0]*x[0]))
            + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
        result[1] = 200.0*(2.0*(1.0 - 6.0*x[1] + 6.0*x[1]*x[1])*(x[0] - 3.0*x[0]*x[0] + 2.0*x[0]*x[0]*x[0])
                           + (12.0*x[0] - 6.0)*(x[1]*x[1] - 2.0*x[1]*x[1]*x[1] + x[1]*x[1]*x[1]*x[1]))
            + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
        result[2] = -32.0*pi*pi*sin(4.0*pi*(x[0] + x[1]));

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Entity& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const
    {
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& x, const Entity& e,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& xi) const
    {
        return velocity(x);
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& x, const Entity& e,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& xi)
    {
        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = 200.0*((x[1]*(1.0 - x[1])*(1.0 - x[1]) - x[1]*x[1]*(1.0 - x[1]))*x[0]*x[0]*(1.0 - x[0])*(1.0 - x[0]));
        result[1] = -200.0*((x[0]*(1.0 - x[0])*(1.0 - x[0]) - x[0]*x[0]*(1.0 - x[0]))*x[1]*x[1]*(1.0 - x[1])*(1.0 - x[1]));

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
        return (sin(4.0*pi*(x[0] + x[1])));
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);

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
