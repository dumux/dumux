// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_SINPROBLEM_HH
#define DUNE_SINPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class SinProblem : public StokesProblem<Grid, Scalar>
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
        result[0] = 8.0*pi*pi*cos(2.0*pi*x[0])*sin(2.0*pi*x[1]) + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
        result[1] = -8.0*pi*pi*sin(2.0*pi*x[0])*cos(2.0*pi*x[1]) + 4.0*pi*cos(4.0*pi*(x[0] + x[1]));
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

    virtual FieldVector<Scalar,dim> J(const FieldVector<Scalar,dim>& x, const Entity& e,
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

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = cos(2.0*pi*x[0])*sin(2.0*pi*x[1]);
        result[1] = -sin(2.0*pi*x[0])*cos(2.0*pi*x[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
        return (sin(4.0*pi*(x[0] + x[1])));
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
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

/** \todo Please doc me! */

template<class Grid, class Scalar>
class SinProblem2 : public StokesProblem<Grid, Scalar>
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
        result[0] = 5.0*pi*pi*cos(2.0*pi*x[0])*sin(pi*x[1]) + 4.0*pi*x[1]*cos(4.0*pi*x[0]);
        result[1] = -10.0*pi*pi*sin(2.0*pi*x[0])*cos(pi*x[1]) + sin(4.0*pi*(x[0]));
        result[2] = -16.0*pi*pi*x[1]*sin(4.0*pi*x[0]);

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Entity& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const
    {
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& x, const Entity& e,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& xi) const
    {
        return velocity(x);
    }

    virtual FieldVector<Scalar,dim> J(const FieldVector<Scalar,dim>& x, const Entity& e,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& xi)
    {
        FieldVector<Scalar,dim> result(0);
        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = cos(2.0*pi*x[0])*sin(pi*x[1]);
        result[1] = -2.0*sin(2.0*pi*x[0])*cos(pi*x[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
        return (x[1]*sin(4.0*pi*x[0]));
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);

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
