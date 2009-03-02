// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_SINPROBLEM_HH
#define DUNE_SINPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class SinProblem : public StokesProblem<Grid, Scalar>
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension, m=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,m> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                    const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,m> result(0);
        result[0] = 8.0*pi*pi*cos(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]) + 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));
        result[1] = -8.0*pi*pi*sin(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]) + 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        return velocity(globalPos);
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = cos(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]);
        result[1] = -sin(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (sin(4.0*pi*(globalPos[0] + globalPos[1])));
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][0] = -2.0*pi*sin(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]);
        result[0][1] = 2.0*pi*cos(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]);
        result[1][0] = -2.0*pi*cos(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]);
        result[1][1] = 2.0*pi*sin(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]);

        return result;
    }

    SinProblem()
    {
        pi = 4.0*atan(1.0);
    }

    double pi;
};

template<class Grid, class Scalar>
class SinProblem2 : public StokesProblem<Grid, Scalar>
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension, m=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,m> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                    const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,m> result(0);
        result[0] = 5.0*pi*pi*cos(2.0*pi*globalPos[0])*sin(pi*globalPos[1]) + 4.0*pi*cos(4.0*pi*(globalPos[0])*globalPos[1]);
        result[1] = -10.0*pi*pi*sin(2.0*pi*globalPos[0])*cos(pi*globalPos[1]) + sin(4.0*pi*(globalPos[0]));

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        return velocity(globalPos);
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = cos(2.0*pi*globalPos[0])*sin(pi*globalPos[1]);
        result[1] = -2.0*sin(2.0*pi*globalPos[0])*cos(pi*globalPos[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (sin(4.0*pi*(globalPos[0]))*globalPos[1]);
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
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
