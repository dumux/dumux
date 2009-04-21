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
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator
    IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = 8.0*pi*pi*cos(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]) + 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));
        result[1] = -8.0*pi*pi*sin(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]) + 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));
        result[2] = -32.0*pi*pi*sin(4.0*pi*(globalPos[0] + globalPos[1]));

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[1] < 1e-6 || globalPos[1] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        return velocity(globalPos);
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,numEq> result(0);

        FieldVector<Scalar, dim-1> localDimM1(0);
        FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

        FieldVector<Scalar,dim> pN = normal;
        pN *= pressure(globalPos);

        FieldVector<Scalar,dim> muGradVN(0);
        velocityGradient(globalPos).umv(normal, muGradVN);
        muGradVN *= mu(globalPos, element, localPos);

        FieldVector<Scalar,dim> temp = muGradVN;
        temp -= pN;

        for (int i=0; i < dim; ++i)
            result[i] = temp[i];

        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = cos(2.0*pi*globalPos[0])*sin(2.0*pi*globalPos[1]);
        result[1] = -sin(2.0*pi*globalPos[0])*cos(2.0*pi*globalPos[1]);

        return result;
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

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (sin(4.0*pi*(globalPos[0] + globalPos[1])));
    }

    virtual FieldVector<Scalar,dim> pressureGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);

        result[0] = 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));
        result[1] = 4.0*pi*cos(4.0*pi*(globalPos[0] + globalPos[1]));

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
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator
    IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = 5.0*pi*pi*cos(2.0*pi*globalPos[0])*sin(pi*globalPos[1]) + 4.0*pi*globalPos[1]*cos(4.0*pi*globalPos[0]);
        result[1] = -10.0*pi*pi*sin(2.0*pi*globalPos[0])*cos(pi*globalPos[1]) + sin(4.0*pi*(globalPos[0]));
        result[2] = -16.0*pi*pi*globalPos[1]*sin(4.0*pi*globalPos[0]);

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[1] < 1e-6 || globalPos[0] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        return velocity(globalPos);
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,numEq> result(0);

        FieldVector<Scalar, dim-1> localDimM1(0);
        FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

        FieldVector<Scalar,dim> pN = normal;
        pN *= pressure(globalPos);

        FieldVector<Scalar,dim> muGradVN(0);
        velocityGradient(globalPos).umv(normal, muGradVN);
        muGradVN *= mu(globalPos, element, localPos);

        FieldVector<Scalar,dim> temp = muGradVN;
        temp -= pN;

        for (int i=0; i < dim; ++i)
            result[i] = temp[i];

        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = cos(2.0*pi*globalPos[0])*sin(pi*globalPos[1]);
        result[1] = -2.0*sin(2.0*pi*globalPos[0])*cos(pi*globalPos[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (globalPos[1]*sin(4.0*pi*globalPos[0]));
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][0] = -2.0*pi*sin(2.0*pi*globalPos[0])*sin(pi*globalPos[1]);
        result[0][1] = pi*cos(2.0*pi*globalPos[0])*cos(pi*globalPos[1]);
        result[1][0] = -4.0*pi*cos(2.0*pi*globalPos[0])*cos(pi*globalPos[1]);
        result[1][1] = 2.0*pi*sin(2.0*pi*globalPos[0])*sin(pi*globalPos[1]);

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
