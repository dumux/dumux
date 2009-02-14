// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_YXPROBLEM_HH
#define DUNE_YXPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class YXProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, m=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,m> q(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,m> result(0);
        result[0] = globalPos[1];
        result[1] = globalPos[0];

        return result;
    }

    virtual BoundaryConditions::Flags bctype(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
            {
        if (globalPos[0] > 1 - 1e-6)
            return BoundaryConditions::process;

        return BoundaryConditions::dirichlet;
            }

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
            {
        return velocity(globalPos);
            }

    virtual Scalar Jn(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
            {
        FieldVector<Scalar, dim-1> localDimM1(0);
        FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
        FieldMatrix<Scalar, dim, dim> gradU = velocityGradient(globalPos);
        FieldVector<Scalar,dim> gradUn(0);
        gradU.umv(normal, gradUn);
        Scalar p = pressure(globalPos);
        Scalar visc = mu(globalPos, element, localPos);

        return (p - visc*(normal*gradUn));
            }

    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
            {
        return -1.0;
            }

    virtual FieldVector<Scalar,dim> Jt(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
            {
        FieldVector<Scalar, dim-1> localDimM1(0);
        FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
        FieldMatrix<Scalar, dim, dim> gradU = velocityGradient(globalPos);
        FieldVector<Scalar,dim> gradUn(0);
        gradU.umv(normal, gradUn);
        Scalar nGradUn = normal*gradUn;
        FieldVector<Scalar,dim> nGradUnN(normal);
        nGradUnN *= nGradUn;

        FieldVector<Scalar,dim> result(gradUn);
        result -= nGradUnN;
        result *= mu(globalPos, element, localPos);

        return result;
            }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = -globalPos[1];
        result[1] = -globalPos[0];

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (globalPos[0]*globalPos[1]);
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][1] = -1.0;
        result[1][0] = -1.0;

        return result;
    }

    YXProblem()
    {}

};

}
#endif
