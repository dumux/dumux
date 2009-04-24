// $Id: exampleproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_EXAMPLEPROBLEM_HH
#define DUNE_EXAMPLEPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class ExampleProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,dim> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = 1.0;
        result[1] = 2.0;

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        //     if (globalPos[1] < 1e-6)
        //       return BoundaryConditions::process;

        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,dim>dirichlet(const FieldVector<Scalar,dim>& globalPos, const Element& element,
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
        if (dim != 2) {
            DUNE_THROW(NotImplemented, "analytic Beavers-Joseph for dim != 2");
        }
        else {
            FieldVector<Scalar, dim-1> localDimM1(0);
            FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
            FieldVector<Scalar,dim> gradUnT = Jt(globalPos, element, intersectionIt, localPos);
            Scalar tGradUnT = gradUnT[0]*normal[1] - gradUnT[1]*normal[0];
            if (fabs(tGradUnT) < 1e-12)
                DUNE_THROW(MathError, "analytic Beavers-Joseph does not work for tGradUnT == 0");

            FieldVector<Scalar,dim> u =dirichlet(globalPos, element, intersectionIt, localPos);
            Scalar uT = u[0]*normal[1] - u[1]*normal[0];

            Scalar c = -uT/tGradUnT;
            if (fabs(c) < 1e-12)
                DUNE_THROW(MathError, "analytic Beavers-Joseph gives c == 0");

            return -uT/tGradUnT;
        }
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
        result[0] = 2.0*globalPos[0]*globalPos[1];
        result[1] = -globalPos[1]*globalPos[1];

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return globalPos[0] - 1.0;
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][0] = 2.0*globalPos[1];
        result[0][1] = 2.0*globalPos[0];
        result[1][1] = -2.0*globalPos[1];

        return result;
    }

    ExampleProblem()
    {}

};

}
#endif
