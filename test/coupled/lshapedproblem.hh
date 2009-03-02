// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_LSHAPEDPROBLEM_HH
#define DUNE_LSHAPEDPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class LShapedProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = globalPos[1];
        result[1] = globalPos[0];

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[0] > 4 - 1e-6 || globalPos[1] < 1e-6 || globalPos[1] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;
        else
            return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        return velocity(globalPos);
    }

    virtual FieldVector<Scalar,dim> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,dim> result(0);

        //      FieldVector<Scalar, dim-1> localDimM1(0);
        //      FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);
        //
        //      FieldVector<Scalar,dim> pN = normal;
        //      pN *= pressure(globalPos);
        //
        //      FieldVector<Scalar,dim> muGradVN(0);
        //      velocityGradient(globalPos).umv(normal, muGradVN);
        //      muGradVN *= mu(globalPos, element, localPos);
        //
        //      Scalar muGradVNN = muGradVN*normal;
        //
        //      result = normal;
        //      result *= muGradVNN;
        //      result -= pN;

        return result;
    }

    // function returns the square root of the permeability (scalar) divided by alpha
    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        //      if (globalPos[0] > 4 - 1e-6)
        //          return 0;
        //      else
        //      {
        // CHANGE also in the porous medium problem!
        double permeability = 1.0;
        double alpha;
        if (globalPos[0] < 1.5 + 1e-6)
            alpha = 2.0/3.0;
        else
            alpha = -2.0;

        //TODO: divide by viscosity - check?
        return sqrt(permeability)/(alpha*mu(globalPos, element, localPos));
        //      }
    }

    //TODO: this is not called yet
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

    LShapedProblem()
    {}

};

}
#endif
