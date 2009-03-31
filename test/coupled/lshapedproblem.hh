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
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[1] < 1e-6 || globalPos[1] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;
        else
            return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
      FieldVector<Scalar, dim> vel = velocity(globalPos);
      FieldVector<Scalar, numEq> res(0);
      for (int i = 0; i < dim; ++i) // TODO: pressure
	res[i] = vel[i];
      
      return res;
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    // function returns the square root of the permeability (scalar) divided by alpha
    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        double alpha;
        if (globalPos[1] < 0.5 + 1e-6)
            alpha = 0.1;
        else
            return(-1.0); // realizes outflow bc

        // CHANGE also in the porous medium problem!
        double permeability = 1.0e-2;

        //TODO: divide by viscosity - check?
        return sqrt(permeability)/(alpha*mu(globalPos, element, localPos));
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 0.01;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);

        result[0] = 4.0*globalPos[1]*(1.0 - globalPos[1]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (globalPos[0]);
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
