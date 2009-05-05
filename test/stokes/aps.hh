// $Id: channelflowproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_APS_HH
#define DUNE_APS_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class APS : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator
    IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = pi/8.0*cos(2*pi*x[1])*sin(2*pi*x[0]) +
                    pi*pi*(1 - 2*cos(2*pi*x[0]))*sin(2*pi*x[1]);

        result[1] = pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1]) -
                    pi*pi*(1 - 2*cos(2*pi*x[1]))*sin(2*pi*x[0]);

        result[2] = 0.5*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1]);

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[1] < 1e-6 || globalPos[1] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;
        return BoundaryConditions::neumann;

//    		return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,numEq> dirichlet(const FieldVector<Scalar,dim>& x, const Element& element,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
        result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);

        return result;
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& x, const Element& element,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& localPos)
    {
    	FieldVector<Scalar,numEq> result(0);

    	FieldVector<Scalar, dim-1> localDimM1(0);
    	FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

        FieldVector<Scalar,dim> pN = normal;
        pN *= pressure(x);

        FieldVector<Scalar,dim> muGradVN(0);
        velocityGradient(x).umv(normal, muGradVN);
        muGradVN *= mu(x, element, localPos);

        FieldVector<Scalar,dim> temp = muGradVN;
        temp -= pN;

        for (int i=0; i < dim; ++i)
            result[i] = temp[i];

        return result;
    }

   virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        return 0;//(-1.0);
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
        result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
    	Scalar result = -cos(2*pi*x[0])*cos(2*pi*x[1])/16.0;
    	return result;
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);

        result[0][0] = 1/2.0*pi*sin(2*pi*x[0])*sin(2*pi*x[1]);
        result[0][1] = pi*cos(2*pi*x[1])*sin(pi*x[0])*sin(pi*x[0]);
        result[1][0] = -pi*cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
        result[1][1] = -1/2.0*pi*sin(2*pi*x[0])*sin(2*pi*x[1]);

        return result;
    }

    APS()
    {
    	pi = 4.0*atan(1.0);
    }

protected:
    double pi;
};

}
#endif
