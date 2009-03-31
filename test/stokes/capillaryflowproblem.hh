// $Id: capillaryflowproblem.hh 733 2009-02-09 08:45:27Z kathinka $

#ifndef DUNE_CAPILLARYFLOWPROBLEM_HH
#define DUNE_CAPILLARYFLOWPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar>
class CapillaryFlowProblem : public StokesProblem<Grid, Scalar>
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

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-10 || globalPos[1] > 0.00001 - 1e-10)// || globalPos[1] < 1e-10)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1.0e-10)
            return velocity(globalPos);
        else
        {
            FieldVector<Scalar,numEq> result(0);
            return result;
        }
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,numEq> result(0);
        FieldVector<Scalar,dim> temp(0);	

        if (globalPos[1] < 1e-10)
        {
            // ASSUMING face-wise constant normal
            FieldVector<Scalar, dim-1> localDimM1(0);
            FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

            FieldVector<Scalar,dim> pN = normal;
            pN *= pressure(globalPos);

            FieldVector<Scalar,dim> muGradVN(0);
            velocityGradient(globalPos).umv(normal, muGradVN);
            muGradVN *= mu(globalPos, element, localPos);

            Scalar muGradVNN = muGradVN*normal;

            //result = normal;
            temp = muGradVN;
            temp -= pN;
            
        }

        // this is some workaround for the boxjacobian
        // and the instantiation of the function assembleBoundaryCondition
        for (int i=0; i < dim; ++i)
            result[i] = temp[i];
            
        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 0.016625;
    }

    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                              const IntersectionIterator& intersectionIt,
                              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[1] < 1e-10)
            return 0;
        else
            return -1.0;
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = -60000000.0 * globalPos[1] * globalPos[1] + 600.0 * globalPos[1]; //entspricht v_m = 1 mm/s
        result[1] = 0.0;


        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (-1995000.0*globalPos[0]);
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][1] = -120000000.0 * globalPos[1] + 600.0;

        return result;
    }

    CapillaryFlowProblem()
    {}
};

}
#endif
