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
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

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

    virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1.0e-10)
            return velocity(globalPos);
        else
        {
            FieldVector<Scalar,dim> result(0);
            return result;
        }
    }

    virtual FieldVector<Scalar,dim> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,dim> result(0);

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
            result = muGradVN;
            result -= pN;
        }

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

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> result(0);
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
