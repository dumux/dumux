// $Id: capillaryflowproblem.hh 733 2009-02-09 08:45:27Z kathinka $

#ifndef DUNE_CAPILLARYSTOKESPROBLEM_HH
#define DUNE_CAPILLARYSTOKESPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class CapillaryStokesProblem : public StokesProblem<Grid, Scalar>
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
        if (globalPos[0] < 1e-10 || globalPos[1] > 8e-6 - 1e-10)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        if (globalPos[0] < 1.0e-10)
        {
            result[0] = 1.5e-3;
            //result = velocity(globalPos);
        }

        return result;
    }

    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& localPos)
    {
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 0.0012;
    }

    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        double alpha;
        if (globalPos[1] > 1e-10 )
            return -1.0;
        else
            return 0;//alpha = 0.1;

        // CHANGE also in the porous medium problem!
        double permeability;
        if (globalPos[0] > 0.2e-3 && globalPos[0] < 0.8e-3)
            permeability = 5e-15*mu(globalPos, element, localPos);
        else
            permeability = 1e-20*mu(globalPos, element, localPos);

        //TODO: divide by viscosity - check?
        return sqrt(permeability)/(alpha*mu(globalPos, element, localPos));
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] = -8.3125e8 * globalPos[1] * globalPos[1] + 6650.0 * globalPos[1]; //entspricht v_m = 1 mm/s
        result[1] = 0.0;


        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (-1995000.0*globalPos[0] + 266.0);
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][1] = -2.0*8.3125e8 * globalPos[1] + 6650.0;

        return result;
    }

    CapillaryStokesProblem()
    {}
};

}
#endif
