// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_YXPROBLEM_HH
#define DUNE_YXPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class YXProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
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

    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const FieldVector<Scalar,dim>& localPos) const
            {
        FieldVector<Scalar,numEq> result(0);
        result[0] = globalPos[1];
        result[1] = globalPos[0];

        return result;
            }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 1.0;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < eps_)// || globalPos[1] < eps_)// || globalPos[1] > 1 - eps_)
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

        // ASSUMING face-wise constant normal
        FieldVector<Scalar, dim-1> localDimM1(0);
        FieldVector<Scalar,dim> normal = intersectionIt->unitOuterNormal(localDimM1);

        FieldVector<Scalar,dim> pN = normal;
        pN *= pressure(globalPos);

        FieldVector<Scalar,dim> muGradVN(0);
        velocityGradient(globalPos).umv(normal, muGradVN);
        muGradVN *= mu(globalPos, element, localPos);

        Scalar muGradVNN = muGradVN*normal;

        result = normal;
        result *= muGradVNN;
        result -= pN;

        return result;
    }

    //! evaluate Beavers-Joseph proportionality constant at given position
    /*! evaluate Beavers-Joseph proportionality constant \f$c = \sqrt(k)/\alpha\f$
    such that \f$u_\tau = - c (\mu \nabla u\cdot n)_\tau\f$
    @param[in]  globalPos    position in global coordinates
    \return     value of the proportionality constant
    */
    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<Scalar,dim>& localPos) const
    {
        return (1.0/mu(globalPos, element, localPos));
    }

    YXProblem()
    {
        eps_= 1e-6;
    }

private:
    Scalar eps_;

};

}
#endif
