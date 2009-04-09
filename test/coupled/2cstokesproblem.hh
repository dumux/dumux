#ifndef DUNE_TWOCSTOKESPROBLEM_HH
#define DUNE_TWOCSTOKESPROBLEM_HH

#include"dumux/stokes/stokestransportproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class TwoCStokesProblem : public StokesTransportProblem<Grid, Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=Grid::dimension+2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,numEq> SolutionVector;


    SolutionVector initial (const GlobalPosition& globalPos, const Element& element,
                            const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        result[0] = 1.0e-3;

        if (globalPos[0] > 0.4 && globalPos[0] < 0.6 && globalPos[1] > 0.4 && globalPos[1] < 0.6 )
            result[dim] = 1;

        return result;
    }

    SolutionVector q(const GlobalPosition& globalPos, const Element& element,
                     const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        return result;
    }

    //TODO: implement bcs as FieldVector with numEq entries
    BoundaryConditions::Flags bctype (const GlobalPosition& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const LocalPosition& localPos) const
    {
        if (globalPos[0] < 1e-6 || globalPos[1] < 1e-6 || globalPos[1] > 1 - 1e-6)
            return BoundaryConditions::dirichlet;
        else
            return BoundaryConditions::neumann;
    }

    SolutionVector g(const GlobalPosition& globalPos, const Element& element,
                     const IntersectionIterator& intersectionIt,
                     const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        result = velocity(globalPos, element, localPos);

        return result;
    }

    SolutionVector J(const GlobalPosition& globalPos, const Element& element,
                     const IntersectionIterator& intersectionIt,
                     const LocalPosition& localPos)
    {
        SolutionVector result(0);

        return result;
    }

    // function returns the square root of the permeability (scalar) divided by alpha times mu
    virtual Scalar beaversJosephC(const GlobalPosition& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const LocalPosition& localPos) const
    {
        Scalar alpha;
        if (globalPos[1] < 0.5 + 1e-6)
            alpha = 0.1;
        else
            return(-1.0); // realizes outflow boundary condition

        //TODO: uses only Kxx, extend to permeability tensor
        Scalar permeability = soil().K(globalPos, element, localPos)[0][0];

        return sqrt(permeability)/(alpha*viscosity(globalPos, element, localPos));
    }

    //TODO: call viscosity from material law
    virtual Scalar viscosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.01;
    }

    virtual SolutionVector velocity(const GlobalPosition& globalPos, const Element& element,
                                    const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        result[0] = 4.0*globalPos[1]*(1.0 - globalPos[1]);

        return result;
    }

    FieldMatrix<Scalar,dim,dim> D (const GlobalPosition& globalPos, const Element& element,
                                   const LocalPosition& localPos) const
    {
        FieldMatrix<Scalar,dim,dim> res(0);

        for (int Dx=0; Dx<dim; Dx++)
            for (int Dy=0; Dy<dim; Dy++)
                if (Dx == Dy)
                    res[Dx][Dy] = 1e-5;

        return res;
    }

    Scalar Qg(const GlobalPosition& globalPos, const Element& element,
                      const LocalPosition& localPos) const
    {
        Scalar result = 0;
        return result;
    }

    Fluid& gasPhase () const
    {
        return gasPhase_;
    }

    MultiComp& multicomp () const
    {
        return multicomp_;
    }

    //TODO: gravity vector instead of scalar
    //    FieldVector<Scalar,dim>& gravity() const
    Scalar gravity () const
    {
        return gravity_;
    }

    Matrix2p<Grid, Scalar>& soil() const
    {
        return soil_;
    }
    /*
      FieldMatrix<Scalar, dim, dim> velocityGradient(const GlobalPosition& globalPos) const
      {
      FieldMatrix<Scalar, dim, dim> result(0);

      return result;
      }
    */

    TwoCStokesProblem(Gas_GL& gasPhase, Matrix2p<Grid, Scalar>& soil, MultiComp& multicomp = *(new CWaterAir))
        :
        StokesTransportProblem<Grid,Scalar>(gasPhase, multicomp),
        gasPhase_(gasPhase),
        soil_(soil),
        multicomp_(multicomp)
    {
        gravity_ = 9.81;
        //    for (int i=0; i<dim; ++i)
        //        gravity_[i] = 0;
        //    gravity_[dim] = -9.81;
    }


protected:
    //    FieldVector<Scalar,dim> gravity_;
    Scalar gravity_;
    Gas_GL& gasPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    MultiComp& multicomp_;
};

}
#endif
