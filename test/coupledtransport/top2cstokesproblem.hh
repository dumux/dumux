#ifndef DUNE_TOPTWOCSTOKESPROBLEM_HH
#define DUNE_TOPTWOCSTOKESPROBLEM_HH

#include"dumux/stokes/stokestransportproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class TopTwoCStokesProblem : public StokesTransportProblem<Grid, Scalar>
{
    enum {velocityXIdx=0, velocityYIdx=1, massFracIdx=2, pressureIdx=3};
    enum {dim=Grid::dimension, numEq=Grid::dimension+2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,numEq> SolutionVector;

public:
    SolutionVector initial (const GlobalPosition& globalPos, const Element& element,
                            const LocalPosition& localPos) const
    {

        SolutionVector result(0);

        result[velocityXIdx] = 0;
        result[velocityYIdx] = -1;
        result[massFracIdx] = 1e-6;
//        result[pressureIdx] = 1e5;

        return result;
    }

    SolutionVector q(const GlobalPosition& globalPos, const Element& element,
                     const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        return result;
    }

    FieldVector<BoundaryConditions::Flags, numEq> bctype (const GlobalPosition& globalPos, const Element& element,
                                                          const IntersectionIterator& intersectionIt,
                                                          const LocalPosition& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (globalPos[1] > 1 - eps_)
//        if (globalPos[0] > left_ - eps_ && globalPos[0] < right_ + eps_ && globalPos[1] > 1 - eps_)
            {
                // both components of the velocity must have same type
                values[velocityXIdx] = BoundaryConditions::dirichlet;
                values[velocityYIdx] = BoundaryConditions::dirichlet;
                values[massFracIdx] = BoundaryConditions::dirichlet;
            }
        if ((globalPos[1] < 0.6 || globalPos[1] > 0.9) && (globalPos[0] < 0 + eps_ || globalPos[0] > 5.5 - eps_))
			{
				values[velocityXIdx] = BoundaryConditions::dirichlet;
				values[velocityYIdx] = BoundaryConditions::dirichlet;
                values[massFracIdx] = BoundaryConditions::dirichlet;
			}

        return values;
    }

    SolutionVector dirichlet(const GlobalPosition& globalPos, const Element& element,
                             const IntersectionIterator& intersectionIt,
                             const LocalPosition& localPos) const
    {
        SolutionVector result(0);

        result[velocityXIdx] = velocity(globalPos, element, localPos)[0];
        result[velocityYIdx] = velocity(globalPos, element, localPos)[1];
        result[massFracIdx] = 1e-5;

        return result;
    }

    SolutionVector neumann(const GlobalPosition& globalPos, const Element& element,
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
        // tangential face of porous media
        if (globalPos[1] < 0.5 + eps_)
            alpha = 0.1;//0.1;
        else
        	return(-1.0); // realizes outflow boundary condition at top


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

        if (globalPos[0] > left_ - eps_ && globalPos[0] < right_ + eps_ && globalPos[1] > 1 - eps_)
        	result[velocityYIdx] = -0.05;

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

    Scalar density(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                   const FieldVector<Scalar,dim>& localPos) const
    {
        //TODO: get density from fluid properties!
        Scalar result = 1;
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
    FieldVector<Scalar,dim> gravity(const FieldVector<Scalar,dim>& globalPos) const
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

    TopTwoCStokesProblem(Gas_GL& gasPhase, Matrix2p<Grid, Scalar>& soil, MultiComp& multicomp = *(new CWaterAir))
        :
        StokesTransportProblem<Grid,Scalar>(gasPhase, multicomp),
        gasPhase_(gasPhase),
        soil_(soil),
        multicomp_(multicomp)
    {
        eps_ = 1e-8;

        left_ = 2.0;
        right_ = 3.5;

        for (int i=0; i<dim; ++i)
            gravity_[i] = 0;
        //        gravity_[dim] = -9.81;
    }


protected:
    FieldVector<Scalar,dim> gravity_;
    Gas_GL& gasPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    MultiComp& multicomp_;
    Scalar eps_;
    Scalar left_, right_;
};

}
#endif
