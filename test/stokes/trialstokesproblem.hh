// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef TRIALSTOKESPROBLEM_HH
#define TRIALSTOKESPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class TrialStokesProblem : public StokesProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    enum {velocityXIdx=0, velocityYIdx=1, velocityZIdx=2, pressureIdx=dim};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

public:
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (globalPos[0] > left_ - eps_ && globalPos[0] < right_ + eps_ && globalPos[1] > 1 - eps_)
//        if (globalPos[1] > 1 - eps_)
			{
                // both components of the velocity must have same type
                values[velocityXIdx] = BoundaryConditions::dirichlet;
                values[velocityYIdx] = BoundaryConditions::dirichlet;
            }
        if (globalPos[0] < 0 + eps_ || globalPos[0] > 5.5 - eps_)
			{
				values[velocityXIdx] = BoundaryConditions::dirichlet;
				values[velocityYIdx] = BoundaryConditions::dirichlet;
			}
//        if (globalPos[1] < 0.5 + eps_) // lower boundary
//			{
//				values[velocityXIdx] = BoundaryConditions::dirichlet;
//				values[velocityYIdx] = BoundaryConditions::dirichlet;
//			}

        return values;
    }

    virtual FieldVector<Scalar,numEq>dirichlet(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
    	FieldVector<Scalar,numEq> result(0);

        result[velocityXIdx] = velocity(globalPos)[0];
        result[velocityYIdx] = velocity(globalPos)[1];

        return result;
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& globalPos, const Element& element,
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
        Scalar alpha;
        // tangential face of porous media
        if (globalPos[1] < 0.5 + eps_)
            alpha = 10.0;
        else
        	return(-1.0); // realizes outflow boundary condition at top


        //TODO: uses only Kxx, extend to permeability tensor
        Scalar permeability = 1e-2;

        return sqrt(permeability)/(alpha*mu(globalPos, element, localPos));
    }

    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
    {
        return 0.01;
    }

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
    	FieldVector<Scalar,numEq> result(0);

        if (globalPos[0] > left_ - eps_ && globalPos[0] < right_ + eps_ && globalPos[1] > 1 - eps_)
        	{
				result[velocityXIdx] = 0;
				result[velocityYIdx] = (globalPos[0] - 2.5)*(globalPos[0] - 3.0);//4.0*globalPos[1]*(1.0 - globalPos[1]);
        	}

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (0);
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldMatrix<Scalar, dim, dim> result(0);
        result[0][1] = -1.0;
        result[1][0] = -1.0;

        return result;
    }

    TrialStokesProblem()
    {
        eps_ = 1e-8;

        left_ = 2.5;
        right_ = 3.0;

        for (int i=0; i<dim; ++i)
            gravity_[i] = 0;
        //        gravity_[dim] = -9.81;
    }
protected:
    FieldVector<Scalar,dim> gravity_;
//    Gas_GL& gasPhase_;
//    Matrix2p<Grid, Scalar>& soil_;
//    MultiComp& multicomp_;
    Scalar eps_;
    Scalar left_, right_;
};

}
#endif
