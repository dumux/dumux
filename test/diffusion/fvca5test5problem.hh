// $Id$

#ifndef FVariableClassA5TEST5PROBLEM_HH
#define FVariableClassA5TEST5PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/material/fluids/uniform.hh"
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/property_baseclasses.hh>
#include <dune/istl/bvector.hh>

namespace Dune
{
FieldMatrix<double,2,2> perm;

template<class Grid, class Scalar>
class FVariableClassA5Test5Soil: public Matrix2p<Grid, Scalar>
{
public:
    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Matrix2p<Grid, Scalar>::modelFlag ModelFlag;
    typedef FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];

        perm[0][0] = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        perm[0][1] = perm[1][0] = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        perm[1][1] = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;

        return perm;
    }

    virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.2;
    }

    virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.0;
    }

    virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.0;
    }

    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {

        std::vector<double> param(2);

        param[0] = 0;
        param[1] = 0;

        return param;
    }

    virtual ModelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid, Scalar>::linear;
    }

    FVariableClassA5Test5Soil(const double delta)
    : Matrix2p<Grid,Scalar>(), delta_(delta)
    {}

private:
	const double delta_;
};


//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class GridView, class Scalar, class VariableClass>
class FVariableClassA5Test5Problem: public DiffusionProblem<GridView, Scalar, VariableClass>
{
protected:
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    FVariableClassA5Test5Problem(VariableClass& variables, Matrix2p<Grid, Scalar>& soil, const double delta)
    : DiffusionProblem<GridView,Scalar,VariableClass>(variables,
    		*(new TwoPhaseRelations<Grid,Scalar>(soil, *(new Uniform), *(new Uniform)))),
    delta_(delta)
    {}

    Scalar sourcePress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        double pi = 4.0*atan(1.0);
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        double ux = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        double uy = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);
        double kxx = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        double kxy = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        double kyy = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
        double f0 = sin(pi*globalPos[0])*sin(pi*globalPos[1])*pi*pi*(1.0 + delta_)*(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])
            + cos(pi*globalPos[0])*sin(pi*globalPos[1])*pi*(1.0 - 3.0*delta_)*globalPos[0]
            + cos(pi*globalPos[1])*sin(pi*globalPos[0])*pi*(1.0 - 3.0*delta_)*globalPos[1]
            + cos(pi*globalPos[1])*cos(pi*globalPos[0])*2.0*pi*pi*(1.0 - delta_)*globalPos[0]*globalPos[1];

        return ((f0 + 2.0*(globalPos[0]*(kxx*ux + kxy*uy) + globalPos[1]*(kxy*ux + kyy*uy)))/rt);
    }

    typename BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return BoundaryConditions::dirichlet;
    }

    Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return (exact(globalPos));
    }

    Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
    {
        FieldVector<Scalar,dim> grad(0);
        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 1.0;
    }

private:
	const double delta_;
};
}

#endif
