#ifndef DUNE_ANALYTICPROBLEMP_HH
#define DUNE_ANALYTICPROBLEMP_HH

#include"dumux/stokes/stokestransportproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class AnalyticProblemP : public StokesTransportProblem<Grid, Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=Grid::dimension+2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;


    virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Element& e,
                                               const FieldVector<Scalar,dim>& xi) const
    {
        return velocitypressuremassfrac(x);
    }

    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> result(0);

        FieldVector<Scalar,dim> gravityVector(gravity(x));

        result[0] = pi/8.0*cos(2*pi*x[1])*sin(2*pi*x[0]) +
            pi*pi*(1 - 2*cos(2*pi*x[0]))*sin(2*pi*x[1]);

        result[1] = pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1]) -
            pi*pi*(1 - 2*cos(2*pi*x[1]))*sin(2*pi*x[0])
            - gravityVector[dim-1]*density(x,e,xi);

        result[2] = 1/4.0*pi*cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0)*
            (2*pi - sin(2*pi*x[0])*sin(pi*x[1]) - sin(pi*x[0])*sin(2*pi*x[1]));

        result[3] = 0;

        return result;
    }

    /////////////////////////////
    // TYPE of the boundaries
    /////////////////////////////
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<Scalar,dim>& xi) const
    {
    	FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::dirichlet);

/*    	if (x[0] > 1.0 - 1e-6)//&&(x[1] >= 17/50.)&&(x[1] <= 18/50.))
    	{
    		for (int i = 0; i < dim; i++)
    			values[i] = BoundaryConditions::neumann;

    		values[dim+1] = BoundaryConditions::neumann;
    	}
*/
        return values;
    }

    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
        {
            return (-1.0);
        }


/*
    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const
    {
//    	if (x[0] < 1e-6 || x[1] < 1e-6 || x[1] > 1 - 1e-6)
//    		return BoundaryConditions::dirichlet;
//    	return BoundaryConditions::neumann;

        return BoundaryConditions::dirichlet;
    }
*/
    virtual FieldVector<Scalar,numEq>dirichlet(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi) const
    {
    	return velocitymassfrac(x);
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi)
    {
        FieldVector<Scalar,numEq> result(0);

//        result[0] = 1/16.0*(cos(2*pi*x[0])*cos(2*pi*x[1]) + 8*pi*sin(2*pi*x[0])*sin(2*pi*x[1]));
//        result[1] = 1/16.0*cos(2*pi*x[0])*(-8*pi + (1+8*pi)*cos(2*pi*x[1]));

        result[0] = 0.5*pi*sin(2*pi*x[0])*sin(2*pi*x[1]);
        result[1] = -pi*cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);


        return result;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim> result(0);
        result[0] = sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
        result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);

        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitypressure(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        result[0] =  sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
        result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);
        result[2] = -cos(2*pi*x[0])*cos(2*pi*x[1])/16.0;

        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
        Scalar result = -cos(2*pi*x[0])*cos(2*pi*x[1])/16.0;

        return result;
    }

    //virtual FieldVector<Scalar,dim> gravity() const= 0;

    virtual FieldVector<Scalar,dim> gravity(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim> result(0);
        result[dim-1] = 0;//-9.81;

        return result;
    }

    virtual Scalar Qg(const FieldVector<Scalar,dim>& x, const Element& e,
                      const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 0.5*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1]);
        return result;
    }

    virtual Scalar Xg(const FieldVector<Scalar,dim>& x) const
    {
        Scalar result = cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);

        return result;
    }

    virtual Scalar density(const FieldVector<Scalar,dim>& x, const Element& e,
                           const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1;

        return result;
    }

    virtual Scalar viscosity(const FieldVector<Scalar,dim>& x, const Element& e,
                             const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1;

        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitypressuremassfrac(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        FieldVector<Scalar,dim> velocityValue(velocity(x));

        result[0] = 0;//velocityValue[0];
        result[1] = 0;//velocityValue[1];
        result[2] = 0;//partialdensity(x);
        result[3] = 0;//pressure(x);

        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitymassfrac(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        FieldVector<Scalar,dim> velocityValue(velocity(x));

        result[0] = velocityValue[0];
        result[1] = velocityValue[1];
        result[2] = partialdensity(x);

        return result;
    }

    virtual FieldMatrix<Scalar,dim,dim> D (const FieldVector<Scalar,dim>& x, const Element& e,
                                           const FieldVector<Scalar,dim>& xi) const
    {
        FieldMatrix<Scalar,dim,dim> res(0);

        for (int kx=0; kx<dim; kx++)
            for (int ky=0; ky<dim; ky++)
                if (kx == ky)
                    res[kx][ky] = 1.0;

        return res;
    }

    /*
      virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
      {
      FieldMatrix<Scalar, dim, dim> result(0);

      return result;
      }
    */
    virtual Fluid& gasPhase () const
    {
        return gasPhase_;
    }

    virtual MultiComp& multicomp () const
    {
        return multicomp_;
    }

    AnalyticProblemP()
    {
        pi = 4.0*atan(1.0);
    }

    AnalyticProblemP(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
        :
        StokesTransportProblem<Grid,Scalar>(gasPhase,multicomp),
        gasPhase_(gasPhase),
        multicomp_(multicomp)
    {
        pi = 4.0*atan(1.0);
    }


protected:
    double pi;
    Fluid& gasPhase_;
    MultiComp& multicomp_;
};

}
#endif
