#ifndef DUNE_STEPROBLEM_HH
#define DUNE_STEPROBLEM_HH

#include"dumux/stokes/stokestrenproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class SteProblem : public StokesTrEnProblem<Grid, Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=Grid::dimension+3};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;


    virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Element& e,
                                               const FieldVector<Scalar,dim>& xi) const
    {
        return velocitypressuremassfractemp(x);
    }

    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = pi/8.0*cos(2*pi*x[1])*sin(2*pi*x[0]) +
            pi*pi*(1 - 2*cos(2*pi*x[0]))*sin(2*pi*x[1]);

        result[1] = pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1]) -
            pi*pi*(1 - 2*cos(2*pi*x[1]))*sin(2*pi*x[0])
            + gravity()*density(x,e,xi);

        result[2] = -1*1/4.0*pi*cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0)*
            (-2*pi + sin(2*pi*x[0])*sin(pi*x[1]) +
             sin(pi*x[0])*sin(2*pi*x[1]));

        result[3] = pi/128.0*(cos(2*pi*(x[0]-2*x[1])) - cos(pi*(4*x[0]-2*x[1])) +
                              1024*pi*cos(2*pi*(x[0]+x[1])) + cos(2*pi*(2*x[0]+x[1])) -
                              cos(2*pi*(x[0]+2*x[1])));

        result[4] = 1/32.0*pi*(-cos(2*pi*x[0])+cos(2*pi*x[1]))*sin(2*pi*x[0])*sin(2*pi*x[1]);

        //    result[4] = 0.5*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1])// 1/32.0*pi*(-cos(2*pi*x[0])+cos(2*pi*x[1]))*sin(2*pi*x[0])*sin(2*pi*x[1]);
        //                + pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1])*gravity();

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const
    {
        //      if ((x[0] == 1.0)&&(x[1] > 0)&&(x[1] < 1))
        //          return BoundaryConditions::neumann;
        //      else
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,numEq> dirichlet(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        return velocitymassfractemp(x);
    }

    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi)
    {
        FieldVector<Scalar,numEq> result(0);
        return result;
        //      return velocitymassfractempN(x);
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

    virtual Scalar temperature(const FieldVector<Scalar,dim>& x) const
    {
        Scalar result = cos(2*pi*(x[0]+x[1]));

        return result;
    }

    virtual Scalar heatConductivity(const FieldVector<Scalar,dim>& x, const Element& e,
                                    const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1;
        return result;
    }

    //virtual FieldVector<Scalar,dim> gravity() const= 0;

    virtual Scalar gravity() const
    {
        return -9.8;
    }

    virtual Scalar Qg(const FieldVector<Scalar,dim>& x, const Element& e,
                      const FieldVector<Scalar,dim>& xi) const
    {
        //last      Scalar result = 1/32.0*pi*(-cos(2*pi*x[0])+cos(2*pi*x[1]))*sin(2*pi*x[0])*sin(2*pi*x[1]);
        Scalar result = 0.5*pi*pi*cos(2*pi*x[0])*cos(2*pi*x[1])// 1/32.0*pi*(-cos(2*pi*x[0])+cos(2*pi*x[1]))*sin(2*pi*x[0])*sin(2*pi*x[1]);
            + pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1])*gravity();

        //Scalar result = 0.5*pi*cos(pi*(x[0]+x[1])/2.0)*sin(pi*x[0])*sin(pi*(x[0]-x[1]))*sin(pi*x[1]);
        //Scalar result = sin(pi*x[0])*sin(pi*(x[0]-x[1]))*sin(pi*x[1]); //x+y
        // result = -pi*sin(pi*x[0])*sin(pi*(x[0]-x[1]))*sin(pi*x[1])*sin(pi*(x[0]+x[1]));
        return result;
    }

    virtual Scalar partialdensity(const FieldVector<Scalar,dim>& x) const
    {
        Scalar result = cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);

        return result;
    }

    virtual Scalar density(const FieldVector<Scalar,dim>& x, const Element& e,
                           const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1+pressure(x);//sin(pi*(x[0]+x[1])/2.0);

        return result;
    }

    virtual Scalar internalenergy(const FieldVector<Scalar,dim>& x, const Element& e,
                                  const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1;

        return result;
    }

    virtual Scalar enthalpy(const FieldVector<Scalar,dim>& x, const Element& e,
                            const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = internalenergy(x,e,xi) + pressure(x)/density(x,e,xi);

        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitypressuremassfractemp(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        FieldVector<Scalar,dim> velocityValue(velocity(x));

        result[0] = velocityValue[0];
        result[1] = velocityValue[1];
        result[2] = partialdensity(x);
        result[3] = temperature(x);
        result[4] = pressure(x);

        return result;
    }

    virtual FieldVector<Scalar,dim+2> velocitymassfrac(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim+2> result(0);
        FieldVector<Scalar,dim> velocityValue(velocity(x));

        result[0] = velocityValue[0];
        result[1] = velocityValue[1];
        result[2] = partialdensity(x);

        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitymassfractemp(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,numEq> result(0);
        FieldVector<Scalar,dim> velocityValue(velocity(x));

        result[0] = velocityValue[0];
        result[1] = velocityValue[1];
        result[2] = partialdensity(x);
        result[3] = temperature(x);

        return result;
    }

    virtual FieldVector<Scalar,dim+2> velocitymassfractempN(const FieldVector<Scalar,dim>& x) const
    {
        FieldVector<Scalar,dim+2> result(0);
        result[0] = 0.5*pi*sin(2*pi*x[0])*sin(2*pi*x[1]);//-1/16.0*(cos(2*pi*x[0])*cos(2*pi*x[1]) + 8*pi*sin(2*pi*x[0])*sin(2*pi*x[1]));
        result[1] = -pi*cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);//-1/16.0*cos(2*pi*x[0])*(-8*pi + (1+8*pi)*cos(2*pi*x[0]));
        result[2] = partialdensity(x);
        result[3] = temperature(x);
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

    SteProblem()
    {
        pi = 4.0*atan(1.0);
    }

    SteProblem(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
        :
        StokesTrEnProblem<Grid,Scalar>(gasPhase,multicomp),
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
