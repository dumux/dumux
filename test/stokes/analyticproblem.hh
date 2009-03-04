#ifndef DUNE_ANALYTICPROBLEM_HH
#define DUNE_ANALYTICPROBLEM_HH

#include"dumux/stokes/stokestransportproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class AnalyticProblem : public StokesTransportProblem<Grid, Scalar>
{
 public:
	 enum {dim=Grid::dimension, numEq=Grid::dimension+2};
	 typedef typename Grid::Traits::template Codim<0>::Entity Element;
	 typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;


  virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Element& e,
		  const FieldVector<Scalar,dim>& xi) const
  {
	  return velocitypressuremassfrac(x);
  }

  virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                const FieldVector<Scalar,dim>& xi) const
  {
    FieldVector<Scalar,numEq> result(0);
    result[0] = pi/8.0*cos(2*pi*x[1])*sin(2*pi*x[0]) +
				pi*pi*(1 - 2*cos(2*pi*x[0]))*sin(2*pi*x[1]);
    result[1] = pi/8.0*cos(2*pi*x[0])*sin(2*pi*x[1]) -
				pi*pi*(1 - 2*cos(2*pi*x[1]))*sin(2*pi*x[0]);

    result[2] = 1/4.0*pi*cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0)*
				(2*pi - sin(2*pi*x[0])*sin(pi*x[1]) - sin(pi*x[0])*sin(2*pi*x[1]));

    result[3] = 0.5*pi*pi*cos(2.0*pi*x[0])*cos(2.0*pi*x[1]);

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<Scalar,dim>& xi) const
  {
	  return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<Scalar,dim+1> g(const FieldVector<Scalar,dim>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<Scalar,dim>& xi) const
  {
    return velocitymassfrac(x);
  }

  virtual FieldVector<Scalar,dim+1> J(const FieldVector<Scalar,dim>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<Scalar,dim>& xi)
  {
      FieldVector<Scalar,dim+1> result(0);
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

   virtual Scalar partialdensity(const FieldVector<Scalar,dim>& x) const
   {
	   Scalar result = cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);

       return result;
     }

   virtual FieldVector<Scalar,numEq> velocitypressuremassfrac(const FieldVector<Scalar,dim>& x) const
   {
	   FieldVector<Scalar,numEq> result(0);
	   result[0] =  sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
	   result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);
	   result[2] =  cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);
	   result[3] =  pressure(x);

	   return result;
   }

   virtual FieldVector<Scalar,dim+1> velocitymassfrac(const FieldVector<Scalar,dim>& x) const
   {
	   FieldVector<Scalar,dim+1> result(0);
	   result[0] =  sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);
	   result[1] = -sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]);
	   result[2] =  cos(pi*x[0]/2.0)*sin(pi*x[1]/2.0);

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

  AnalyticProblem()
  {
      pi = 4.0*atan(1.0);
  }

  AnalyticProblem(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
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
