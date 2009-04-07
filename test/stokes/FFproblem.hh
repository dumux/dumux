#ifndef DUNE_FFPROBLEM_HH
#define DUNE_FFPROBLEM_HH

#include"dumux/stokes/stokestrenproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class FFProblem : public StokesTrEnProblem<Grid, Scalar>
{
 public:
	 enum {dim=Grid::dimension, numEq=Grid::dimension+3};
	 typedef typename Grid::Traits::template Codim<0>::Entity Element;
	 typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;


  virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Element& e,
		  const FieldVector<Scalar,dim>& xi) const
  {
	  FieldVector<Scalar,numEq> result(0);

	  Scalar rX = 0.05;
	  if (x[1] >= 0.1)
		  rX = 0.01;

	  result[0] = 0.01;
	  result[1] = 0;
	  result[2] = rX;
	  result[3] = 283.15;
	  result[4] = 1e+5;

	 return result;
  }

  virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                const FieldVector<Scalar,dim>& xi) const
  {
    FieldVector<Scalar,numEq> result(0);

    if ((x[0] >= 0.2)&&(x[0] <= 0.8)&&(x[1] <= 0.02))
    	result[3] = 1e-3;

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                        const IntersectionIterator& intersectionIt,
                        const FieldVector<Scalar,dim>& xi) const
  {
//	  if (x[0] == 1.0)
//		  return BoundaryConditions::neumann;
//	  else
		  return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<Scalar,dim+2> g(const FieldVector<Scalar,dim>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<Scalar,dim>& xi) const
  {
	  FieldVector<Scalar,dim+2> result(0);
	  Scalar xV = 0;

	  if ((x[0] == 0)||(x[0] == 1.0))
			  xV = 0.01;

	  Scalar rX = 0.05;
	  if (x[1] >= 0.1)
		  rX = 0.01;

	  result[0] = xV;
	  result[2] = rX;
	  result[3] = 283.15;

	  return result;
  }

  virtual FieldVector<Scalar,dim+2> J(const FieldVector<Scalar,dim>& x, const Element& e,
                const IntersectionIterator& intersectionIt,
                const FieldVector<Scalar,dim>& xi)
  {
      FieldVector<Scalar,dim+2> result(0);
      return result;
  }

  virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& x) const
  {
    FieldVector<Scalar,dim> result(0);

    result[0] = 0.01;

    return result;
  }

  virtual FieldVector<Scalar,numEq> velocitypressure(const FieldVector<Scalar,dim>& x) const
  {
	FieldVector<Scalar,numEq> result(0);
	FieldVector<Scalar,dim> vel(velocity(x));
	result[0] = vel[0];
	result[1] = vel[1];
	result[2] = 1e+5;

	return result;
  }

  virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
  {
    Scalar result = 1e+5;

    return result;
  }

  virtual Scalar temperature(const FieldVector<Scalar,dim>& x) const
  {
      Scalar result = 283.15;

      return result;
  }

  virtual Scalar heatConductivity(const FieldVector<Scalar,dim>& x, const Element& e,
  				    const FieldVector<Scalar,dim>& xi) const
  {
  	  Scalar result = 0.025;
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
	  Scalar result = 0;
	  return result;
  }

   virtual Scalar partialdensity(const FieldVector<Scalar,dim>& x) const
   {
	   Scalar result = 0.01;

	   if (x[1] >= 0.1)
		   result = 0.05;

       return result;
     }

   virtual Scalar density(const FieldVector<Scalar,dim>& x, const Element& e,
           const FieldVector<Scalar,dim>& xi) const
   {
	   Scalar result = 0.15;

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

	   Scalar rX = 0.01;
	   if (x[1] >= 0.1)
		   rX = 0.05;

	   result[0] = 0.01;
	   result[2] = rX;
	   result[3] = 283.15;
	   result[4] = 1e+5;

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

   virtual FieldVector<Scalar,dim+2> velocitymassfractemp(const FieldVector<Scalar,dim>& x) const
      {
   	   FieldVector<Scalar,dim+2> result(0);
   	   FieldVector<Scalar,dim> velocityValue(velocity(x));

   	   result[0] = velocityValue[0];
   	   result[1] = velocityValue[1];
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
				res[kx][ky] = 1.6e-4;

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

  FFProblem()
  {
      pi = 4.0*atan(1.0);
  }

  FFProblem(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
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
