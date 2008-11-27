/// $Id: co2_soilproperties.hh 610 2008-09-18 17:04:34Z melanie $
#ifndef LEAKYWELL_SOILPROPERTIES
#define LEAKYWELL_SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class G, class RT>
class LeakyWellSoil: public HomogeneousSoil<G, RT>
{
public:
typedef	typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::ctype DT;
	enum
	{	dim=G::dimension, m=3};

	FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{
	     if (x[0] > -0.2 && x[0] < 0.2 && x[1] > -0.2 && x[1] < 0.2)
		   return permlocWell_;
//	     else if (x[2] > 30.+1.e-2 && x[2] < 130-1.e-2)
//	       return permlocAquitard_;
	     else
		   return permloc_;
	}

	double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		return 0.15;
	}

	double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.2;
	}

	double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.05;
	}

	virtual double heatCap(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		return 	(800 /* spec. heat cap. of sediment */
						* 2650 /* density of sediment */
						* (1-porosity(x, e, xi)));
	}

	virtual double heatCond(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double sat) const
	{
		static const double ldry = 0.32;
		static const double lsat = 2.7;
		double l_wurz, Sw;
		Sw = sat;
		if(Sw<0.0) Sw = 0.0; /* ACHTUNG Regularisierung */
		if(Sw>1.) Sw = 1.; /* ACHTUNG Regularisierung */

		l_wurz = ldry + sqrt(Sw)*(lsat-ldry);

		if(isnan(l_wurz)) {
			std::cout <<"isnan heatcondwurzel \n"<<std::endl;
			l_wurz = 0.0;
		}
		return(l_wurz);
	}

	std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		std::vector<double> param(2);

		//Brooks-Corey parameters
		param[0] = 2.; // lambda
		param[1] = 10000.; // entry-pressure

		return param;
	}

	typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		return Matrix2p<G,RT>::brooks_corey;
	}


	LeakyWellSoil()
	:HomogeneousSoil<G,RT>()
	{
	  permloc_ = 0;
	  permlocWell_ = 0;
	  for (int i = 0; i < dim; i++)
		permloc_[i][i] = 2.0e-14;
      for (int i = 0; i < dim; i++)
        permlocWell_[i][i] = 1.0e-12;
      for (int i = 0; i < dim; i++)
    	permlocAquitard_[i][i] = 1.0e-18;
	}

	private:
	Dune::FieldMatrix<DT,dim,dim> permloc_, permlocWell_, permlocAquitard_;
};
} // end namespace
#endif
