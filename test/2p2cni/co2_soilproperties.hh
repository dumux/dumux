/// $Id$
#ifndef CO2_SOILPROPERTIES_H
#define CO2_SOILPROPERTIES_H

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class Grid, class Scalar>
class CO2Soil: public HomogeneousSoil<Grid, Scalar>
{
public:
typedef	typename Grid::Traits::template Codim<0>::Entity Element;
	typedef typename Grid::ctype CoordScalar;

    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    
    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

	const FieldMatrix<CoordScalar,dim,dim> &K(const GlobalPosition &x, const Element& e, const LocalPosition &xi)
	{
	     if (x[0] > -0.2 && x[0] < 0.2 && x[1] > -0.2 && x[1] < 0.2)
		   return permlocWell_;

		   return permloc_;
	}

	double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
	{
		return 0.15;
	}

	double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		return 0.2;
	}

	double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		return 0.05;
	}

	virtual double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
	{
		return 	(800 /* spec. heat cap. of sediment */
						* 2650 /* density of sediment */
						* (1-porosity(x, e, xi)));
	}

	virtual double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
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

	std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		std::vector<double> param(2);

		//Brooks-Corey parameters
		param[0] = 2.; // lambda
		param[1] = 10000.; // entry-pressure

		return param;
	}

	typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
	{
		return Matrix2p<Grid,Scalar>::brooks_corey;
	}


	CO2Soil()
	:HomogeneousSoil<Grid,Scalar>()
	{
	  permTest_ = 0;
	  permloc_ = 0;
	  permlocWell_ = 0;
	  for (int i = 0; i < dim; i++)
		permTest_[i][i] = 1.0e-15;
	  for (int i = 0; i < dim; i++)
		permloc_[i][i] = 2.0e-14;
      for (int i = 0; i < dim; i++)
        permlocWell_[i][i] = 1.0e-12;
	}

	private:
	Dune::FieldMatrix<Scalar,dim,dim> permloc_, permlocWell_, permTest_;
};
} // end namespace
#endif
