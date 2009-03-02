/// $Id: co2_soilproperties.hh 610 2008-09-18 17:04:34Z melanie $
#ifndef MINCCO2_SOILPROPERTIES
#define MINCCO2_SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class Grid, class Scalar>
class MincCO2Soil: public HomogeneousSoil<Grid, Scalar>
{
public:
typedef	typename Grid::Traits::template Codim<0>::Entity Entity;
	typedef typename Grid::ctype DT;
	enum
	{	dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef Dune::FieldVector<DT,dim> LocalPosition;
    typedef Dune::FieldVector<DT,dimWorld> GlobalPosition;

//	FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
//	{
////	     if (x[0] > -0.2 && x[0] < 0.2 && x[1] > -0.2 && x[1] < 0.2)
////		   return permlocWell_;
//
//		   return permloc_;
//	}
	FieldMatrix<DT,dim,dim> KFracture (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{
//	     if (x[0] > -0.2 && x[0] < 0.2 && x[1] > -0.2 && x[1] < 0.2)
//		   return permlocWell_;

		   return K_Fracture;
	}
	FieldMatrix<DT,dim,dim> KMatrix (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{
//	     if (x[0] > -0.2 && x[0] < 0.2 && x[1] > -0.2 && x[1] < 0.2)
//		   return permlocWell_;

		   return K_Matrix;
	}


	double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		return 0.15;
	}
	double porosityFracture(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
		{
			return 0.25;
		}
	double porosityMatrix(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
		{
			return 0.15;
		}

	double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.2;
	}
	double Sr_wFracture(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.2;
		}
	double Sr_wMatrix(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.2;
		}

	double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.05;
	}
	double Sr_nFracture(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.05;
		}
	double Sr_nMatrix(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
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

	typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		return Matrix2p<Grid,Scalar>::brooks_corey;
	}


	MincCO2Soil()
	:HomogeneousSoil<Grid,Scalar>()
	{
	  permloc_ = 0;
	  permlocWell_ = 0;
	  	K_Fracture = 0.;
	  	K_Matrix  = 0.;
	  for (int i = 0; i < dim; i++)
		{permloc_[i][i] = 2.0e-12;
	    K_Fracture[i][i]= 2.0e-12;
	    K_Matrix[i][i]= 2.0e-12;}
          for (int i = 0; i < dim; i++){
               permlocWell_[i][i] = 1.0e-10;
          K_FractureWell[i][i] = 1.0e-10;
          K_MatrixWell[i][i] = 1.0e-10;}
	}

	private:
	Dune::FieldMatrix<DT,dim,dim> permloc_, permlocWell_;
	Dune::FieldMatrix<DT,dim,dim> K_Fracture, K_FractureWell, K_Matrix, K_MatrixWell;

};
} // end namespace
#endif
