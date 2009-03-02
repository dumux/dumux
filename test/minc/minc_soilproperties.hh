/// $Id: co2_soilproperties.hh 610 2008-09-18 17:04:34Z melanie $
#ifndef MINCCO2_SOILPROPERTIES
#define MINCCO2_SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class Grid, class Scalar>
class MincLensSoil: public Matrix2p<Grid,Scalar>
{
public:
typedef	typename Grid::Traits::template Codim<0>::Entity Entity;
	typedef typename Grid::ctype DT;
	enum
	{	dim=Grid::dimension};

    virtual const FieldMatrix<DT,dim,dim> &K (const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi)
    {
      if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
       && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
          return Kin_;
        else
            return Kout_;
    }
	
	virtual const FieldMatrix<DT,dim,dim> KFracture (const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi)
	{
		   return K_Fracture;
	}
	virtual const FieldMatrix<DT,dim,dim> KMatrix (const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi)
	{
		   return K_Matrix;
	}


	double porosity(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
	{
		return 0.15;
	}
	virtual double porosityFracture(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
		{
			return 0.25;
		}
	double porosityMatrix(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
		{
			return 0.15;
		}

	double Sr_w(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.2;
	}
	double Sr_wFracture(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.2;
		}
	double Sr_wMatrix(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.2;
		}

	double Sr_n(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0.05;
	}
	double Sr_nFracture(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.05;
		}
	double Sr_nMatrix(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
		{
			return 0.05;
		}

	typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
	{
		return Matrix2p<Grid,Scalar>::brooks_corey;
	}
	
	std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		std::vector<double> param(2);

		//Brooks-Corey parameters
		param[0] = 2.; // lambda
		param[1] = 10000.; // entry-pressure

		return param;
	}

	MincLensSoil(const FieldVector<DT,dim>& outerLowerLeft = 0., const FieldVector<DT,dim>& outerUpperRight = 0,
            const FieldVector<DT,dim>& innerLowerLeft = 0., const FieldVector<DT,dim>& innerUpperRight = 0)
    : Matrix2p<Grid,Scalar>(),    outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
    innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight)
	{
          Kin_ = Kout_ = 0;
          K_Fracture = K_Matrix = 0;
          for(int i = 0; i < dim; i++)
          {
              Kin_[i][i] = 1e-13;
              Kout_[i][i] = 5e-10;
              K_Fracture[i][i] = 1e-8;
              K_Matrix[i][i]=1e-10;
          }
	}

	private:
	FieldMatrix<DT,dim,dim> K_Fracture, K_FractureWell, K_Matrix, K_MatrixWell;
    FieldMatrix<DT,dim,dim> Kin_;
    FieldMatrix<DT,dim,dim> Kout_;
    FieldVector<DT,dim> outerLowerLeft_;
    FieldVector<DT,dim> outerUpperRight_;
    FieldVector<DT,dim> innerLowerLeft_;
    FieldVector<DT,dim> innerUpperRight_;
};
} // end namespace
#endif
