#ifndef DUNE_LENSSOIL_HH
#define DUNE_LENSSOIL_HH

/**
 * @file
 * @brief  Class for defining an instance of a Matrix2p soil
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{

  template<class G, class RT>
  class LensSoil: public Matrix2p<G,RT>
  {
  public:
  	typedef typename G::Traits::template Codim<0>::Entity Entity;
  	typedef typename G::ctype DT;
  	enum {dim=G::dimension, m=1};

  	// define PERMEABILITY tensor
  	virtual FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
  	{
		if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
		 && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
			return Kin_;
  		else
  			return Kout_;
  	}
  	virtual double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return 0.4;
  	}

  	virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
		if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
		 && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
  			return 0.18;
  		else
  			return 0.05;
  	}

  	virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		return 0.0;
  	}

  	virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return Matrix2p<G,RT>::van_genuchten;
  	}

  	virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		// example for Brooks-Corey parameters
		std::vector<double> param(5);

		if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
		 && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
  		{
			param[0] = 1-1/4.7;
			param[1] = 4.7;
			param[2] = 0.5;
			param[3] = 1/3.;
			param[4] = 0.0037;
  		}
		else
		{
			param[0] = 1-1/7.3;
			param[1] = 7.3;
			param[2] = 1/2.;
			param[3] = 1/3.;
			param[4] = 0.00045;

		}

  		return param;
  	}

  	LensSoil(const FieldVector<DT,dim>& outerLowerLeft = 0., const FieldVector<DT,dim>& outerUpperRight = 0,
  			 const FieldVector<DT,dim>& innerLowerLeft = 0., const FieldVector<DT,dim>& innerUpperRight = 0)
  	: Matrix2p<G,RT>(),	outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
		innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight)
  	{
  		Kin_ = Kout_ = 0;
  		for(int i = 0; i < dim; i++)
  		{
  			Kin_[i][i] = 1e-13;
			Kout_[i][i] = 5e-10;
  		}
  	}

  	~LensSoil()
  	{}

  private:
  	FieldMatrix<DT,dim,dim> Kin_;
  	FieldMatrix<DT,dim,dim> Kout_;
	FieldVector<DT,dim> outerLowerLeft_;
	FieldVector<DT,dim> outerUpperRight_;
	FieldVector<DT,dim> innerLowerLeft_;
	FieldVector<DT,dim> innerUpperRight_;
  };

} // end namespace
#endif

