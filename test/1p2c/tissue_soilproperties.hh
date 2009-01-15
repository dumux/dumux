/// $Id$
#ifndef TISSUE_SOILPROPERTIES
#define TISSUE_SOILPROPERTIES

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

template<class G, class RT>
class TissueSoil: public Matrix2p<G, RT>
{
public:
typedef	typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::ctype DT;
	enum
	{	dim=G::dimension, m=2};

	const FieldMatrix<DT,dim,dim> &K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) 
	{
	     if (10<x[0] && x[0]<12 && 10<x[1] && x[1]<12)
		   return permloc_;						//tumor tissue
	     else
		   return permlocWell_;					//healthy tissue
	}

	double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		if (10<x[0] && x[0]<12 && 10<x[1] && x[1]<12)
		   return 0.31;						//tumor tissue
		else
		   return 0.13;					//healthy tissue
	}

	double tortuosity (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
		{
			if (10<x[0] && x[0]<12 && 10<x[1] && x[1]<12)
			   return 0.706;						//tumor tissue
			else
			   return 0.280;					//healthy tissue
		}
	
	virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		return 0;
	}
	
	virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const 
	{
		return 0;
	}
	
	virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
	{
		std::vector<double> param(0);
		return param;
	}
	
	TissueSoil()
	:Matrix2p<G,RT>()
	{
	  permloc_ = 0;
	  permlocWell_ = 0;
	  for (int i = 0; i < dim; i++)
		permloc_[i][i] = 2.142e-11; 			//tumor tissue
      for (int i = 0; i < dim; i++)
        permlocWell_[i][i] = 4.424e-12;		//healthy tissue
	}

	private:
	Dune::FieldMatrix<DT,dim,dim> permloc_, permlocWell_;
};
} // end namespace
#endif
