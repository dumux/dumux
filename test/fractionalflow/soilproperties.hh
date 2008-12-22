#ifndef MATRIXPROPERTIES
#define MATRIXPROPERTIES

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

template<class Grid, class Scalar>
class HeterogeneousSoil: public Matrix2p<Grid, Scalar>
{
public:
	enum
	{dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
	typedef	typename Grid::Traits::template Codim<0>::Entity Element;
	typedef typename Matrix2p<Grid, Scalar>::modelFlag modelFlag;
	typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
	typedef Dune::FieldVector<Scalar, dim> LocalPosition;
	typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

	virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos)
	{
		return K_;
	}
	virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
	{
		return 0.2;
	}

	virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
	{
		return 0.2;
	}

	virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
	{
		return 0.2;
	}

	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
	 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
	virtual double heatCap(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
	{
		return 790 /* spec. heat cap. of granite */
		* 2700 /* density of granite */
		* porosity(globalPos, element, localPos);
	}

	virtual double heatCond(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double sat) const
	{
		static const double lWater = 0.6;
		static const double lGranite = 2.8;
		double poro = porosity(globalPos, element, localPos);
		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
		double ldry = pow(lGranite, (1-poro));
		return ldry + sqrt(sat) * (ldry - lsat);
	}

	virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
	{

		std::vector<double> param(2);
		if (globalPos[0]<=300)
		{
			//linear parameters
			param[0] = 0.2;
			param[1] = 0.;
		}
		else
		{
			//Brooks-Corey parameters
			param[0] = 3; // lambda
			param[1] = 0.; // entry-pressure
		}
		return param;
	}

	virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
	{
//		if (x[0]<=300)
//			return 1;
//		else
			return Matrix2p<Grid, Scalar>::brooks_corey;
	}

	HeterogeneousSoil()
        :Matrix2p<Grid,Scalar>(),K_(1e-7)//,permeability(g, name, create)
	{}

	~HeterogeneousSoil()
	{}

private:
    FieldMatrix K_;

public:

//	RandomPermeability<Grid> permeability;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
