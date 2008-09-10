#ifndef SOILPROPERTIES
#define SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

	template<class G, class RT>
	class TutorialSoil: public HomogeneousSoil<G, RT>
	{
	public:
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		enum
		{	n=G::dimension, m=1};

		FieldMatrix<DT,n,n> K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi)
		{
			return 1e-7;
		}
		double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
		{
			return 0.2;
		}

		double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
		{
			return 0.2;
		}

		double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
		{
			return 0.2;
		}

		std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
		{
			std::vector<double> param(2);

			//Brooks-Corey parameters
			param[0] = 2; // lambda
			param[1] = 0.; // entry-pressure

			return param;
		}

		typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
		{
			return Matrix2p<G,RT>::Brooks_Corey;
		}

		TutorialSoil()
		:HomogeneousSoil<G,RT>()
		{}
	} // end namespace
#endif
