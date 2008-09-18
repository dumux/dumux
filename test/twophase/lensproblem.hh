// $Id: lensproblem.hh 566 2008-09-11 11:38:31Z bernd $

#ifndef DUNE_LENSPROBLEM_HH
#define DUNE_LENSPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/property_baseclasses.hh>
#include<dumux/twophase/twophaseproblem.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
  //! base class that defines the parameters of a diffusion equation
  /*! An interface for defining parameters for the stationary diffusion equation
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
   * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
   * on \f$\Gamma_2\f$. Here,
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
   * and \f$\lambda\f$ the total mobility, possibly depending on the
   * saturation.
   *
   *	Template parameters are:
   *
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values
   */
  template<class G, class RT>
  class LensProblem : public TwoPhaseProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {dim=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	enum {pWIdx = 0, sNIdx = 1};

	virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
					const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<DT,dim>& xi) const
	{
		FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

		if (x[0] < outerLowerLeft_[0] + eps_ || x[0] > outerUpperRight_[0] - eps_) {
			//std::cout << "Dirichlet: " << x << std::endl;
			values = BoundaryConditions::dirichlet;
		}
//		if (values[0] == BoundaryConditions::dirichlet)
//			std::cout << "Dirichlet: " << x[0] << ", " << x[1] << std::endl;
//		else
//			std::cout << "Neumann: " << x[0] << ", " << x[1] << std::endl;

		return values;
	}

	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		if (x[0] < outerLowerLeft_[0] + eps_) {
			RT a = -(1 + 0.5/height_);
			RT b = -a*outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = 0.0;//outerSnr_;
		}
		else {
			RT a = -1;
			RT b = outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = 0.0;//outerSnr_;
		}

		return values;
	}

	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		RT lambda = (outerUpperRight_[0] - x[0])/width_;
		if (lambda > 0.5 && lambda < 2.0/3.0 && x[1] > outerUpperRight_[1] - eps_) {
			values[sNIdx] = -0.04;
		}

		return values;
	}

	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{

		FieldVector<RT,m> values;

		values[pWIdx] = -densityW_*gravity_[1]*(height_ - x[1]);

		if (x[0] < outerLowerLeft_[0] + eps_) {
			RT a = -(1 + 0.5/height_);
			RT b = -a*outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
		}

		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
			values[sNIdx] = 0.0;//innerSnr_;
		else
			values[sNIdx] = 0.0;//outerSnr_;

		return values;
	}

	virtual FieldVector<RT,dim> gravity () const
	{
		return gravity_;
	}

	LensProblem(Fluid& liq1, Fluid& liq2, Matrix2p<G, RT>& soil,
			const FieldVector<DT,dim> outerLowerLeft = 0., const FieldVector<DT,dim> outerUpperRight = 0,
			const FieldVector<DT,dim> innerLowerLeft = 0., const FieldVector<DT,dim> innerUpperRight = 0,
			TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>))
	: TwoPhaseProblem<G,RT>(liq1, liq2, soil, law),
		outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
		innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
		eps_(1e-8*outerUpperRight[0]),
		densityW_(liq1.density()), densityN_(liq2.density())
	{
		height_ = outerUpperRight[1] - outerLowerLeft[1];
		width_ = outerUpperRight[0] - outerLowerLeft[0];

		gravity_[0] = 0;
		gravity_[1] = -9.81;
	}

	private:
		FieldVector<DT,dim> outerLowerLeft_;
		FieldVector<DT,dim> outerUpperRight_;
		FieldVector<DT,dim> innerLowerLeft_;
		FieldVector<DT,dim> innerUpperRight_;
		DT width_, height_;
		DT eps_;
		RT densityW_, densityN_;
		FieldVector<DT,dim> gravity_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

  template<class G, class RT>
  class LensSoil: public Matrix2p<G,RT>
  {
  public:
  	typedef typename G::Traits::template Codim<0>::Entity Entity;
  	typedef typename G::ctype DT;
  	enum {dim=G::dimension, m=1};

  	virtual FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
  	{
  		if (((x[0]<4.) && (x[0]>1.)) && ((x[1]<3.) && (x[1]>2.)))
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
  		if (((x[0]<4.) && (x[0]>1.)) && ((x[1]<3.) && (x[1]>2.)))
  			return 0.18;
  		else
  			return 0.05;
  	}

  	virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		return 0.0;
  	}

  	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
  			 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
  	virtual double heatCap(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return 	790 /* spec. heat cap. of granite */
  						* 2700 /* density of granite */
  						* porosity(x, e, xi);
  	}

  	virtual double heatCond(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double sat) const
  	{
  		static const double lWater = 0.6;
  		static const double lGranite = 2.8;
  		double poro = porosity(x, e, xi);
  		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
  		double ldry = pow(lGranite, (1-poro));
  		return ldry + sqrt(sat) * (ldry - lsat);
  	}

  	virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		// example for Brooks-Corey parameters
		std::vector<double> param(5);

		if (((x[0]<4.) && (x[0]>1.)) && ((x[1]<3.) && (x[1]>2.)))
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
//		if (((x[0]<4.) && (x[0]>1.)) && ((x[1]<3.) && (x[1]>2.)))
//  		{
//  			param[0] = 2.; // lambda
//  			param[1] = 1e4; // entry-pressures
//  		}
//		else
//		{
//			param[0] = 2.; // lambda
//  			param[1] = 1e4; // entry-pressures
//		}

  		return param;
  	}

  	virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return Matrix2p<G,RT>::van_genuchten;
  	}

  	LensSoil():Matrix2p<G,RT>()
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

  };

} // end namespace
#endif
