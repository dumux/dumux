// $Id: lensproblem.hh 566 2008-09-11 11:38:31Z bernd $

#ifndef DUNE_RICHARDSPROBLEM_HH
#define DUNE_RICHARDSPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/property_baseclasses.hh>
#include"dumux/richards/richardsproblem.hh"

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
  class RichardsSwProblem : public RichardsProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {dim=G::dimension, m=1};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	enum {sWIdx =0 };

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

		switch (intersectionIt.boundaryId()) {
		case 1: case 2: case 3: case 4: case 5:
			values = Dune::BoundaryConditions::neumann;
			break;
		case 6:
			values = Dune::BoundaryConditions::dirichlet;
			break;
		}

		return values;
	}

	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		switch (intersectionIt.boundaryId()) {
		case 6:
//			values[pWIdx] = 1.0e+5 - densityW_*gravity_[2]*(height_-x[2]);
			values[sWIdx] = 0.05;
			break;
			}

		return values;
	}

	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		switch (intersectionIt.boundaryId()) {
		case 1: case 2: case 3: case 4:
			values[sWIdx] = 0;
			break;
		case 5:
			values[sWIdx] = -0.4;
			break;
		}


		return values;
	}

	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{

		FieldVector<RT,m> values;

//		values[pWIdx] = 1.0e+5 - densityW_*gravity_[2]*(height_-x[2]);
		values[sWIdx] = 0.05;

		return values;
	}

	virtual FieldVector<RT,dim> gravity () const
	{
		return gravity_;
	}

	RichardsSwProblem(Fluid& liq1, Matrix2p<G, RT>& soil, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>))
	: RichardsProblem<G,RT>(liq1, soil, law),
		densityW_(liq1.density())
	{
		height_ = 20;
		width_ = 20;

		gravity_[0] = 0;
		gravity_[1] = 0;
		gravity_[dim-1] = -9.81;
	}

	private:
		DT width_, height_;
		RT densityW_, densityN_;
		FieldVector<DT,dim> gravity_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

  template<class G, class RT>
  class RichardsSoil: public Matrix2p<G,RT>
  {
  public:
  	typedef typename G::Traits::template Codim<0>::Entity Entity;
  	typedef typename G::ctype DT;
  	enum {dim=G::dimension, m=1};

  	virtual FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
  	{
  			return Kout_;
  	}
  	virtual double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return 0.4;
  	}

  	virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
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
  		// VanGenuchten
//		std::vector<double> param(5);
//			param[0] = 1-1/7.3; // n
//			param[1] = 7.3;     // m
//			param[2] = 1/2.;    // epsilon
//			param[3] = 1/3.;    // gama
//			param[4] = 0.00045; // alpha

  		// linear
			std::vector<double> param(2);
				param[0] = 0; 	 // pCmin
				param[1] = 1e+4; // pCMax


  		return param;
  	}

  	virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
//  		return Matrix2p<G,RT>::van_genuchten;
  		return Matrix2p<G,RT>::linear;
  	}

  	RichardsSoil():Matrix2p<G,RT>()
  	{
  		Kout_ = 0;
  		for(int i = 0; i < dim; i++)
  		{
			Kout_[i][i] = 5e-10;
  		}
  	}

  	~RichardsSoil()
  	{}

  private:
  	FieldMatrix<DT,dim,dim> Kout_;

  };

} // end namespace
#endif
