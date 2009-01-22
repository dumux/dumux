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
   *	- Scalar    type used for return values
   */
  template<class Grid, class Scalar>
  class RichardsPwProblem : public RichardsProblem<Grid, Scalar> {
	typedef typename Grid::ctype Scalar;
	enum {dim=Grid::dimension, m=1};
	typedef typename Grid::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	enum {pWIdx =0 };

	virtual FieldVector<Scalar,m> q (const FieldVector<Scalar,dim>& x, const Entity& e,
					const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<Scalar,m> values(0);

		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<Scalar,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<Scalar,dim>& xi) const
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

	virtual FieldVector<Scalar,m> g (const FieldVector<Scalar,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<Scalar,m> values(0);

		switch (intersectionIt.boundaryId()) {
		case 6:
//			values[pWIdx] = 1.0e+5 - densityW_*gravity_[2]*(height_-x[2]);
			values[pWIdx] = -1.0e+6; //- densityW_*gravity_[2]*(height_-x[2]); //-1.0e+6 - densityW_*gravity_[2]*(height_-x[2]);
			break;
			}

		return values;
	}

	virtual FieldVector<Scalar,m> J (const FieldVector<Scalar,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<Scalar,m> values(0);

		switch (intersectionIt.boundaryId()) {
		case 1: case 2: case 3: case 4:
			values[pWIdx] = 0;
			break;
		case 5:
			values[pWIdx] = -0.4;
			break;
		}


		return values;
	}

	virtual FieldVector<Scalar,m> initial (const FieldVector<Scalar,dim>& x, const Entity& e,
				  const FieldVector<Scalar,dim>& xi) const
	{

		FieldVector<Scalar,m> values;

//		values[pWIdx] = 1.0e+5 - densityW_*gravity_[2]*(height_-x[2]);
		values[pWIdx] = -1.0e+6;// - densityW_*gravity_[2]*(height_-x[2]);

		return values;
	}

	virtual FieldVector<Scalar,dim> gravity () const
	{
		return gravity_;
	}

	RichardsPwProblem(Fluid& liq1, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>))
	: TwoPhaseProblem<Grid,Scalar>(liq1, soil, law),
		densityW_(liq1.density())
	{
		height_ = 20;
		width_ = 20;

		gravity_[0] = 0;
		gravity_[1] = 0;
		gravity_[dim-1] = -9.81; //-9.81;
	}

	private:
		Scalar width_, height_;
		Scalar densityW_, densityN_;
		FieldVector<Scalar,dim> gravity_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

  template<class Grid, class Scalar>
  class RichardsSoil: public Matrix2p<Grid,Scalar>
  {
  public:
  	typedef typename Grid::Traits::template Codim<0>::Entity Entity;
  	typedef typename Grid::ctype Scalar;
  	enum {dim=Grid::dimension, m=1};

  	virtual FieldMatrix<Scalar,dim,dim> K (const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi)
  	{
  			return Kout_;
  	}
  	virtual double porosity(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		return 0.4;
  	}

  	virtual double Sr_w(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi, const double T) const
  	{
  			return 0.05;
  	}

  	virtual double Sr_n(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi, const double T) const
  	{
  		return 0.0;
  	}

  	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
  			 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
  	virtual double heatCap(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		return 	790 /* spec. heat cap. of granite */
  						* 2700 /* density of granite */
  						* porosity(x, e, xi);
  	}

  	virtual double heatCond(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi, const double sat) const
  	{
  		static const double lWater = 0.6;
  		static const double lGranite = 2.8;
  		double poro = porosity(x, e, xi);
  		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
  		double ldry = pow(lGranite, (1-poro));
  		return ldry + sqrt(sat) * (ldry - lsat);
  	}

  	virtual std::vector<double> paramRelPerm(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi, const double T) const
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
				param[1] = 1e+6; // pCMax

  		return param;
  	}

  	virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<Scalar,dim>& x, const Entity& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		//  		return Matrix2p<Grid,Scalar>::van_genuchten;
  		  		return Matrix2p<Grid,Scalar>::linear;
  	}

  	RichardsSoil():Matrix2p<Grid,Scalar>()
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
  	FieldMatrix<Scalar,dim,dim> Kout_;

  };

} // end namespace
#endif
