#ifndef DUNE_INJECTIONPROBLEM_HH
#define DUNE_INJECTIONPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/2p2c/2p2cproblem.hh>

/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
  //! class that defines the parameters of an air injection under a low permeable layer
  /*! Problem definition of an air injection under a low permeable layer. Air enters the domain
   * at the right boundary and migrates upwards.
   * Problem was set up using the rect2d.dgf grid.
   *
   *	Template parameters are:
   *
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values
   */
  template<class Grid, class RT>
  class InjectionProblem : public TwoPTwoCProblem<Grid, RT>
  {
	  enum {dim=Grid::dimension, m=2};
	  typedef typename Grid::ctype Scalar;
	  typedef typename Grid::Traits::template Codim<0>::Entity Element;
	  typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	  enum {pWIdx = 0, switchIdx = 1}; // phase index
	  enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // phase state


/////////////////////////////
// TYPE of the boundaries
/////////////////////////////
	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<Scalar,dim>& x, const Element& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

		if (x[0] < eps_)
			values = BoundaryConditions::dirichlet;
//		if (x[1] < eps_)
//			values = BoundaryConditions::dirichlet;

		return values;
	}

/////////////////////////////
// DIRICHLET boundaries
/////////////////////////////
	virtual FieldVector<RT,m> g (const FieldVector<Scalar,dim>& x, const Element& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<RT,m> values(0);
		RT densityW_ = 1000.0;

		values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - x[1]);
		values[switchIdx] = 1e-8;  // may be Sn, Xaw or Xwn!!

//		if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//		 && x[0] >= innerLowerLeft_[0])
//			values[switchIdx] = 0.2;
//		else
//			values[switchIdx] = 1e-6;

		return values;
	}

/////////////////////////////
// NEUMANN boundaries
/////////////////////////////
	virtual FieldVector<RT,m> J (const FieldVector<Scalar,dim>& x, const Element& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		//RT lambda = (x[1])/height_;

		if (x[1] < 15 && x[1] > 5)
			values[switchIdx] = -1e-3;

		return values;
	}

/////////////////////////////
// INITIAL values
/////////////////////////////
	virtual FieldVector<RT,m> initial (const FieldVector<Scalar,dim>& x, const Element& e,
				const FieldVector<Scalar,dim>& xi) const
	{

		FieldVector<RT,m> values;
		RT densityW_ = 1000.0;

		values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - x[1]);
		values[switchIdx] = 1e-8;

//		if ((x[0] > 60.0 - eps_) && (x[1] < 10 && x[1] > 5))
//			values[switchIdx] = 0.05;

//			if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//			 && x[0] >= innerLowerLeft_[0])
//				values[switchIdx] = 0.2;
//			else
//				values[switchIdx] = 1e-6;

		return values;
	}


	int initialPhaseState (const FieldVector<Scalar,dim>& x, const Element& e,
				  const FieldVector<Scalar,dim>& xi) const
	{

		enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
		int state;

//		state = bothPhases;
		state = waterPhase;

//		if ((x[0] > 60.0 - eps_) && (x[1] < 10 && x[1] > 5))
//			state = bothPhases;
//			if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//			      && x[0] >= innerLowerLeft_[0])
//				state = 2;
//			else

		return state;
	}

/////////////////////////////
// sources and sinks
/////////////////////////////
	virtual FieldVector<RT,m> q (const FieldVector<Scalar,dim>& x, const Element& e,
					const FieldVector<Scalar,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

//////////////////////////////

	virtual FieldVector<RT,dim> gravity () const
	{
		return gravity_;
	}

	double depthBOR () const
	{
		return depthBOR_;
	}

	InjectionProblem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, RT>& soil,
			const FieldVector<Scalar,dim> outerLowerLeft = 0., const FieldVector<Scalar,dim> outerUpperRight = 0.,
			const FieldVector<Scalar,dim> innerLowerLeft = 0., const FieldVector<Scalar,dim> innerUpperRight = 0.,
			const RT depthBOR = 0., TwoPhaseRelations<Grid, RT>& law = *(new TwoPhaseRelations<Grid, RT>),
			MultiComp& multicomp = *(new CWaterAir))
	: TwoPTwoCProblem<Grid,RT>(liq, gas, soil, multicomp, law),
	  outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
	  depthBOR_(depthBOR), eps_(1e-8*outerUpperRight[0])
	  {
		height_ = outerUpperRight[1] - outerLowerLeft[1];
		width_ = outerUpperRight[0] - outerLowerLeft[0];

		gravity_[0] = 0;
		gravity_[1] = -9.81;
	  }

  private:
	  FieldVector<Scalar,dim> outerLowerLeft_, outerUpperRight_;
	  FieldVector<Scalar,dim> innerLowerLeft_, innerUpperRight_;
	  Scalar width_, height_;
	  Scalar depthBOR_, eps_;
//	  RT densityW_, densityN_;
	  FieldVector<Scalar,dim> gravity_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

  template<class Grid, class RT>
  class InjectionSoil: public Matrix2p<Grid,RT>
  {
  public:
  	typedef typename Grid::Traits::template Codim<0>::Entity Element;
  	typedef typename Grid::ctype Scalar;
  	enum {dim=Grid::dimension, m=1};

  	virtual FieldMatrix<Scalar,dim,dim> K (const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi)
  	{
  		if (x[1] < layerBottom_)
  			return highK_;
  		else
  			return lowK_;
  	}
  	virtual double porosity(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		return 0.3;
  	}

  	virtual double Sr_w(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi, const double T) const
  	{
  		return 0.2;
  	}

  	virtual double Sr_n(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi, const double T) const
  	{
  		return 0.0;
  	}

  	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
  			 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
  	virtual double heatCap(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		return 	790 /* spec. heat cap. of granite */
  						* 2700 /* density of granite */
  						* porosity(x, e, xi);
  	}

  	virtual double heatCond(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi, const double sat) const
  	{
  		static const double lWater = 0.6;
  		static const double lGranite = 2.8;
  		double poro = porosity(x, e, xi);
  		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
  		double ldry = pow(lGranite, (1-poro));
  		return ldry + sqrt(sat) * (ldry - lsat);
  	}

  	virtual std::vector<double> paramRelPerm(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi, const double T) const
  	{
  		// example for Brooks-Corey parameters
  		std::vector<double> param(2);
  		param[0] = 2.; // lambda
  		param[1] = 1e4; // entry-pressures

  		return param;
  	}

  	virtual typename Matrix2p<Grid,RT>::modelFlag relPermFlag(const FieldVector<Scalar,dim>& x, const Element& e, const FieldVector<Scalar,dim>& xi) const
  	{
  		return Matrix2p<Grid,RT>::brooks_corey;
  	}

  	InjectionSoil():Matrix2p<Grid,RT>()
  	{
  		lowK_ = highK_ = 0.;
  		for(int i = 0; i < dim; i++){
  			lowK_[i][i] = 1e-13;
  			highK_[i][i] = 1e-12;
  		}
  		layerBottom_ = 22.0;
  	}

  	~InjectionSoil()
  	{}

  private:
  	FieldMatrix<Scalar,dim,dim> lowK_, highK_;
	double layerBottom_;
  };


} //end namespace
#endif
