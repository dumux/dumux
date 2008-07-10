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
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/linearlaw.hh>
#include<dumux/twophase/twophaseproblem.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch
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
	enum {n=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	enum {pWIdx = 0, sNIdx = 1};
	enum {swrIdx = 0, snrIdx = 1, alphaIdx = 2, nIdx = 3};

	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi)
	{
		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0] 
		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
			return innerK_;
		else
			return outerK_;
	}

	virtual FieldVector<RT,m> q (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e, 
					const IntersectionIterator& intersectionIt, 
					   const FieldVector<DT,n>& xi) const 
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

	virtual FieldVector<RT,m> g (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);

		if (x[0] < outerLowerLeft_[0] + eps_) {
			RT a = -(1 + 0.5/height_);
			RT b = -a*outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = outerSnr_;
		}
		else {
			RT a = -1;
			RT b = outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = outerSnr_;			
		}
		
		return values;
	}
	  
	virtual FieldVector<RT,m> J (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);

		RT lambda = (outerUpperRight_[0] - x[0])/width_;
		if (lambda > 0.5 && lambda < 2.0/3.0 && x[1] > outerUpperRight_[1] - eps_) {
			values[sNIdx] = -0.04;
		}
		
		return values;
	}
	  
	virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const 
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
			values[sNIdx] = innerSnr_;
		else
			values[sNIdx] = outerSnr_;
	
		return values;
	}

	double porosity (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0] 
		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
			return innerPorosity_;
		else 
			return outerPorosity_;
	}
	
	virtual FieldVector<RT,n> gravity () const 
	{
		return gravity_;
	}
	  
	virtual FieldVector<RT,4> materialLawParameters (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,4> values;
		
		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0] 
		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]) {
			values[swrIdx] = innerSwr_;
			values[snrIdx] = innerSnr_;
			values[alphaIdx] = innerAlpha_;
			values[nIdx] = innerN_;
		}
		else {
			values[swrIdx] = outerSwr_;
			values[snrIdx] = outerSnr_;
			values[alphaIdx] = outerAlpha_;
			values[nIdx] = outerN_;
		}
		
		return values;
	}

	LensProblem(TwoPhaseRelations& law = *(new LinearLaw), 
			const FieldVector<DT,n> outerLowerLeft = 0, const FieldVector<DT,n> outerUpperRight = 0, 
			const FieldVector<DT,n> innerLowerLeft = 0, const FieldVector<DT,n> innerUpperRight = 0, 
//			RT outerK = 4.6e-10, RT innerK = 4.6e-10, 
//			RT outerSwr = 0.05, RT outerSnr = 0, RT innerSwr = 0.05, RT innerSnr = 0, 
//			RT outerPorosity = 0.4, RT innerPorosity = 0.4, 
//			RT outerAlpha = 0.0037, RT innerAlpha = 0.0037, 
//			RT outerN = 4.7, RT innerN = 4.7)
RT outerK = 4.6e-10, RT innerK = 9.05e-13, 
RT outerSwr = 0.05, RT outerSnr = 0, RT innerSwr = 0.18, RT innerSnr = 0, 
RT outerPorosity = 0.4, RT innerPorosity = 0.4, 
RT outerAlpha = 0.0037, RT innerAlpha = 0.00045, 
RT outerN = 4.7, RT innerN = 7.3)
	: TwoPhaseProblem<G, RT>(law), 
	  outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight), 
	  innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight), 
	  eps_(1e-8*outerUpperRight[0]), 
	  densityW_(law.wettingPhase.density()), densityN_(law.nonwettingPhase.density()), 
	  outerSwr_(outerSwr), outerSnr_(outerSnr), innerSwr_(innerSwr), innerSnr_(innerSnr), 
	  outerPorosity_(outerPorosity), innerPorosity_(innerPorosity), 
	  outerAlpha_(outerAlpha), innerAlpha_(innerAlpha), 
	  outerN_(outerN), innerN_(innerN)
	{	
		outerK_[0][0] = outerK_[1][1] = outerK;
		outerK_[0][1] = outerK_[1][0] = 0;
		
		innerK_[0][0] = innerK_[1][1] = innerK;
		innerK_[0][1] = innerK_[1][0] = 0;
		
		height_ = outerUpperRight[1] - outerLowerLeft[1];
		width_ = outerUpperRight[0] - outerLowerLeft[0];
		
		gravity_[0] = 0;
		gravity_[1] = -9.81;
	}
	
	private:
		FieldMatrix<DT,n,n> outerK_;
		FieldMatrix<DT,n,n> innerK_;
		FieldVector<DT,n> outerLowerLeft_;
		FieldVector<DT,n> outerUpperRight_;
		FieldVector<DT,n> innerLowerLeft_;
		FieldVector<DT,n> innerUpperRight_;
		DT width_, height_;
		DT eps_;
		RT densityW_, densityN_;
		FieldVector<DT,n> gravity_;
		RT outerSwr_, outerSnr_, innerSwr_, innerSnr_;
		RT outerPorosity_, innerPorosity_;
		RT outerAlpha_, innerAlpha_;
		RT outerN_, innerN_;
  };

}
#endif
