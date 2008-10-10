// $Id$ 

#ifndef DUNE_FRACTIONALFLOW_HH
#define DUNE_FRACTIONALFLOW_HH

#include "dumux/diffusion/diffusion.hh"
#include "dumux/transport/transport.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical two phase flow model
 * @author Bernd Flemisch, last changed by Markus Wolff
 */

/**
 * \defgroup fracflow Fractional Flow Formulation
 */

namespace Dune
{
	/*!\ingroup fracflow 
	 * \brief Standard two phase model. 
	 *
	 * This class implements the standard two phase model 
	 * for the pressure \f$p\f$ and the 
	 * wetting phase saturation \f$S\f$, namely,  
	 * \f{align*}
	 * - \text{div}\, (\lambda (S) K \text{grad}\, p ) &= 0, \\
	 * S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_t(p, S)) &= 0,
	 * \f}
	 * supplemented by appropriate initial and boundary conditions. 
	 */

  template<class G, class Diffusion, class Transport, class VC>
  class FractionalFlow  : public Transport, public Diffusion {
  public:
	  typedef typename Transport::RepresentationType RepresentationType;
	  typedef typename Diffusion::RepresentationType PressType;
	  typedef typename Diffusion::NumberType RT;
	  
//	  Diffusion& diffusion;
  
	//! \brief Calculate the pressure.
	/*!
	 *  \param t time 
	 *  
	 *  Calculates the pressure \f$p\f$ as solution of the boundary value problem 
	 *  \f[ - \text{div}\, (\lambda(S) K \text{grad}\, p ) = 0, \f]
	 *  subject to appropriate boundary and initial conditions. 
	 *  Employ the method \a pressure of Diffusion. 
	 */
//	void pressure(const RT t=0)
//	{
//		Diffusion::pressure(t);
//	}
	
	  
		virtual void initial() = 0;
		
		//! return const reference to saturation vector
		const RepresentationType& operator* () const
		{
		  return this->transproblem.variables.saturation;
		}

		//! return reference to saturation vector
		RepresentationType& operator* ()
		{
		  return this->transproblem.variables.saturation;
		}

	 
	//! \brief Calculate the total velocity.
	/*!
	 *  \param t time 
	 * 
	 *  Given the piecewise constant pressure \f$p\f$ in form of the vector \a pressure, 
	 *  this method calculates the total velocity according to the formula 
	 *  \f$\boldsymbol{v}_\text{t} = - \lambda(S) K \text{grad}\, p\f$. 
	 *  Employ the method \a totalVelocity of Diffusion. 
	 */
	virtual void totalVelocity(const RT t=0) = 0;
	
	//! \brief Calculate the update vector.
	/*!
	 *  \param[in]  t         time 
	 *  \param[out] dt        time step size
	 *  \param[out] updateVec vector for hte update values
	 *  
	 *  Calculate the update vector, i.e., the discretization 
	 *  of \f$\text{div}\, (f_\text{w}(S) \boldsymbol{v}_t)\f$.
	 */
	virtual int update(const RT t, RT& dt, RepresentationType& updateVec, RT cFLFactor = 1) = 0;  
	
	virtual void vtkout (const char* name, int k) const = 0;
		
	//! Construct a FractionalFlow object.
	FractionalFlow (Diffusion& diff, Transport& trans)
            : Transport(trans), Diffusion(diff)
	{ 
		if (trans.level() > diff.level()) 
		  DUNE_THROW(Exception,"from class Twophase (or derived): transport class level is higher than diffusion class level!");
	}

	//! always define virtual destructor in abstract base class
	virtual ~FractionalFlow () {}
  };
}
#endif
