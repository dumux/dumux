#ifndef DUNE_DIFFUSION_HH
#define DUNE_DIFFUSION_HH

#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * 
 * \defgroup diffusion Diffusion
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Base class for defining an instance of a numerical diffusion model.
  /*! An interface for defining a numerical diffusion model for the 
   *  solution of equations of the form 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation. 
	Template parameters are:

	- Grid      a DUNE grid type
	- RT        type used for return values 
	- RepresentationType type of the vector holding the pressure values 
	- SatType   type of the vector holding the saturation values 
	- VelType   type of the vector holding the velocity values 

   */
  template<class G, class RT, class RepresentationType, class SatType, class VelType>
  class Diffusion {
  public:
	RepresentationType press; //!< vector of pressure values
	DiffusionProblem<G, RT>& problem; //!< problem data
	typedef RT NumberType;
	
	//! \brief Calculate the pressure.
	/*!
	 *  \param saturation vector containing the saturation values 
	 *  \param t time 
	 *  
	 *  Calculates the pressure \f$p\f$ as solution of the boundary value problem 
	 *  \f[ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f]
	 *  subject to appropriate boundary and initial conditions. 
	 */
	virtual void pressure(const SatType& saturation=0, const RT t=0) = 0;  
	
	//! \brief Calculate the total velocity.
	/*!
	 *  \param saturation vector containing the saturation values 
	 *  \param t time 
	 * 
	 *  \return block vector to store the total velocity 
	 *  
	 *  Given the piecewise constant pressure \f$p\f$ in form of the vector \a pressure, 
	 *  this method calculates the total velocity according to the formula 
	 *  \f$\boldsymbol{v}_\text{t} = - \lambda K \text{grad}\, p\f$. 
	 *  The method is used in FractionalFlow to provide the velocity field required for the saturation equation. 
	 */
	virtual void totalVelocity(VelType& velocity, const SatType& saturation=0, const RT t=0) const = 0;  
	
	//! generate vtk output
	virtual void vtkout (const char* name, int k) const = 0;
	
	//! return const reference to pressure vector
	const RepresentationType& operator* () const
	{
	  return press;
	}

	//! return reference to permeability vector
	RepresentationType& operator* ()
	{
	  return press;
	}

	//! always define virtual destructor in abstract base class
	virtual ~Diffusion () {}
	
	Diffusion(const G& g, DiffusionProblem<G, RT>& prob) 
	: grid(g), problem(prob), level_(g.maxLevel())
	{ 
	}
	Diffusion(const G& g, DiffusionProblem<G, RT>& prob, int lev) 
	: problem(prob), grid(g), level_(lev)
	{ 
	}
	
	int level() const
	{
		return level_;
	}
	
  protected:
	  const G& grid;
	  const int level_;
  };

}
#endif
