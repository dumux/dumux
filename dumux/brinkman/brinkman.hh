#ifndef DUNE_BRINKMAN_HH
#define DUNE_BRINKMAN_HH

#include "dumux/brinkman/brinkmanproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * 
 * \defgroup diffusion Brinkman
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
	- VelType   type of the vector holding the velocity values 

   */
  template<class G, class RT, class RepresentationType, class VelType>
  class Brinkman {
  public:
	RepresentationType pressure; //!< vector of pressure values
	RepresentationType pressureCorrection; //!< vector of pressure correction values
	VelType velocity;
	VelType velocityCorrection;
	BrinkmanProblem<G, RT>& problem; //!< problem data
	typedef RT NumberType;
	
	virtual void computeVelocity() = 0;
	
	void SIMPLE() 
	{
		double tolerance = 1e-5; 
		int maxIter = 10;
		double error = 1e100;
		int iter = 0;
		while (error > tolerance && iter <= maxIter)
		{
			iter++;
			computeVelocity();
			//computePressureCorrection();
			//computeVelocityCorrection();
			//correctVelocityAndPressure();
			error = pressureCorrection.two_norm()/pressure.two_norm();
		}
	}
	
	//! generate vtk output
	virtual void vtkout (const char* name, int k) const = 0;
	
	//! return const reference to pressure vector
	const RepresentationType& operator* () const
	{
	  return pressure;
	}

	//! return reference to permeability vector
	RepresentationType& operator* ()
	{
	  return pressure;
	}

	//! always define virtual destructor in abstract base class
	virtual ~Brinkman () {}
	
	//! without specification of a level, the class works on the leaf grid.
	/**
	 * \param g grid object of type G
	 * \param prob a problem class object derived from BrinkmanProblem
	*/
	Brinkman(const G& g, BrinkmanProblem<G, RT>& prob) 
	: grid(g), problem(prob)
	{ 
	}
	
  protected:
	  const G& grid;
  };

}
#endif
