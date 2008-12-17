// $Id$

#ifndef DUNE_DIFFUSION_HH
#define DUNE_DIFFUSION_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

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
	- Scalar        type used for return values

   */
  template<class Grid, class Scalar,class VC>
  class Diffusion {
  public:
	FractionalFlowProblem<Grid, Scalar, VC>& diffproblem; //!< problem data

	//! \brief Calculate the pressure.
	/*!
	 *  \param saturation vector containing the saturation values
	 *  \param t time
	 *
	 *  Calculates the pressure \f$p\f$ as solution of the boundary value problem
	 *  \f[ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f]
	 *  subject to appropriate boundary and initial conditions.
	 */
	virtual void pressure(const Scalar t=0) = 0;

	//! \brief Calculate the total velocity.
	/*!
	 *  \param t time
	 *
	 *  \return block vector to store the total velocity
	 *
	 *  Given the piecewise constant pressure \f$p\f$ in form of the vector \a pressure,
	 *  this method calculates the total velocity according to the formula
	 *  \f$\boldsymbol{v}_\text{t} = - \lambda K \text{grad}\, p\f$.
	 *  The method is used in FractionalFlow to provide the velocity field required for the saturation equation.
	 */

	virtual void calcTotalVelocity(const Scalar t=0) const {
		return;
	}

	//! always define virtual destructor in abstract base class
	virtual ~Diffusion () {}

	//! without specification of a level, the class works on the leaf grid.
	/**
	 * \param grid grid object of type Grid
	 * \param prob a problem class object derived from DiffusionProblem
	*/
	Diffusion(const Grid& grid, FractionalFlowProblem<Grid, Scalar, VC>& prob)
	: grid(grid), diffproblem(prob), level_(grid.maxLevel())
	{
	}

	//! Constructor with possibility to specify a level for the diffusion class to work on.
	/**
	 * \param grid grid object of type Grid
	 * \param prob a problem class object derived from DiffusionProblem
	 * \param lev the grid level to work on
	 */
	Diffusion(const Grid& grid, FractionalFlowProblem<Grid, Scalar, VC>& prob, int lev)
	: diffproblem(prob), grid(grid), level_(lev)
	{
	}

	//! Returns the grid level on which the class works
	int level() const
	{
		return level_;
	}
	  const Grid& grid;
  protected:
	  const int level_;
  };

}
#endif
