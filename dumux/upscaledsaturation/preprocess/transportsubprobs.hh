// $Id$
#ifndef DUNE_TRANSPORTSUBPROBS_HH
#define DUNE_TRANSPORTSUBPROBS_HH

#include "dumux/upscaledsaturation/preprocess/fractionalflowproblemsubprobs.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * \defgroup transport Transport
 */

namespace Dune
{
  //! \ingroup transport
  //! Base class for defining an instance of a numerical transport model.
  /*! An interface for defining a numerical transport model for the
   *  solution of equations of the form
   *  \f$S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_\text{total}) = 0\f$,
   * \f$S = g\f$ on \f$\Gamma_1\f$, and \f$S(t = 0) = S_0\f$. Here,
   * \f$S\f$ denotes the wetting phase saturation,
   * \f$\boldsymbol{v}_\text{total}\f$ the total velocity,
   * and \f$f_\text{w}\f$ the wetting phase fractional flow function.

	- Grid      a DUNE grid type
	- RT        type used for return values
	- RepresentationType   type of the vector holding the saturation values
	- VelType   type of the vector holding the velocity values

   */
  template<class G, class RT, class VC>
  class TransportSubProbs {
  public:

	typedef typename VC::ScalarType RepresentationType;
	FractionalFlowProblemSubProbs<G, RT, VC>& transproblem; //!< problem data

	//! \brief Calculate the update vector.
	/*!
	 *  \param[in]  t         time
	 *  \param[out] dt        time step size
	 *  \param[out] updateVec vector for hte update values
	 *
	 *  Calculate the update vector, i.e., the discretization
	 *  of \f$\text{div}\, (f_\text{w}(S) \boldsymbol{v}_t)\f$.
	 */
	virtual int update(const RT t, RT& dt, RepresentationType& updateVec, RT& CLFFac) = 0;

	void initial()
	{
		initialTransport();
		return;
	}

	//! \brief Sets the initial solution \f$S_0\f$.
	virtual void initialTransport() = 0;

	//! return const reference to saturation vector
	virtual const RepresentationType& operator* () const
	{
	  return transproblem.variables.saturation;
	}

	//! return reference to saturation vector
	virtual RepresentationType& operator* ()
	{
	  return transproblem.variables.saturation;
	}

	//! always define virtual destructor in abstract base class
	virtual ~TransportSubProbs () {}

	/*! @brief constructor
	 *
	 *  This constructor gives the additional possibility to specify a grid level on which
	 *  the transport equation shall be solved. This especially important for multiscale methods.
	 *  @param g a DUNE grid object
	 *  @param prob an object of class TransportProblem or derived
	 *  @param lev the grid level on which the Transport equation is to be solved.
	 */
	TransportSubProbs(const G& g, FractionalFlowProblemSubProbs<G, RT, VC>& prob, int lev)
	: transproblem(prob), grid(g), level_(lev)
	{ }

	//! returns the level on which the transport eqution is solved.
	int level() const
	{
		return level_;
	}

	const G& grid;

  protected:
	  const int level_;
  };

}
#endif
