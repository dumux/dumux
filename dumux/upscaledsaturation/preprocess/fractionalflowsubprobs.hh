// $Id$
#ifndef DUNE_FRACTIONALFLOWSUBPROBS_HH
#define DUNE_FRACTIONALFLOWSUBPROBS_HH

#include "dumux/upscaledsaturation/preprocess/diffusionsubprobs.hh"
#include "dumux/upscaledsaturation/preprocess/transportsubprobs.hh"

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

  template<class Grid, class DiffusionSubProbs, class TransportSubProbs, class VC>
  class FractionalFlowSubProbs  : public TransportSubProbs, public DiffusionSubProbs {
  public:
      typedef typename TransportSubProbs::RepresentationType RepresentationType;
      typedef typename VC::ScalarVectorType PressType;
      typedef typename DiffusionSubProbs::ScalarType Scalar;

//      Diffusion& diffusion;

    //! \brief Calculate the pressure.
    /*!
     *  \param t time
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value problem
     *  \f[ - \text{div}\, (\lambda(S) K \text{grad}\, p ) = 0, \f]
     *  subject to appropriate boundary and initial conditions.
     *  Employ the method \a pressure of Diffusion.
     */

        virtual void initial() = 0;

        //! return const reference to saturation vector
        const RepresentationType& operator* () const
        {
          return this->transProblem.variables.saturation;
        }

        //! return reference to saturation vector
        RepresentationType& operator* ()
        {
          return this->transProblem.variables.saturation;
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
    virtual void totalVelocity(const Scalar t=0) = 0;

    //! \brief Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[out] dt        time step size
     *  \param[out] updateVec vector for hte update values
     *
     *  Calculate the update vector, i.e., the discretization
     *  of \f$\text{div}\, (f_\text{w}(S) \boldsymbol{v}_t)\f$.
     */
    virtual void update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar cFLFactor = 1) = 0;


    //! Construct a FractionalFlow object.
    FractionalFlowSubProbs (DiffusionSubProbs& diffusionProblem, TransportSubProbs& transportProblem)
    : DiffusionSubProbs(diffusionProblem), TransportSubProbs(transportProblem)
    {
        if (transportProblem.level() > diffusionProblem.level())
          DUNE_THROW(Exception,"from class Twophase (or derived): transport class level is higher than diffusion class level!");
    }

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlowSubProbs () {}
  };
}
#endif
