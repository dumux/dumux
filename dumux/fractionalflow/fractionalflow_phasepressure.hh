// $Id$

#ifndef DUNE_FRACTIONALFLOW_PHASEPRESSURE_HH
#define DUNE_FRACTIONALFLOW_PHASEPRESSURE_HH

#include "dumux/diffusion/diffusion_phasepressure.hh"
#include "dumux/transport/transport_phasepressure.hh"

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

template<class Grid, class Diffusion, class Transport, class VC>
class FractionalFlow  : public Transport, public Diffusion {
public:
    typedef typename VC::ScalarVectorType RepresentationType;
    typedef typename Diffusion::ScalarType Scalar;


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

    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar cFLFactor = 1) = 0;

    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    virtual void vtkout (const char* name, int k) const = 0;

    //! Construct a FractionalFlow object.
    FractionalFlow (Diffusion& diffusion, Transport& transport)
        : Transport(transport), Diffusion(diffusion)
    {
        if (transport.transProblem.variables.levelTransport > diffusion.diffProblem.variables.levelDiffusion)
            DUNE_THROW(Exception,"from class Twophase (or derived): transport class level is higher than diffusion class level!");
    }

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlow () {}

    const Grid& grid() const
    { return Transport::grid(); }

};
}
#endif
