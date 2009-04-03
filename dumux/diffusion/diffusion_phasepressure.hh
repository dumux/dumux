// $Id$

#ifndef DUNE_DIFFUSION_PHASEPRESSURE_HH
#define DUNE_DIFFUSION_PHASEPRESSURE_HH

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
 - Scalar        type used for return values

*/
template<class Grid, class Scalar, class VC, class Problem = DiffusionProblem<Grid, Scalar, VC> >
class Diffusion
{
public:
    typedef Scalar ScalarType;

    //! \brief Calculate the pressure.
    /*!
     *  \param saturation vector containing the saturation values
     *  \param t time
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value problem
     *  \f[ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f]
     *  subject to appropriate boundary and initial conditions.
     */
    virtual void pressure(const Scalar t = 0)
    {
        return;
    }

    virtual void pressure(bool first, const Scalar t = 0)
    {
        return;
    }

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

    virtual void calcTotalVelocity(const Scalar t = 0) const
    {
        return;
    }

    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    //! always define virtual destructor in abstract base class
    virtual ~Diffusion()
    {
    }

    //! without specification of a level, the class works on the leaf grid.
    /**
     * \param grid grid object of type Grid
     * \param prob a problem class object derived from DiffusionProblem
     */
    Diffusion(const Grid& grid, Problem& prob) :
        grid(grid), diffProblem(prob)
    {
    }

    const Grid& grid;
    Problem& diffProblem; //!< problem data
};

}
#endif
