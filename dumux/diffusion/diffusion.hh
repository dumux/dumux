// $Id:$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_DIFFUSION_HH
#define DUNE_DIFFUSION_HH

#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch, Markus Wolff
 */

/*!
 * \ingroup fracflow
 * \defgroup diffusion Diffusion
 */

namespace Dune
{
//! \ingroup diffusion
//! Base class for defining an instance of a numerical diffusion model.
/*! An interface for defining a numerical diffusion model for the
 *  solution of equations of the form
 *  \f[\text{div}\, \boldsymbol{v} = q,\f]
 *  where, the velocity \f$\boldsymbol{v} \sim \boldsymbol{K} \nabla p \f$,
 *  \f$p\f$ is a pressure and \f$q\f$ a source/sink term

 Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem
*/
template<class GridView, class Scalar, class VC, class Problem = DiffusionProblem<GridView, Scalar, VC> >
class Diffusion
{
public:
    typedef Scalar ScalarType;

    //! Calculate the pressure.
    /*!
     *  \param t time
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value problem
     *  \f[  \text{div}\, \boldsymbol{v} = q, \f]
     *  subject to appropriate boundary conditions.
     */
    virtual void pressure(bool first, const Scalar t = 0)
    {
        return;
    }

    //! Calculate the velocity.
    /*!
     *  \param t time
     *
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity
     *  The method is needed in the IMPES (Implicit Pressure Explicit Saturation) algorithm which is used for a fractional flow formulation
     *  to provide the velocity field required for the solution of the saturation equation.
     */
    virtual void calculateVelocity(const Scalar t = 0) const
    {
        return;
    }

    //! Start a post-processing procedure at the end of a timestep.
    /*!
     *  \param t time
     *  \param dt time step
     *
     *  If an explicit <it>Euler</it> time discretization is used this function will be called
     *  at the end of each time step.
     */
    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    //! Returns a reference to the problem
    virtual Problem& problem()
    {
        return diffProblem;
    }


    //! always define virtual destructor in abstract base class
    virtual ~Diffusion()
    {
    }

    //! Constructs a Diffusion object
    /**
     * \param grid gView object of type GridView
     * \param prob a problem class object
     */
    Diffusion(GridView& gView, Problem& prob) :
        gridView(gView), diffProblem(prob)
    {
    }

protected:
    const GridView& gridView; //!< object of type Dune::GridView
    Problem& diffProblem; //!< problem data
};

}
#endif
