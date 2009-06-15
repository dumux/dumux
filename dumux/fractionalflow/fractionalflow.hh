// $Id$
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
#ifndef DUNE_FRACTIONALFLOW_HH
#define DUNE_FRACTIONALFLOW_HH

#include "dumux/diffusion/diffusion.hh"
#include "dumux/transport/transport.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical multi-phase flow model
 * @author Bernd Flemisch, Markus Wolff
 */

/**
 * \defgroup fracflow Decoupled and Fractional Flow
 */
/*!
 * \ingroup fracflow
 * \defgroup impes IMPES (IMplicit Pressure Explicit Saturation)
 */

namespace Dune
{
/*!
 * \ingroup impes
 * \ingroup fracflow
 * \brief Base class for fractional flow formulation of multi-phase flow.
 * An interface to combine a diffusive model solving a pressure equation with a transport model solving a saturation equation
 * in the fractional flow context.
 *
 Template parameters are:

 - GridView      a DUNE gridview type
 - Diffusion     class defining the diffusion model
 - Transport     class defining the transport model
 - VC            type of a class containing different variables of the model
 */
template<class GridView, class Diffusion, class Transport, class VC>
class FractionalFlow
{
typedef    typename Diffusion::ScalarType Scalar;
public:
    typedef typename VC::ScalarVectorType RepresentationType;//!< Data type for a Vector of Scalars

    //! Set initial solution and initialize parameters
    virtual void initial() = 0;

    //! return const reference to saturation vector
    const RepresentationType& operator* () const
    {
        return transport.problem().variables().saturation();
    }

    //! return reference to saturation vector
    RepresentationType& operator* ()
    {
        return transport.problem().variables().saturation();
    }

    //! Calculate the update.
    /*!
     *  \param  t         time
     *  \param dt         time step size
     *  \param updateVec  vector for the update values
     *  \param CLFFac     security factor for the time step criterion (0 < CLFFac <= 1)
     *
     *  Calculates the update. Called from Dune::Timeloop.
     */
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar cFLFactor = 1) = 0;

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

    //! \brief Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    virtual void vtkout (const char* name, int k) const = 0;

    //! Constructs a FractionalFlow object
    /**
     * \param diff an object of type Diffusion
     * \param trans an object of type Transport
     */
    FractionalFlow (Diffusion& diff, Transport& trans)
    : diffusion(diff), transport(trans)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlow ()
    {}
protected:
    Diffusion& diffusion;//!< object of type Diffusion including the diffusion model
    Transport& transport;//!< object of type Transport indluding the transport model
};
}
#endif
