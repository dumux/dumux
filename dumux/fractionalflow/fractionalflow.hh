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
 * @brief  Base class for defining an instance of a numerical two phase flow model
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

template<class GridView, class Diffusion, class Transport, class VC>
class FractionalFlow
{
public:
typedef    typename VC::ScalarVectorType RepresentationType;
    typedef typename Diffusion::ScalarType Scalar;

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
    FractionalFlow (Diffusion& diff, Transport& trans)
    : diffusion(diff), transport(trans)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlow ()
    {}
protected:
    Diffusion& diffusion;
    Transport& transport;
};
}
#endif
