// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Bernd Flemisch, Jochen Fritz                      *
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

#ifndef NUMERICALFLUX_HH
#define NUMERICALFLUX_HH

namespace Dune
{


/*! \ingroup transport
 *  \brief numerical flux.
 *
 *  Given a specific fractional flow function \f$f\f$,
 *  this class is responsible for calculating the numerical flux required in Transport
 *  for the discretization of the scalar transport equation.
 */
template<class Scalar>
class NumericalFlux {
public:
    /*! \brief Realizes the numerical flux function.
     *
     *  \param Sa first argument, usually the saturation value from one cell
     *  \param Sb second argument, usually the saturation value from another cell
     *  \param fa fractional flow function evaluated at \a Sa
     *  \param fb fractional flow function evaluated at \a Sb
     *
     *  \return Given a specific fractional flow function \f$f\f$,
     *  this numerical flux function returns
     *  - \f$\min (f(a), f(b))\f$ if \f$a < b\f$,
     *  - \f$\max (f(a), f(b))\f$ else.
     */
    virtual Scalar operator() (Scalar Sa, Scalar Sb, Scalar fa, Scalar fb) const = 0;

    virtual ~NumericalFlux()
    { }
};


//! Upwind numerical flux.
template<class Scalar>
class Upwind : public NumericalFlux<Scalar> {
public:
    /*!
     *  \return Given a specific fractional flow function \f$f\f$,
     *  this numerical flux function returns
     *  - \f$\min (f(a), f(b))\f$ if \f$a < b\f$,
     *  - \f$\max (f(a), f(b))\f$ else.
     *  This results in an approximative Godunov scheme, which coincides with the upwind scheme in case
     *  of a monotone function \f$f\f$.
     */
    Scalar operator() (Scalar Sa, Scalar Sb, Scalar fa, Scalar fb) const
    {
        return Sa < Sb ? std::min(fa, fb) : std::max(fa, fb);
    }
};

//! \brief Lax-Friedrichs numerical flux.
template<class Scalar>
class LaxFriedrichs : public NumericalFlux<Scalar> {
public:
    /*!
     *  \return Given a specific fractional flow function \f$f\f$,
     *  this numerical flux function returns \f$\frac 12 (f(a) + f(b)) + \frac 12 (a - b)\f$.
     */
    Scalar operator() (Scalar Sa, Scalar Sb, Scalar fa, Scalar fb) const
    {
        return 0.5*(fa + fb) + 0.5*(Sa - Sb);
    }
};
}

#endif

