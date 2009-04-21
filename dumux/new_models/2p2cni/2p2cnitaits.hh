/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUMUX_2P2CNITRAITS_HH
#define DUMUX_2P2CNITRAITS_HH

#include <dumux/new_models/2p2c/2p2ctraits.hh>

namespace Dune
{
/*!
 * \brief Traits for the non-isothermal 2-phase 2-component model
 */
template <class Scalar,
          class BaseTraits = TwoPTwoCPwSnTraits<Scalar> >
class TwoPTwoCNITraits : public BaseTraits
{
public:
    static const int numEq = 3;  //!< Override the mumber of primary variables: We also have temperature
    // Primary variable indices
    static const int temperatureIdx = 2; //! The index for temperature in solution vectors.
};

}

#endif
