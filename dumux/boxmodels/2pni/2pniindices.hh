// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 Bernd Flemisch                                       *
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
/*!
 * \file
 *
 * \brief Defines the indices used by the non-isotherm two-phase BOX model.
 */
#ifndef DUMUX_2PNI_INDICES_HH
#define DUMUX_2PNI_INDICES_HH

namespace Dumux
{
/*!
 * \addtogroup TwoPNIBoxModel
 */
// \{

/*!
 * \brief Enumerations for the non-isothermal two-phase model
 */
template <int PVOffset = 0>
class TwoPNIIndices : public TwoPIndices<PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 2; //! The primary variable index for temperature
    static const int energyEqIdx = PVOffset + 2; //! The equation index of the energy equation
    // Phase indices
       static const int wPhaseIdx = 0; //!< Index of the wetting phase in a phase vector
       static const int nPhaseIdx = 1; //!< Index of the non-wetting phase in a phase vector
};
// \}
}

#endif
