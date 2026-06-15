// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \copydoc Dumux::BedloadFormula
 */
#ifndef DUMUX_MATERIAL_BEDLOADFORMULA_HH
#define DUMUX_MATERIAL_BEDLOADFORMULA_HH

namespace Dumux {
/*!
 * \ingroup BedloadTransport
 * \brief Abstract base class for bedload transport formulas.
 */

template <typename VolumeVariables>
class BedloadFormula
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:

    /*!
     * \brief Compute the bedload transport rate.
     *
     * \param volVars Volume Variables.
     *
     * Compute the bedload transport rate in x- and y-direction.
     *
     * \return bedload transport rate [m^3/s/m]. First entry is the x-component, the second the y-component.
     */
    virtual Dune::FieldVector<Scalar, 2> bedloadDischarge(const VolumeVariables& volVars) const = 0;

    virtual ~BedloadFormula() {}
};

} // end namespace Dumux

#endif // DUMUX_MATERIAL_BEDLOADFORMULA_HH
