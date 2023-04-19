// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief The Kozeny-Carman relationship for the calculation of a porosity-dependent permeability.
 */
#ifndef DUMUX_PERMEABILITY_KOZENY_CARMAN_HH
#define DUMUX_PERMEABILITY_KOZENY_CARMAN_HH

#include <cmath>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The Kozeny-Carman relationship for the calculation of a porosity-dependent permeability.
 *        When the porosity is implemented as solution-independent, using this relationship for the
 *        permeability leads to unnecessary overhead.
 *
 * \tparam PermeabilityType The type used for the intrinsic permeability
 */
template<class PermeabilityType>
class PermeabilityKozenyCarman
{
public:
    /*!
     * \brief Calculates the permeability for a given sub-control volume
     * \param refPerm Reference permeability before porosity changes
     * \param refPoro The poro corresponding to the reference permeability
     * \param poro The porosity for which permeability is to be evaluated
     */
    template<class Scalar>
    PermeabilityType evaluatePermeability(PermeabilityType refPerm, Scalar refPoro, Scalar poro) const
    {
        using Dune::power;
        auto factor = power((1.0 - refPoro)/(1.0 - poro), 2) * power(poro/refPoro, 3);
        refPerm *= factor;
        return refPerm;
    }
};

} // namespace Dumux

#endif
