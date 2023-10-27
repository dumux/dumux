// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MineralizationModel
 * \brief Defines the properties required for compositional porous medium flow
 *        models considering mineralization processes of one or more of the
 *        components.
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consists of only a storage and a source term,
 * \f[
 * \frac{\partial ( \varrho_\lambda \phi_\lambda )} {\partial t} = q_\lambda,
 * \f]
 *
 * where:
 * * \f$ \varrho_\lambda \f$ is the mass density of the solid phase \f$ \lambda \f$,
 * * \f$ \phi_\lambda \f$ is the porosity of the solid,
 * * \f$ q_\lambda \f$ is a source or sink term.
 */

#ifndef DUMUX_MINERALIZATION_MODEL_HH
#define DUMUX_MINERALIZATION_MODEL_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup MineralizationModel
 * \brief Specifies a number properties of
 *        models that consider mineralization processes.
 *
 * \ţparam NonMinTraits traits class of the underlying model
 *                      not considering mineralization.
 * \tparam numPS number of solid phases to be considered.
 * \tparam numInertSP number of inert solid phases to be considered.
 */
template<class NonMinTraits, int numSC, int numInertSC>
struct MineralizationModelTraits : public NonMinTraits
{
    //! the number of mineral phases
    static constexpr int numSolidComps() { return numSC; }
     //! the number of inert mineral phases
    static constexpr int numInertSolidComps() { return numInertSC; }
    //! we additionally solve one equation per precipitating mineral phase
    static constexpr int numEq() { return NonMinTraits::numEq() + numSC - numInertSC; }
};
} // end namespace Dumux

#endif
