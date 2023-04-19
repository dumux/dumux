// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
#ifndef DUMUX_NONEQUILIBRIUM_MASS_FLUID_STATE_HH
#define DUMUX_NONEQUILIBRIUM_MASS_FLUID_STATE_HH

#include <dune/common/exceptions.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system
 *
 *        This fluidstate ought to be used for the case of:
 *        - local thermal equilibrium
 *        - local chemical non-equilibrium
 */
template <class ScalarType, class FluidSystem>
class NonEquilibriumMassFluidState
: public NonEquilibriumFluidState<ScalarType, FluidSystem>
{
    using ParentType = NonEquilibriumFluidState<ScalarType, FluidSystem>;
public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    //! export the scalar type
    using Scalar = ScalarType;

    using ParentType::ParentType;

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    using ParentType::setTemperature;
    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
     * \note only the setTemperature(value) overload in the NonEquilibriumFluidState makes sense for thermal equilibrium
     */
    void setTemperature(const int phaseIdx, const Scalar value)
    { DUNE_THROW(Dune::NotImplemented, "This is a fluidstate for *chemical* non-equilibrium, not thermal!"); }
};

} // end namespace Dumux

#endif
