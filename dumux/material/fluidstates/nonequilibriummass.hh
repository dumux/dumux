// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
