// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
#ifndef DUMUX_NONEQUILIBRIUM_ENERGY_FLUID_STATE_HH
#define DUMUX_NONEQUILIBRIUM_ENERGY_FLUID_STATE_HH

#include <dumux/common/valgrind.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/common/propertysystem.hh>

#include <cmath>
#include <algorithm>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(FluidSystem);
}

/*!
 * \ingroup FluidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system
 *
 *        This fluidstate ought to be used for the case of:
 *        - local thermal non-equilibrium
 *        - local chemical equilibrium
 */
template <class TypeTag>
class NonEquilibriumEnergyFluidState :
    public NonEquilibriumFluidState<typename GET_PROP_TYPE(TypeTag, Scalar),
                                    typename GET_PROP_TYPE(TypeTag, FluidSystem)>
    {
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
        enum { numPhases = FluidSystem::numPhases };
        enum { numComponents = FluidSystem::numComponents };

    NonEquilibriumEnergyFluidState()
        : NonEquilibriumFluidState<Scalar, FluidSystem>()
    {
    }

    /*****************************************************
     * Access to fluid properties which only make sense
     * if assuming chemical equilibrium
     *****************************************************/
    /*!
     * \brief The fugacity of a component
     *
     * This assumes chemical equilibrium.
     */
    Scalar fugacity(int compIdx) const
    { return fugacity(0, compIdx); }

    /*!
     * \brief The fugacity of a component in a phase \f$\mathrm{[Pa]}\f$
    */
    Scalar fugacity(int phaseIdx, int compIdx) const
    {
        // Unfortunately throw does not work when triggered from a constructor
        std::cout <<"file: "<< __FILE__ << ", line: " << __LINE__ <<". This is a fluidstate for *thermal* non-equilibrium, not chemical! \n ";
        DUNE_THROW(Dune::NotImplemented, "This is a fluidstate for *thermal* non-equilibrium, not chemical!");
        return 0.;
    }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
     *        Both versions of the function need to be here.
     *        Otherwise the compiler gets confused.
     *        Thus, this is just forwarding to the Parent
     *        (unclear why this is necessary).
     */
    void setTemperature(int phaseIdx, Scalar value)
    {
        this->setTemperature(phaseIdx, value);
    }

    /*!
     * \brief Set the temperature \f$\mathrm{[K]}\f$ of a fluid phase
     *        Both versions of the function need to be here.
     *        Otherwise the compiler gets confused.
     *        Thus, this is just presenting the signature to the compiler.
     */
    void setTemperature(Scalar value)
    {
        // Unfortunately throw does not work when triggered from a constructor
        std::cout <<"file: "<< __FILE__ << ", line: " << __LINE__ <<". This is a fluidstate for *thermal* non-equilibrium, not chemical! \n ";
        DUNE_THROW(Dune::NotImplemented, "This is a fluidstate for *thermal* non-equilibrium, not chemical!");
    }

};

} // end namespace Dumux

#endif
