// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \brief @copydoc Dumux::FluidSystems::H2ON2Kinetic
 */
#ifndef DUMUX_H2O_N2_FLUID_SYSTEM_KINETIC_HH
#define DUMUX_H2O_N2_FLUID_SYSTEM_KINETIC_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/idealgas.hh>

#include <dumux/material/binarycoefficients/h2o_n2.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup FluidSystems
 * \brief A two-phase fluid system with two components water \f$(\mathrm{H_2O})\f$
 *        Nitrogen \f$(\mathrm{N_2})\f$ for non-equilibrium models. TODO: Is this fluid system necessary??
 */
template <class Scalar, class Policy = H2ON2DefaultPolicy<>>
class H2ON2Kinetic :
    public FluidSystems::H2ON2<Scalar, Policy>
{
private:
    using ParentType = FluidSystems::H2ON2<Scalar, Policy>;

    using IdealGas = Dumux::IdealGas<Scalar>;
public:
    //! The type of parameter cache objects
    using ParameterCache = NullParameterCache;

    /*!
     * \brief Return the enthalpy of a component in a phase.
     * \param fluidState A container with the current (physical) state of the fluid
     * \param phaseIdx The index of the phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar componentEnthalpy(FluidState &fluidState,
                                    const int phaseIdx,
                                    const int compIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        switch (phaseIdx) {
            case ParentType::liquidPhaseIdx:
                switch(compIdx) {
                    case ParentType::H2OIdx:
                        return ParentType::H2O::liquidEnthalpy(T, p);
                    case ParentType::N2Idx:
                        return ParentType::N2::gasEnthalpy(T, p); // TODO: should be liquid enthalpy
                    default:
                        DUNE_THROW(Dune::NotImplemented, "wrong index");
                }
            case ParentType::gasPhaseIdx:
                switch(compIdx) {
                    case ParentType::H2OIdx:
                        return ParentType::H2O::gasEnthalpy(T, p);
                    case ParentType::N2Idx:
                        return ParentType::N2::gasEnthalpy(T, p);
                    default:
                        DUNE_THROW(Dune::NotImplemented, "wrong index");
                }
            default:
                DUNE_THROW(Dune::NotImplemented, "wrong index");
        }
    }

    /*!
     * \brief Return the Henry constant for a component in a phase. \f$\mathrm{[Pa]}\f$
     * \param temperature The given temperature
     */
    static Scalar henry(Scalar temperature)
    {
        return BinaryCoeff::H2O_N2::henry(temperature);
    }

    /*!
     * \brief Return the vapor pressure of a component above one phase. \f$\mathrm{[Pa]}\f$
     * \param temperature The given temperature
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        return ParentType::H2O::vaporPressure(temperature);
    }
};
} // end namespace Fluidsystem
} // end namespace Dumux

#endif
