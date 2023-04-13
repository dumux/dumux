// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Makes the capillary pressure-saturation relations available under the M-phase API for material laws
 *
 * Makes the capillary pressure-saturation relations
 * available under the M-phase API for material laws
 */
#ifndef DUMUX_MP_ADAPTER_HH
#define DUMUX_MP_ADAPTER_HH

#include <dune/common/fvector.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief An adapter for mpnc to use the capillary pressure-saturation relationships
 */
template <class MaterialLaw, int numFluidPhases = std::decay_t<MaterialLaw>::numFluidPhases()>
class MPAdapter
{
    static_assert(AlwaysFalse<MaterialLaw>::value, "Adapter not implemented for the specified number of phases");
};


template<class MaterialLaw>
class MPAdapter<MaterialLaw, 2>
: public Adapter<MPAdapter<MaterialLaw, 2>, MultiPhasePcKrSw>
{
public:
    using Scalar = typename std::decay_t<MaterialLaw>::Scalar;

    MPAdapter(MaterialLaw&& pcKrS)
    : pcKrS_(std::forward<MaterialLaw>(pcKrS))
    {}

    /*!
     * \brief The capillary pressure-saturation curve.
     * \param state Fluidstate
     * \param wPhaseIdx the phase index of the wetting phase
     */
    template <class FluidState>
    auto capillaryPressures(const FluidState& state, int wPhaseIdx) const
    {
        Dune::FieldVector<typename FluidState::Scalar, 2> values;

        const int nPhaseIdx = 1 - wPhaseIdx;
        // non-wetting phase gets the capillary pressure added
        values[nPhaseIdx] = 0;
        // wetting phase does not get anything added
        values[wPhaseIdx] = - pcKrS_.pc(state.saturation(wPhaseIdx));

        return values;
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param state Fluidstate
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class FluidState>
    auto relativePermeabilities(const FluidState& state, int wPhaseIdx) const
    {
        Dune::FieldVector<typename FluidState::Scalar, 2> values;

        const int nPhaseIdx = 1 - wPhaseIdx;
        values[wPhaseIdx] = pcKrS_.krw(state.saturation(wPhaseIdx));
        values[nPhaseIdx] = pcKrS_.krn(state.saturation(wPhaseIdx));

        return values;
    }
private:
    MaterialLaw pcKrS_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the MPAdapter class.
 *        Makes sure that MPAdapter stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
MPAdapter(T&&) -> MPAdapter<T>;

} // end namespace Dumux::FluidMatrix


#endif
