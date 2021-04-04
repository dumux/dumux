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
