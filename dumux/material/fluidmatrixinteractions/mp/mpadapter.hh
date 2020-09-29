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

// remove from here after release 3.3 /////////////

#include <algorithm>
#include <cassert>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {


/*!
 * \ingroup Fluidmatrixinteractions
 * \brief An adapter for mpnc to use the capillary pressure-saturation relationships
 */
template <class MaterialLaw, int numFluidPhases>
class MPAdapter
{
    static_assert(AlwaysFalse<MaterialLaw>::value, "Adapter not implemented for the specified number of phases");
};

template <class MaterialLaw>
class [[deprecated("Use new material laws and FluidMatrix::MPAdapter instead!")]] MPAdapter<MaterialLaw, 2 /*numFluidPhases*/>
{
public:
    using Params = typename MaterialLaw::Params;
    /*!
     * \brief The capillary pressure-saturation curve.
     * \param values Container for the return values
     * \param params Array of parameters
     * \param state Fluidstate
     * \param wPhaseIdx the phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state,
                                   int wPhaseIdx)
    {
        assert(values.size() == 2);
        const int nPhaseIdx = 1 - wPhaseIdx;
        // nonwetting phase gets the capillary pressure added
        values[nPhaseIdx] = 0;
        // wetting phase does not get anything added
        values[wPhaseIdx] = - MaterialLaw::pc(params, state.saturation(wPhaseIdx));
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param values Container for the return values
     * \param params Array of parameters
     * \param state Fluidstate
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state,
                                       int wPhaseIdx)
    {
        assert(values.size() == 2);
        const int nPhaseIdx = 1 - wPhaseIdx;
        values[wPhaseIdx] = MaterialLaw::krw(params, state.saturation(wPhaseIdx));
        values[nPhaseIdx] = MaterialLaw::krn(params, state.saturation(wPhaseIdx));
    }
};


} // end namespace Dumux

// remove until here after release 3.3 /////////////

#include <dune/common/fvector.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief An adapter for mpnc to use the capillary pressure-saturation relationships
 */
template <class MaterialLaw, int numFluidPhases = MaterialLaw::numFluidPhases()>
class MPAdapter
{
    static_assert(AlwaysFalse<MaterialLaw>::value, "Adapter not implemented for the specified number of phases");
};


template<class MaterialLaw>
class MPAdapter<MaterialLaw, 2>
: public Adapter<MPAdapter<MaterialLaw, 2>, MultiPhasePcKrSw>
{
public:
    using Scalar = typename MaterialLaw::Scalar;

    MPAdapter(const MaterialLaw& pcKrS)
    : pcKrS_(pcKrS)
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
    const MaterialLaw& pcKrS_;
};

} // end namespace Dumux::FluidMatrix


#endif
