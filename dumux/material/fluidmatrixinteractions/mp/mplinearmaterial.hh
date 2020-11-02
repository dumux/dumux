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
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 */
#ifndef DUMUX_MP_LINEAR_MATERIAL_HH
#define DUMUX_MP_LINEAR_MATERIAL_HH

#include "mplinearmaterialparams.hh"

#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 * \sa MpLinearMaterialParams
 */
template <int numPhasesV, class ScalarT, class ParamsT = MpLinearMaterialParams<numPhasesV, ScalarT> >
class [[deprecated("Use FluidMatrix::MPLinearMaterial. Will be removed after 3.3")]] MpLinearMaterial
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;
    enum { numPhases = numPhasesV };

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param values Container for the return values
     * \param params Array of Parameters
     * \param state The fluid state
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state,
                                   int wPhaseIdx = 0)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar S = state.saturation(phaseIdx);
            values[phaseIdx] =
                S*params.pcMaxSat(phaseIdx) +
                (1 - S)*params.pcMinSat(phaseIdx);
        }
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param values Container for the return values
     * \param params Array of Parameters
     * \param state The fluid state
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state,
                                       int wPhaseIdx = 0)
    {
        using std::clamp;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            values[phaseIdx] = clamp(state.saturation(phaseIdx), 0.0, 1.0);
    }
};
} // end namespace Dumux

#include <dune/common/fvector.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 */
template<class S, int numFluidPhases>
class MPLinearMaterial
: public Adapter<MPLinearMaterial<S, numFluidPhases>, MultiPhasePcKrSw>
{
    struct Params
    {
        Params(const std::array<S, numFluidPhases>& pcMaxSat,
               const std::array<S, numFluidPhases>& pcMinSat)
        : pcMaxSat_(pcMaxSat), pcMinSat_(pcMinSat) {}

        S pcMaxSat(int phaseIdx) const { return pcMaxSat_[phaseIdx]; }
        S pcMinSat(int phaseIdx) const { return pcMinSat_[phaseIdx]; }
    private:
        std::array<S, numFluidPhases> pcMaxSat_;
        std::array<S, numFluidPhases> pcMinSat_;
    };

public:
    using BasicParams = Params;
    using Scalar = S;

    MPLinearMaterial(const BasicParams& basicParams)
    : basicParams_(basicParams)
    {}

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param state The fluid state
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class FluidState>
    auto capillaryPressures(const FluidState& state, int wPhaseIdx = 0) const
    {
        static_assert(FluidState::numPhases == numFluidPhases, "FluidState doesn't match the number of fluid phases!");
        Dune::FieldVector<typename FluidState::Scalar, numFluidPhases> values;
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
        {
            const Scalar saturation = state.saturation(phaseIdx);
            values[phaseIdx] =
                saturation*basicParams_.pcMaxSat(phaseIdx) +
                (1 - saturation)*basicParams_.pcMinSat(phaseIdx);
        }
        return values;
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param state The fluid state
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class FluidState>
    auto relativePermeabilities(const FluidState& state, int wPhaseIdx = 0) const
    {
        static_assert(FluidState::numPhases == numFluidPhases, "FluidState doesn't match the number of fluid phases!");
        using std::clamp;
        Dune::FieldVector<typename FluidState::Scalar, FluidState::numPhases> values;
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            values[phaseIdx] = clamp(state.saturation(phaseIdx), 0.0, 1.0);
        return values;
    }
private:
    BasicParams basicParams_;
};

} // end namespace Dumux::FluidMatrix

#endif
