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
 * \brief This file contains the diffusion module for the vertex data
 *        of the fully coupled MpNc model
 */
#ifndef DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH

#include <dumux/common/valgrind.hh>
#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

namespace Dumux {

/*!
 * \brief Variables for the diffusive fluxes in the MpNc model within
 *        a finite volume.
 */
template<class TypeTag, bool enableDiffusion>
class MPNCVolumeVariablesDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };

public:
    /*!
     * \brief The constructor
     */
    MPNCVolumeVariablesDiffusion()
    {}
    /*!
     * \brief update
     *
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param volVars The volume variables
     * \param problem  The problem
     */
    void update(FluidState &fluidState,
                ParameterCache &paramCache,
                const VolumeVariables &volVars,
                const Problem &problem)
    {
        Valgrind::SetUndefined(*this);

        // diffusion coefficents in liquid
        diffCoeffL_[0] = 0.0;
        for (int compIdx = 1; compIdx < numComponents; ++compIdx) {
            diffCoeffL_[compIdx] =
                FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                        paramCache,
                                                        wPhaseIdx,
                                                        0,
                                                        compIdx);
        }
        Valgrind::CheckDefined(diffCoeffL_);

        // diffusion coefficents in gas
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            diffCoeffG_[compIIdx][compIIdx] = 0;
            for (int compJIdx = compIIdx + 1; compJIdx < numComponents; ++compJIdx) {
                diffCoeffG_[compIIdx][compJIdx] =
                        FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                                paramCache,
                                                                nPhaseIdx,
                                                                compIIdx,
                                                                compJIdx);

                // fill the symmetric part of the diffusion coefficent
                // matrix
                diffCoeffG_[compJIdx][compIIdx] = diffCoeffG_[compIIdx][compJIdx];
            }
        }
        Valgrind::CheckDefined(diffCoeffG_);
    }

    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     * \param phaseIdx The local index of the phases
     * \param compIIdx The local index of the first component in the phase
     * \param compJIdx The local index of the second component in the phase
     */
    Scalar diffCoeff(const unsigned int phaseIdx,
                     const unsigned int compIIdx,
                     const unsigned int compJIdx) const
    {
        if (phaseIdx == nPhaseIdx)
            // TODO: tensorial diffusion coefficients
            return diffCoeffG_[compIIdx][compJIdx];

        using std::max;
        using std::min;
        const unsigned int i = min(compIIdx, compJIdx);
        const unsigned int j = max(compIIdx, compJIdx);
        if (i != 0)
            return 0;
        return diffCoeffL_[j];
    }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(diffCoeffL_);
        Valgrind::CheckDefined(diffCoeffG_);
    }


protected:
    // the diffusion coefficients for the porous medium for the
    // liquid phase
    Scalar diffCoeffL_[numComponents];

    // the diffusion coefficients for the porous medium for the
    // gas phase
    Scalar diffCoeffG_[numComponents][numComponents];
};


/*!
 * \brief Variables for the disabled diffusive fluxes in the MpNc model within
 *        a finite volume.
 */
// dummy class for the case where diffusion is disabled
template<class TypeTag>
class MPNCVolumeVariablesDiffusion<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief The constructor
     */
    MPNCVolumeVariablesDiffusion()
    {}
    /*!
     * \brief update
     *
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param volVars The volume variables
     * \param problem  The problem
     */
    void update(FluidState &fluidState,
                ParameterCache &paramCache,
                const VolumeVariables &volVars,
                const Problem &problem)
    { }
    /*!
     * \brief The binary diffusion coefficient for each component in the fluid phase.
     * \param compIdx The local index of the components
     */
    Scalar diffCoeffL(const unsigned int compIdx) const
    { return 0; }
    /*!
     * \brief The binary diffusion coefficient for each component in the gas phase.
     * \param compIIdx The local index of the first component in the phase
     * \param compJIdx The local index of the second component in the phase
     */
    Scalar diffCoeffG(const unsigned int compIIdx, const unsigned int compJIdx) const
    { return 0; }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    { }
};

}

#endif // DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH
