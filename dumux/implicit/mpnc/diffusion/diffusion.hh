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
 * \brief This file contains parts to calculate the diffusive flux in
 *        the fully coupled MpNc model
 */
#ifndef DUMUX_MPNC_DIFFUSION_HH
#define DUMUX_MPNC_DIFFUSION_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>


#include <dumux/implicit/mpnc/mpncproperties.hh>

namespace Dumux {
/*!
 * \brief Calculates the diffusive flux in the fully coupled MpNc model.
 *        Called from the mass module.
 */
template <class TypeTag, bool enableDiffusion>
class MPNCDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nCompIdx = FluidSystem::nCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> ComponentMatrix;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:

    /*!
     * \brief Calls the function for the diffusion in the gas and liquid phases, respectively.
     *
     *        In the gas phase, the mutual influence of mole fractions can be considered.
     *
     * \param fluxes The Diffusive flux over the sub-control-volume face for each component
     * \param phaseIdx The index of the phase we are calculating the diffusive flux for
     * \param fluxVars The flux variables at the current sub control volume face
     * \param molarDensity The (molar) density of the phase
     */
    static void flux(ComponentVector &fluxes,
                     const unsigned int phaseIdx,
                     const FluxVariables &fluxVars,
                     const Scalar molarDensity )
    {
        if ( not FluidSystem::isLiquid(phaseIdx) )
            gasFlux_(fluxes, fluxVars, molarDensity);
        else if ( FluidSystem::isLiquid(phaseIdx) ){
            #if MACROSCALE_DIFFUSION_ONLY_GAS
                    return ; // in the case that only the diffusion in the gas phase is considered,
                             // the liquidFlux should not be called
            #endif
            liquidFlux_(fluxes, fluxVars, molarDensity);
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Invalid phase index: " << phaseIdx);
    }

protected:

    /*!
     * \brief Calculates the diffusive flux in the liquid phase: only Fick diffusion (no mutual influence) is considered.
     *
     * \param fluxes The Diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub control volume face
     * \param molarDensity The (molar) density of the phase
     */
    static void liquidFlux_(ComponentVector &fluxes,
                            const FluxVariables &fluxVars,
                            const Scalar molarDensity)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // TODO: tensorial diffusion coefficients
            const Scalar xGrad = fluxVars.moleFractionGrad(wPhaseIdx, compIdx)*fluxVars.face().normal;
            fluxes[compIdx] =
                - xGrad *
                molarDensity *
                fluxVars.porousDiffCoeffL(compIdx) ;
        }
    }

    /*!
     * \brief Calculates the diffusive flux in the gas phase:
     *        The property UseMaxwellDiffusion selects whether Maxwell or Fick diffusion is applied.
     *        Has to be the same for a 2-component system. However, Maxwells does not work always.
     *
     * \param fluxes The Diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub control volume face
     * \param molarDensity The (molar) density of the phase
     *
     * Reid et al. (1987, p. 596) \cite reid1987 <BR>
     */
    static void gasFlux_(ComponentVector &fluxes,
                         const FluxVariables &fluxVars,
                         const Scalar molarDensity)
    {
        // Alternative: use Fick Diffusion: no mutual influence of diffusion
        if (GET_PROP_VALUE(TypeTag, UseMaxwellDiffusion) ){
            // Stefan-Maxwell equation
            //
            // See: R. Reid, et al.: "The Properties of Liquids and
            // Gases", 4th edition, 1987, McGraw-Hill, p 596

            // TODO: tensorial diffusion coefficients
            ComponentMatrix M(0);

            for (int compIIdx = 0; compIIdx < numComponents - 1; ++compIIdx) {
                for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                    Scalar Dij = fluxVars.porousDiffCoeffG(compIIdx, compJIdx);
                    if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(Dij, 0.0, 1.0e-30)) {
                        M[compIIdx][compJIdx] += fluxVars.moleFraction(nPhaseIdx, compIIdx) / Dij;
                        M[compIIdx][compIIdx] -= fluxVars.moleFraction(nPhaseIdx, compJIdx) / Dij;
                    }
                }
            }

            for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
                M[numComponents - 1][compIIdx] = 1.0;
            }

            ComponentVector rightHandSide ; // see source cited above
            for (int compIIdx = 0; compIIdx < numComponents - 1; ++compIIdx) {
                rightHandSide[compIIdx] = molarDensity*(fluxVars.moleFractionGrad(nPhaseIdx, compIIdx)*fluxVars.face().normal);
            }
            rightHandSide[numComponents - 1] = 0.0;

            M.solve(fluxes, rightHandSide);
        }
        else{// Fick Diffusion
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                // TODO: tensorial diffusion coefficients
                const Scalar xGrad = fluxVars.moleFractionGrad(nPhaseIdx, compIdx)*fluxVars.face().normal;
                fluxes[compIdx] =
                    - xGrad *
                    molarDensity
                    * fluxVars.porousDiffCoeffG(compIdx, nCompIdx) ; // this is == 0 for nComp==comp,
                                                                     // i.e. no diffusion of the main component of the phase
                }
        }
    }

    // return whether a concentration can be assumed to be a trace
    // component in the context of diffusion
    static Scalar isTraceComp_(Scalar x)
    { return x < 0.5/numComponents; }
};

/*!
 * \brief Specialization of the diffusion module for the case where
 *        diffusion is disabled.
 *
 * This class just does nothing.
 */
template <class TypeTag>
class MPNCDiffusion<TypeTag, false>
{
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef Dune::FieldVector<Scalar, numComponents>        ComponentVector;

public:
    static void flux(ComponentVector &fluxes,
                     const unsigned int phaseIdx,
                     const FluxVariables &fluxVars,
                     const Scalar molarDensity)
    { fluxes = 0;  }
};

}

#endif // DUMUX_MPNC_DIFFUSION_HH
