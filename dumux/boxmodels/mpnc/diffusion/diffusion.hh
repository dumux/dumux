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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains parts to calculate the diffusive flux in
 *        the fully coupled two-phase N-component model
 */
#ifndef DUMUX_MPNC_DIFFUSION_HH
#define DUMUX_MPNC_DIFFUSION_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/boxmodels/mpnc/mpncproperties.hh>

namespace Dumux {

template <class TypeTag, bool enableDiffusion>
class MPNCDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> ComponentMatrix;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    static void flux(ComponentVector &fluxes,
                     const unsigned int phaseIdx,
                     const FluxVariables &fluxVars,
                     const Scalar molarDensity)
    {
        if (phaseIdx == nPhaseIdx)
            gasFlux_(fluxes, fluxVars, molarDensity);
        else if (phaseIdx == wPhaseIdx){
            #if MACROSCALE_DIFFUSION_ONLY_GAS
                    return ; // in the case that only the diffusion in the gas phase is considered, the liquidFlux should not be called
            #endif
            liquidFlux_(fluxes, fluxVars, molarDensity);
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Invalid phase index: " << phaseIdx);
    }

protected:
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
                fluxVars.face().normal.two_norm() * // because we want a mole flux and not an area specific flux
                fluxVars.porousDiffCoeffL(compIdx) ;
        }
    }

    static void gasFlux_(ComponentVector &fluxes,
                         const FluxVariables &fluxVars,
                         const Scalar molarDensity)
    {
        // Stefan-Maxwell equation
        //
        // See: R. Reid, et al.: "The Properties of Liquids and
        // Gases", 4th edition, 1987, McGraw-Hill, p 596

        // TODO: tensorial diffusion coefficients
        ComponentMatrix M(0);

        for (int compIIdx = 0; compIIdx < numComponents - 1; ++compIIdx) {
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                Scalar Dij = fluxVars.porousDiffCoeffG(compIIdx, compJIdx);
                if (Dij) {
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
