// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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

namespace Dumux {

template <class TypeTag, bool enableDiffusion>
class MPNCDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> DiffMatrix;
    typedef Dune::FieldVector<Scalar, numComponents> DiffVector;
    typedef Dune::FieldVector<Scalar, numComponents> CompVector;

public:
    static void flux(CompVector &fluxes,
                     int phaseIdx,
                     const FluxVariables &fluxDat,
                     Scalar molarDensity)
    {
        if (phaseIdx == gPhaseIdx)
            gasFlux_(fluxes, fluxDat, molarDensity);
        else if (phaseIdx == lPhaseIdx){
            #if MACROSCALE_DIFFUSION_ONLY_GAS
                    return ; // in the case that only the diffusion in the gas phase is considered, the liquidFlux should not be called
            #endif
            liquidFlux_(fluxes, fluxDat, molarDensity);
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Invalid phase index: " << phaseIdx);
    }

protected:
    static void liquidFlux_(CompVector &fluxes,
                            const FluxVariables &fluxDat,
                            Scalar molarDensity)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // TODO: tensorial diffusion coefficients
            Scalar xGrad = fluxDat.moleFracGrad(lPhaseIdx, compIdx)*fluxDat.face().normal;
            fluxes[compIdx] =
                - xGrad *
                molarDensity *
                fluxDat.face().normal.two_norm() * // because we want a mole flux and not an area specific flux
                fluxDat.porousDiffCoeffL(compIdx) ;
        }
    }

    static void gasFlux_(CompVector &fluxes,
                         const FluxVariables &fluxDat,
                         Scalar molarDensity)
    {
        // Stefan-Maxwell equation
        //
        // See: R. Reid, et al.: "The Properties of Liquids and
        // Gases", 4th edition, 1987, McGraw-Hill, p 596

        // TODO: tensorial diffusion coefficients
        DiffMatrix M(0);

        for (int i = 0; i < numComponents - 1; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                Scalar Dij = fluxDat.porousDiffCoeffG(i, j);
                if (Dij) {
                    M[i][j] += fluxDat.moleFrac(gPhaseIdx, i) / Dij;
                    M[i][i] -= fluxDat.moleFrac(gPhaseIdx, j) / Dij;
                }
            }
        };

        for (int i = 0; i < numComponents; ++i) {
            M[numComponents - 1][i] = 1.0;
        }

        DiffVector rightHandSide ; // see source cited above
        for (int i = 0; i < numComponents - 1; ++i) {
            rightHandSide[i] = molarDensity*(fluxDat.moleFracGrad(gPhaseIdx, i)*fluxDat.face().normal);
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
    typedef Dune::FieldVector<Scalar, numComponents>        CompVector;

public:
    static void flux(CompVector &fluxes,
                     int phaseIdx,
                     const FluxVariables &fluxVars,
                     Scalar totalConcentration)
    { fluxes = 0;  }
};

}

#endif // DUMUX_MPNC_DIFFUSION_HH
