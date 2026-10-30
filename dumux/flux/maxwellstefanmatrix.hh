// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief The matrix of the Maxwell-Stefan equations
 */
#ifndef DUMUX_FLUX_MAXWELL_STEFAN_MATRIX_HH
#define DUMUX_FLUX_MAXWELL_STEFAN_MATRIX_HH

#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Flux
 * \brief The matrix \f$\mathbf{B}\f$ of the Maxwell-Stefan equations
 *        \f$ \mathbf{B} \mathbf{j} = -\varrho \nabla \mathbf{x} \f$ of a phase in the mass-averaged
 *        reference system, for the mass fluxes \f$\mathbf{j}\f$ of the first \f$n-1\f$ components
 *        (see R. Taylor and R. Krishna, Multicomponent Mass Transfer, 1993)
 *
 * \f[
 *   B_{ii} = \frac{x_i M}{D_{in} M_n} + \sum_{j \neq i} \frac{x_j M}{D_{ij} M_i}, \quad
 *   B_{ij} = x_i \left( \frac{M}{D_{in} M_n} - \frac{M}{D_{ij} M_j} \right),
 * \f]
 * where \f$M\f$ is the average molar mass of the phase and \f$M_i\f$ the molar mass of component \f$i\f$.
 *
 * \param moleFraction the mole fractions of all components
 * \param averageMolarMass the average molar mass of the phase
 * \param diffusionCoefficient the binary diffusion coefficient of two components, `diffusionCoefficient(compIIdx, compJIdx)`
 */
template<class FluidSystem, class Scalar, int numComponents, class DiffusionCoefficient>
Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>
maxwellStefanMatrix(const Dune::FieldVector<Scalar, numComponents>& moleFraction,
                    const std::type_identity_t<Scalar> averageMolarMass,
                    const DiffusionCoefficient& diffusionCoefficient)
{
    Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1> msMatrix(0.0);

    const auto lastCompIdx = numComponents-1;
    const auto molarMassLast = FluidSystem::molarMass(lastCompIdx);
    for (int compIIdx = 0; compIIdx < numComponents-1; ++compIIdx)
    {
        const auto xi = moleFraction[compIIdx];
        const auto molarMassI = FluidSystem::molarMass(compIIdx);
        const auto dIn = diffusionCoefficient(compIIdx, lastCompIdx);
        msMatrix[compIIdx][compIIdx] += xi*averageMolarMass/(dIn*molarMassLast);

        for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            if (compJIdx == compIIdx)
                continue;

            const auto xj = moleFraction[compJIdx];
            const auto dIj = diffusionCoefficient(compIIdx, compJIdx);
            msMatrix[compIIdx][compIIdx] += xj*averageMolarMass/(dIj*molarMassI);
            if (compJIdx < numComponents-1)
                msMatrix[compIIdx][compJIdx] += xi*(averageMolarMass/(dIn*molarMassLast)
                                                    - averageMolarMass/(dIj*FluidSystem::molarMass(compJIdx)));
        }
    }

    return msMatrix;
}

} // end namespace Dumux::Detail

#endif
