// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Maxwell-Stefan's law for control-volume finite element schemes, evaluated at an interpolation point
 */
#ifndef DUMUX_FLUX_CVFE_MAXWELL_STEFANS_LAW__HH
#define DUMUX_FLUX_CVFE_MAXWELL_STEFANS_LAW__HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/cvfe/coefficients_.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Maxwell-Stefan's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * The mass fluxes \f$\mathbf{j}\f$ of the first \f$n-1\f$ components solve
 * \f$ \mathbf{B} \mathbf{j} = -\varrho \nabla \mathbf{x} \f$ with the mole fraction gradients
 * \f$\nabla \mathbf{x}\f$ and the Maxwell-Stefan matrix
 * \f[
 *   B_{ii} = \frac{x_i M}{D_{in} M_n} + \sum_{j \neq i} \frac{x_j M}{D_{ij} M_i}, \quad
 *   B_{ij} = x_i \left( \frac{M}{D_{in} M_n} - \frac{M}{D_{ij} M_j} \right),
 * \f]
 * where \f$M\f$ is the average molar mass of the phase. The flux of the last component closes the
 * mass balance, \f$ \mathbf{j}_n = -\sum_{i<n} \mathbf{j}_i \f$. Density, average molar mass and mole
 * fractions are interpolated with the shape functions.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 * \tparam FluidSystem the fluid system
 * \tparam referenceSystem the reference system of the diffusive fluxes (mass averaged only)
 */
template<class Scalar, class GridDiscretization, class FluidSystem,
         ReferenceSystemFormulation referenceSystem = ReferenceSystemFormulation::massAveraged>
class CVFEMaxwellStefansLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int dimWorld = GridView::dimensionworld;

    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    static_assert(referenceSystem == ReferenceSystemFormulation::massAveraged,
                  "Only the mass averaged reference system is supported for the Maxwell-Stefan formulation");

public:
    //! The flux terms of all components of a phase, one row per component
    using ComponentFluxTerms = Dune::FieldMatrix<Scalar, numComponents, dimWorld>;

    //! The reference system of the diffusive fluxes
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    /*!
     * \brief Returns the diffusive flux terms of all components of a phase at an interpolation point
     * \note Looks up the cache of the interpolation point, see the overload below.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class DiffusionCoefficient>
    static ComponentFluxTerms fluxTerm(const Problem& problem,
                                       const Element& element,
                                       const ElementDiscretization& elemDisc,
                                       const ElementVariables& elemVars,
                                       const int phaseIdx,
                                       const IpData& ipData,
                                       const DiffusionCoefficient& diffusionCoefficient)
    {
        return fluxTerm(problem, element, elemDisc, elemVars, phaseIdx,
                        ipData, cache(elemVars, ipData), diffusionCoefficient);
    }

    /*!
     * \brief Returns the diffusive flux terms of all components of a phase at an interpolation point
     * \note The binary diffusion coefficient of components i and j is given per local dof,
     *       `diffusionCoefficient(vars, i, j)`, and then interpolated, or per interpolation point,
     *       `diffusionCoefficient(elemDisc, elemVars, ipData, i, j)`, for instance
     *       Experimental::CVFE::faceTensorAverageAtScvf. The fluxes vanish where a coefficient
     *       vanishes, for instance in the absence of the phase.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class IpCache, class DiffusionCoefficient>
    static ComponentFluxTerms fluxTerm(const Problem& problem,
                                       const Element& element,
                                       const ElementDiscretization& elemDisc,
                                       const ElementVariables& elemVars,
                                       const int phaseIdx,
                                       const IpData& ipData,
                                       const IpCache& ipCache,
                                       const DiffusionCoefficient& diffusionCoefficient)
    {
        const auto& shapeValues = ipCache.shapeValues();

        Scalar density(0.0);
        Scalar averageMolarMass(0.0);
        ComponentVector moleFraction(0.0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof];
            const auto shapeValue = shapeValues[localDof.index()][0];
            density += vars.density(phaseIdx)*shapeValue;
            averageMolarMass += vars.averageMolarMass(phaseIdx)*shapeValue;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                moleFraction[compIdx] += vars.moleFraction(phaseIdx, compIdx)*shapeValue;
        }

        ComponentFluxTerms fluxTerms(0.0);

        ReducedComponentMatrix msMatrix(0.0);
        const auto lastCompIdx = numComponents-1;
        const auto molarMassLast = FluidSystem::molarMass(lastCompIdx);
        for (int compIIdx = 0; compIIdx < numComponents-1; ++compIIdx)
        {
            const Scalar dIn = Dumux::Detail::CVFE::coefficientAtIp(elemDisc, elemVars, ipData, ipCache, diffusionCoefficient, compIIdx, lastCompIdx);
            if (dIn == 0.0)
                return fluxTerms;

            const auto xi = moleFraction[compIIdx];
            const auto molarMassI = FluidSystem::molarMass(compIIdx);
            msMatrix[compIIdx][compIIdx] += xi*averageMolarMass/(dIn*molarMassLast);

            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                if (compJIdx == compIIdx)
                    continue;

                const Scalar dIj = Dumux::Detail::CVFE::coefficientAtIp(elemDisc, elemVars, ipData, ipCache, diffusionCoefficient, compIIdx, compJIdx);
                if (dIj == 0.0)
                    return fluxTerms;

                const auto xj = moleFraction[compJIdx];
                msMatrix[compIIdx][compIIdx] += xj*averageMolarMass/(dIj*molarMassI);
                if (compJIdx < numComponents-1)
                    msMatrix[compIIdx][compJIdx] += xi*(averageMolarMass/(dIn*molarMassLast)
                                                        - averageMolarMass/(dIj*FluidSystem::molarMass(compJIdx)));
            }
        }

        msMatrix.invert();

        for (int compIdx = 0; compIdx < numComponents-1; ++compIdx)
        {
            GlobalPosition gradX(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                gradX.axpy(elemVars[localDof].moleFraction(phaseIdx, compIdx), ipCache.gradN(localDof.index()));

            for (int compIIdx = 0; compIIdx < numComponents-1; ++compIIdx)
                fluxTerms[compIIdx].axpy(-density*msMatrix[compIIdx][compIdx], gradX);
        }

        for (int compIdx = 0; compIdx < numComponents-1; ++compIdx)
            fluxTerms[lastCompIdx] -= fluxTerms[compIdx];

        return fluxTerms;
    }
};

} // end namespace Dumux::Experimental

#endif
