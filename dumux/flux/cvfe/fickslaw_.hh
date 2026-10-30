// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Fick's law for control-volume finite element schemes, evaluated at an interpolation point
 */
#ifndef DUMUX_FLUX_CVFE_FICKS_LAW__HH
#define DUMUX_FLUX_CVFE_FICKS_LAW__HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/cvfe/coefficients_.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Fick's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * The diffusive flux of component \f$\kappa\f$ in phase \f$\alpha\f$ is
 * \f$ -\varrho_\alpha \mathbf{D}^\kappa_\alpha \nabla x^\kappa_\alpha \f$, with mass density and mass
 * fraction or molar density and mole fraction depending on the reference system. The density and
 * the fraction are interpolated with the shape functions.
 *
 * The flux of the main component of the phase is not computed (unless the fluid system is a tracer
 * fluid system); a model that balances the main component closes it with the negative sum of the
 * other component fluxes.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 * \tparam FluidSystem the fluid system
 * \tparam referenceSystem the reference system of the diffusive fluxes
 */
template<class Scalar, class GridDiscretization, class FluidSystem,
         ReferenceSystemFormulation referenceSystem = ReferenceSystemFormulation::massAveraged>
class CVFEFicksLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int dimWorld = GridView::dimensionworld;

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
     * \note The diffusion coefficient of a component is given per local dof, `diffusionCoefficient(vars, compIdx)`,
     *       and then interpolated, or per interpolation point, `diffusionCoefficient(elemDisc, elemVars, ipData, compIdx)`,
     *       for instance Experimental::CVFE::faceTensorAverageAtScvf. effectiveDiffusionCoefficient() gives the
     *       coefficient of the variables.
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
        for (const auto& localDof : localDofs(elemDisc))
            density += massOrMolarDensity(elemVars[localDof], referenceSystem, phaseIdx)*shapeValues[localDof.index()][0];

        ComponentFluxTerms fluxTerms(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            GlobalPosition gradX(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                gradX.axpy(massOrMoleFraction(elemVars[localDof], referenceSystem, phaseIdx, compIdx),
                           ipCache.gradN(localDof.index()));

            fluxTerms[compIdx] = mv(Dumux::Detail::CVFE::coefficientAtIp(elemDisc, elemVars, ipData, ipCache, diffusionCoefficient, compIdx), gradX);
            fluxTerms[compIdx] *= -density;
        }

        return fluxTerms;
    }

    //! The effective diffusion coefficient of a component in a phase, as given by the variables
    template<class Variables>
    static decltype(auto) effectiveDiffusionCoefficient(const Variables& vars, const int phaseIdx, const int compIdx)
    {
        if constexpr (FluidSystem::isTracerFluidSystem())
            return vars.effectiveDiffusionCoefficient(0, 0, compIdx);
        else
            return vars.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
    }
};

} // end namespace Dumux::Experimental

#endif
