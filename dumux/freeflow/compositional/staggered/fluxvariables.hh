// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup FreeflowNCModel
  * \copydoc Dumux::FreeflowNCFluxVariablesImpl
  */
#ifndef DUMUX_FREEFLOW_NC_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_FREEFLOW_NC_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief The flux variables class for the multi-component free-flow model using the staggered grid discretization.
 */

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FreeflowNCFluxVariablesImpl;

template<class TypeTag>
class FreeflowNCFluxVariablesImpl<TypeTag, DiscretizationMethods::Staggered>
: public NavierStokesFluxVariables<TypeTag>
{
    using ParentType = NavierStokesFluxVariables<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    static constexpr auto numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    template<class ElementVolumeVariables, class ElementFaceVariables, class FluxVariablesCache>
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace& scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto diffusiveFluxes = MolecularDiffusionType::flux(problem, element, fvGeometry, elemVolVars, scvf);

        static constexpr auto referenceSystemFormulation = MolecularDiffusionType::referenceSystemFormulation();

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            auto upwindTerm = [compIdx](const auto& volVars)
            {
                const auto density = useMoles ? volVars.molarDensity() : volVars.density();
                const auto fraction =  useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
                return density * fraction;
            };

            flux[compIdx] = ParentType::advectiveFluxForCellCenter(problem, fvGeometry, elemVolVars, elemFaceVars, scvf, upwindTerm);

            // check for the reference system and adapt units of the diffusive flux accordingly.
            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
            {
                flux[compIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx) : diffusiveFluxes[compIdx];
            }
            else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                flux[compIdx] += useMoles ? diffusiveFluxes[compIdx] : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            else
                DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
        }


        // in case one balance is substituted by the total mass balance
        if (ModelTraits::replaceCompEqIdx() < numComponents)
        {
            // accumulate fluxes to a total mass based flux
            Scalar totalMassFlux = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                totalMassFlux += useMoles ? flux[compIdx]*FluidSystem::molarMass(compIdx) : flux[compIdx];

            flux[ModelTraits::replaceCompEqIdx()] = totalMassFlux;
        }

        return flux;
    }
};

} // end namespace

#endif
