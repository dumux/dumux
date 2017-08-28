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
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 */
#ifndef DUMUX_TRACER_LOCAL_RESIDUAL_HH
#define DUMUX_TRACER_LOCAL_RESIDUAL_HH

namespace Dumux
{

/*!
 * \ingroup Implicit
 * \ingroup TracerLocalResidual
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 *
 */
template<class TypeTag>
class TracerLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param scv The sub control volume
     *  \param volVars The primary and secondary varaibles on the scv
     *  \param useMoles If mole or mass fractions are used
     */
    PrimaryVariables computeStorage(const Problem& problem,
                                    const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        PrimaryVariables storage(0.0);

        // formulation with mole balances
        if (useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                storage[compIdx] += volVars.porosity()
                                    * volVars.molarDensity(0)
                                    * volVars.moleFraction(0, compIdx);
        }
        // formulation with mass balances
        else
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                storage[compIdx] += volVars.porosity()
                                    * volVars.density(0)
                                    * volVars.massFraction(0, compIdx);
        }

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face
     * \param elemFluxVarsCache The cache related to flux compuation
     * \param useMoles If mole or mass fractions are used
     */
    PrimaryVariables computeFlux(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        // get upwind weights into local scope
        PrimaryVariables flux(0.0);
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(0);
        // formulation with mole balances
        if (useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const VolumeVariables& volVars)
                { return volVars.molarDensity()*volVars.moleFraction(0, compIdx); };

                // advective fluxes
                flux[compIdx] += fluxVars.advectiveFlux(0, upwindTerm);
                // diffusive fluxes
                flux[compIdx] += diffusiveFluxes[compIdx];
            }
        }
        // formulation with mass balances
        else
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const VolumeVariables& volVars)
                { return volVars.density()*volVars.massFraction(0, compIdx); };

                // advective fluxes
                flux[compIdx] += fluxVars.advectiveFlux(0, upwindTerm);
                // diffusive fluxes
                flux[compIdx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            }
        }

        return flux;
    }

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars) const
    {
        // we know that these values are constant throughout the simulation
        static const auto phi = curVolVars.porosity();
        static const auto phi_rho = phi*curVolVars.density();

        const auto volume = element.geometry().volume();
        partialDerivatives[0][0] += volume*phi_rho/this->timeLoop().timeStepSize();
    }

    // TODO: IMPLICIT ANALYTICAL DERIVATIVE CONTRIBUTIONS
    // template<class PartialDerivativeMatrix>
    // void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
    //                           const Problem& problem,
    //                           const Element& element,
    //                           const FVElementGeometry& fvGeometry,
    //                           const VolumeVariables& curVolVars) const {}

    // TODO: IMPLICIT ANALYTICAL DERIVATIV CONTRIBUTIONS
    // template<class PartialDerivativeMatrices>
    // void addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
    //                         const Problem& problem,
    //                         const Element& element,
    //                         const FVElementGeometry& fvGeometry,
    //                         const ElementVolumeVariables& curElemVolVars,
    //                         const ElementFluxVariablesCache& elemFluxVarsCache,
    //                         const SubControlVolumeFace& scvf) const
    // {}

    // TODO: IMPLICIT ANALYTICAL DERIVATIV CONTRIBUTIONS
    // template<class PartialDerivativeMatrices>
    // void addDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
    //                                  const Problem& problem,
    //                                  const Element& element,
    //                                  const FVElementGeometry& fvGeometry,
    //                                  const ElementVolumeVariables& curElemVolVars,
    //                                  const ElementFluxVariablesCache& elemFluxVarsCache,
    //                                  const SubControlVolumeFace& scvf) const
    // {}

    // TODO: IMPLICIT ANALYTICAL DERIVATIV CONTRIBUTIONS
    // template<class PartialDerivativeMatrices>
    // void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
    //                              const Problem& problem,
    //                              const Element& element,
    //                              const FVElementGeometry& fvGeometry,
    //                              const ElementVolumeVariables& curElemVolVars,
    //                              const ElementFluxVariablesCache& elemFluxVarsCache,
    //                              const SubControlVolumeFace& scvf) const
    // {}
};

} // end namespace Dumux

#endif
