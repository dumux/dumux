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
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 */

#ifndef DUMUX_TRACER_LOCAL_RESIDUAL_HH
#define DUMUX_TRACER_LOCAL_RESIDUAL_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

/*!
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 */
template<class TypeTag>
class TracerLocalResidual: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int phaseIdx = 0;

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The primary and secondary varaibles on the scv
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        // regularize saturation so we don't get singular matrices when the saturation is zero
        // note that the fluxes will still be zero (zero effective diffusion coefficient),
        // and we still solve the equation storage = 0 yielding the correct result
        using std::max;
        const Scalar saturation = max(1e-8, volVars.saturation(phaseIdx));

        // formulation with mole balances
        if (useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                storage[compIdx] += volVars.porosity()
                                    * volVars.molarDensity(phaseIdx)
                                    * volVars.moleFraction(phaseIdx, compIdx)
                                    * saturation;
        }
        // formulation with mass balances
        else
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                storage[compIdx] += volVars.porosity()
                                    * volVars.density(phaseIdx)
                                    * volVars.massFraction(phaseIdx, compIdx)
                                    * saturation;
        }

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face
     * \param elemFluxVarsCache The cache related to flux compuation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        NumEqVector flux(0.0);
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);

        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        // formulation with mole balances
        if (useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const auto& volVars)
                { return volVars.molarDensity()*volVars.moleFraction(phaseIdx, compIdx); };

                // advective fluxes
                flux[compIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                // diffusive fluxes
                if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    flux[compIdx] += diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx);
                else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    flux[compIdx] += diffusiveFluxes[compIdx];
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            }
        }
        // formulation with mass balances
        else
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const auto& volVars)
                { return volVars.density()*volVars.massFraction(phaseIdx, compIdx); };

                // advective fluxes
                flux[compIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                // diffusive fluxes
                if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    flux[compIdx] += diffusiveFluxes[compIdx];
                else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    flux[compIdx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            }
        }

        return flux;
    }

    /*!
     * \brief TODO docme!
     *
     * \param partialDerivatives TODO docme!
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param curVolVars The current volume variables
     * \param scv The sub control volume
     */
    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const
    {
        // regularize saturation so we don't get singular matrices when the saturation is zero
        // note that the fluxes will still be zero (zero effective diffusion coefficient),
        // and we still solve the equation storage = 0 yielding the correct result
        using std::max;
        const auto saturation = max(1e-8, curVolVars.saturation(phaseIdx));

        const auto porosity = curVolVars.porosity();
        const auto rho = useMoles ? curVolVars.molarDensity() : curVolVars.density();
        const auto d_storage = Extrusion::volume(scv)*porosity*rho*saturation/this->timeLoop().timeStepSize();

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            partialDerivatives[compIdx][compIdx] += d_storage;
    }

    /*!
     * \brief TODO docme!
     *
     * \param partialDerivatives TODO docme!
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param curVolVars The current volume variables
     * \param scv The sub control volume
     */
    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    {
        // TODO maybe forward to the problem? -> necessary for reaction terms
    }

    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod != DiscretizationMethod::box, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        if constexpr (FVElementGeometry::GridGeometry::discMethod != DiscretizationMethod::cctpfa)
            DUNE_THROW(Dune::NotImplemented, "Analytic flux differentiation only implemented for tpfa");

        // advective term: we do the same for all tracer components
        auto rho = [](const VolumeVariables& volVars)
        { return useMoles ? volVars.molarDensity() : volVars.density(); };

        // the volume flux
        const auto volFlux = problem.spatialParams().volumeFlux(element, fvGeometry, curElemVolVars, scvf);

        // the upwind weight
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

        // get the inside and outside volvars
        const auto& insideVolVars = curElemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];

        const Scalar insideWeight = std::signbit(volFlux) ? (1.0 - upwindWeight) : upwindWeight;
        const Scalar outsideWeight = 1.0 - insideWeight;
        const auto advDerivII = volFlux*rho(insideVolVars)*insideWeight;
        const auto advDerivIJ = volFlux*rho(outsideVolVars)*outsideWeight;

        // diffusive term
        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        const auto& fluxCache = elemFluxVarsCache[scvf];
        const Scalar rhoInside = massOrMolarDensity(insideVolVars, referenceSystemFormulation, phaseIdx);
        const Scalar rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystemFormulation, phaseIdx);
        const Scalar massOrMolarDensity = 0.5*(rhoInside + rhoOutside);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // diffusive term
            Scalar diffDeriv = 0.0;
            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
            {
                diffDeriv = useMoles ? massOrMolarDensity*fluxCache.diffusionTij(phaseIdx, compIdx)/FluidSystem::molarMass(compIdx)
                                            : massOrMolarDensity*fluxCache.diffusionTij(phaseIdx, compIdx);
            }
            else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
            {
                diffDeriv = useMoles ? massOrMolarDensity*fluxCache.diffusionTij(phaseIdx, compIdx)
                                            : massOrMolarDensity*fluxCache.diffusionTij(phaseIdx, compIdx)*FluidSystem::molarMass(compIdx);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");

            derivativeMatrices[scvf.insideScvIdx()][compIdx][compIdx] += (advDerivII + diffDeriv);
            if (!scvf.boundary())
                derivativeMatrices[scvf.outsideScvIdx()][compIdx][compIdx] += (advDerivIJ - diffDeriv);
        }
    }

    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::box, void>
    addFluxDerivatives(JacobianMatrix& A,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {

        // advective term: we do the same for all tracer components
        auto rho = [](const VolumeVariables& volVars)
        { return useMoles ? volVars.molarDensity() : volVars.density(); };

        // the volume flux
        const auto volFlux = problem.spatialParams().volumeFlux(element, fvGeometry, curElemVolVars, scvf);

        // the upwind weight
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

        // get the inside and outside volvars
        const auto& insideVolVars = curElemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];

        const auto insideWeight = std::signbit(volFlux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight = 1.0 - insideWeight;
        const auto advDerivII = volFlux*rho(insideVolVars)*insideWeight;
        const auto advDerivIJ = volFlux*rho(outsideVolVars)*outsideWeight;

        // diffusive term
        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        using DiffusionType = GetPropType<T, Properties::MolecularDiffusionType>;
        const auto ti = DiffusionType::calculateTransmissibilities(problem,
                                                                   element,
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   scvf,
                                                                   elemFluxVarsCache[scvf],
                                                                   phaseIdx);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            for (const auto& scv : scvs(fvGeometry))
            {
                // diffusive term
                auto diffDeriv = 0.0;
                if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    diffDeriv += useMoles ? ti[compIdx][scv.indexInElement()]/FluidSystem::molarMass(compIdx)
                                        : ti[compIdx][scv.indexInElement()];
                else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    diffDeriv += useMoles ? ti[compIdx][scv.indexInElement()]
                                            : ti[compIdx][scv.indexInElement()]*FluidSystem::molarMass(compIdx);
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
                A[insideScv.dofIndex()][scv.dofIndex()][compIdx][compIdx] += diffDeriv;
                A[outsideScv.dofIndex()][scv.dofIndex()][compIdx][compIdx] -= diffDeriv;
            }

            A[insideScv.dofIndex()][insideScv.dofIndex()][compIdx][compIdx] += advDerivII;
            A[insideScv.dofIndex()][outsideScv.dofIndex()][compIdx][compIdx] += advDerivIJ;
            A[outsideScv.dofIndex()][outsideScv.dofIndex()][compIdx][compIdx] -= advDerivII;
            A[outsideScv.dofIndex()][insideScv.dofIndex()][compIdx][compIdx] -= advDerivIJ;
        }
    }

    template<class PartialDerivativeMatrices>
    void addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                       const Problem& problem,
                                       const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& curElemVolVars,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const SubControlVolumeFace& scvf) const
    {
        // do the same as for inner facets
        addFluxDerivatives(derivativeMatrices, problem, element, fvGeometry,
                           curElemVolVars, elemFluxVarsCache, scvf);
    }

    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        // TODO maybe forward to the problem?
    }
};

} // end namespace Dumux

#endif
