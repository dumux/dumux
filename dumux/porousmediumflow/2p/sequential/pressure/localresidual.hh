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
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 */
#ifndef DUMUX_PRESSURE_LOCAL_RESIDUAL_HH
#define DUMUX_PRESSURE_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 *
 */
template<class TypeTag>
class PressureLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    static constexpr int numPhases = ModelTraits::numPhases();
    static constexpr int pressureEqIdx = ModelTraits::Indices::pressureEqIdx; //!< first index for the mass balance
    static constexpr int pressureIdx = ModelTraits::Indices::pressureIdx; //!< first index for the mass balance


public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem TODO docme!
     * \param scv The sub control volume
     * \param volVars The primary and secondary varaibles on the scv
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        // Todo what about a slightly compressible case?

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem TODO docme!
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
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // the physical quantities for which we perform upwinding
            auto upwindTerm = [phaseIdx](const auto& volVars)
                              { return volVars.mobility(phaseIdx); };

            flux[pressureEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
        }

        return flux;
    }

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const {}

    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    {
        //problem.addSourceDerivatives(partialDerivatives, element, fvGeometry, curVolVars, scv);
    }

    //! flux derivatives for the cell-centered tpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GET_PROP_TYPE(T, FVGridGeometry)::discMethod == DiscretizationMethod::cctpfa, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "2p/sequential/pressure/localresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(1),
                      "2p/sequential/pressure/localresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "2p/sequential/pressure/localresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(1),
                      "2p/sequential/pressure/localresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(ModelTraits::numPhases() == 2,
                      "2p/sequential/pressure/localresidual.hh: Only two-phase models are allowed!");
        //static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1,
        //              "2p/sequential/pressure/incompressiblelocalresidual.hh: Analytic differentiation has to be checked for p1-s0 formulation!");

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 1, elemFluxVarsCache);

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto outsideElement = fvGeometry.fvGridGeometry().element(outsideScvIdx);
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[outsideScvIdx];
        const auto& insideMaterialParams = problem.spatialParams().materialLawParams(element,
                                                                                     insideScv,
                                                                                     elementSolution<FVElementGeometry>(insideVolVars.priVars()));
        const auto& outsideMaterialParams = problem.spatialParams().materialLawParams(outsideElement,
                                                                                      outsideScv,
                                                                                      elementSolution<FVElementGeometry>(outsideVolVars.priVars()));

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        const auto mobW_inside = insideVolVars.mobility(0);
        const auto mobN_inside = insideVolVars.mobility(1);
        const auto mobW_outside = outsideVolVars.mobility(0);
        const auto mobN_outside = outsideVolVars.mobility(1);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // precalculate values
        const auto mobup_w = (flux_w > 0.) ? mobW_inside : mobW_outside;
        const auto mobup_n = (flux_n > 0.) ? mobN_inside : mobN_outside;
        const auto tij_up_t = tij*(mobup_w + mobup_n);

        // partial derivative of the wetting phase flux w.r.t. p_w
        // add partial derivatives to the respective given matrices
        derivativeMatrices[scvf.insideScvIdx()][pressureEqIdx][pressureIdx] += tij_up_t;
        derivativeMatrices[scvf.outsideScvIdx()][pressureEqIdx][pressureIdx] -= tij_up_t;
    }

    //! Dirichlet flux derivatives for the cell-centered tpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GET_PROP_TYPE(T, FVGridGeometry)::discMethod == DiscretizationMethod::cctpfa, void>
    addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const SubControlVolumeFace& scvf) const
    {
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 1, elemFluxVarsCache);

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];
        const auto& insideMaterialParams = problem.spatialParams().materialLawParams(element,
                                                                                     insideScv,
                                                                                     elementSolution<FVElementGeometry>(insideVolVars.priVars()));

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        const auto mobW_inside = insideVolVars.mobility(0);
        const auto mobN_inside = insideVolVars.mobility(1);
        const auto mobW_outside = outsideVolVars.mobility(0);
        const auto mobN_outside = outsideVolVars.mobility(1);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // precalculate values
        const auto mobup_w = (flux_w > 0.) ? mobW_inside : mobW_outside;
        const auto mobup_n = (flux_n > 0.) ? mobN_inside : mobN_outside;
        const auto tij_up_t = tij*(mobup_w + mobup_n);

        // partial derivative of the wetting phase flux w.r.t. p_w
        // add partial derivatives to the respective given matrices
        derivativeMatrices[scvf.insideScvIdx()][pressureEqIdx][pressureIdx] += tij_up_t;
    }

    //! Robin-type flux derivatives
    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        //! Robin-type boundary conditions are problem-specific.
        //! We can't put a general implementation here - users defining Robin-type BCs
        //! while using analytical Jacobian assembly must overload this function!
    }

};

} // end namespace Dumux

#endif
