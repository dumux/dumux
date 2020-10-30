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
 * \ingroup TwoPModel
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a two-phase, incompressible test problem.
 */

#ifndef DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH
#define DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH

#include <cmath>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a two-phase, incompressible test problem.
 */
template<class TypeTag>
class TwoPIncompressibleLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    using ParentType = ImmiscibleLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int pressureIdx = ModelTraits::Indices::pressureIdx;
    static constexpr int saturationIdx = ModelTraits::Indices::saturationIdx;
    static constexpr int conti0EqIdx = ModelTraits::Indices::conti0EqIdx;

public:
    using ParentType::ParentType;

    /*!
     * \brief Adds storage derivatives for wetting and nonwetting phase
     *
     * Compute storage derivatives for the wetting and the nonwetting phase with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param partialDerivatives The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
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
        static_assert(!FluidSystem::isCompressible(0),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(1),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(ModelTraits::numFluidPhases() == 2,
                      "2p/incompressiblelocalresidual.hh: Only two-phase models are allowed!");
        static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1,
                      "2p/incompressiblelocalresidual.hh: Analytic differentiation has to be checked for p1-s0 formulation!");

        const auto poreVolume = Extrusion::volume(scv)*curVolVars.porosity();
        const auto dt = this->timeLoop().timeStepSize();

        // partial derivative of the phase storage terms w.r.t. S_n
        partialDerivatives[conti0EqIdx+0][saturationIdx] -= poreVolume*curVolVars.density(0)/dt;
        partialDerivatives[conti0EqIdx+1][saturationIdx] += poreVolume*curVolVars.density(1)/dt;
    }

    /*!
     * \brief Adds source derivatives for wetting and nonwetting phase.
     *
     * \param partialDerivatives The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
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
    { /* TODO maybe forward to problem for the user to implement the source derivatives?*/ }

    /*!
     * \brief Adds flux derivatives for wetting and nonwetting phase for cell-centered FVM using TPFA
     *
     * Compute derivatives for the wetting and the nonwetting phase flux with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param derivativeMatrices The partial derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::cctpfa, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(1),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(1),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(ModelTraits::numFluidPhases() == 2,
                      "2p/incompressiblelocalresidual.hh: Only two-phase models are allowed!");
        static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1,
                      "2p/incompressiblelocalresidual.hh: Analytic differentiation has to be checked for p1-s0 formulation!");

        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 1, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto outsideElement = fvGeometry.gridGeometry().element(outsideScvIdx);
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[outsideScvIdx];

        // old material law interface is deprecated: Replace this by
        // const auto& insidefluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element,
        //                                                                                           insideScv,
        //                                                                                           elementSolution<FVElementGeometry>(insideVolVars.priVars()));
        // const auto& outsidefluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(outsideElement,
        //                                                                                            outsideScv,
        //                                                                                            elementSolution<FVElementGeometry>(outsideVolVars.priVars()));
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto& insidefluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                          element,
                                                                          insideScv,
                                                                          elementSolution<FVElementGeometry>(insideVolVars.priVars()));
        const auto& outsidefluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                           outsideElement,
                                                                           outsideScv,
                                                                           elementSolution<FVElementGeometry>(outsideVolVars.priVars()));

        // get references to the two participating derivative matrices
        auto& dI_dI = derivativeMatrices[insideScvIdx];
        auto& dI_dJ = derivativeMatrices[outsideScvIdx];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(0);
        static const auto rho_n = insideVolVars.density(1);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(0);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(1);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(0);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(1);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(0);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(1);

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(0);
        const auto outsideSw = outsideVolVars.saturation(0);
        const auto dKrw_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrw_dsw(insideSw);
        const auto dKrw_dSn_outside = -1.0*outsidefluidMatrixInteraction.dkrw_dsw(outsideSw);
        const auto dKrn_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrn_dsw(insideSw);
        const auto dKrn_dSn_outside = -1.0*outsidefluidMatrixInteraction.dkrn_dsw(outsideSw);
        const auto dpc_dSn_inside = -1.0*insidefluidMatrixInteraction.dpc_dsw(insideSw);
        const auto dpc_dSn_outside = -1.0*outsidefluidMatrixInteraction.dpc_dsw(outsideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // precalculate values
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        const auto rho_mu_flux_w = rhow_muw*flux_w;
        const auto rho_mu_flux_n = rhon_mun*flux_n;
        const auto tij_up_w = tij*up_w;
        const auto tij_up_n = tij*up_n;

        // partial derivative of the wetting phase flux w.r.t. p_w
        dI_dI[conti0EqIdx+0][pressureIdx] += tij_up_w;
        dI_dJ[conti0EqIdx+0][pressureIdx] -= tij_up_w;

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dI[conti0EqIdx+0][saturationIdx] += rho_mu_flux_w*dKrw_dSn_inside*insideWeight_w;
        dI_dJ[conti0EqIdx+0][saturationIdx] += rho_mu_flux_w*dKrw_dSn_outside*outsideWeight_w;

        // partial derivative of the nonwetting phase flux w.r.t. p_w
        dI_dI[conti0EqIdx+1][pressureIdx] += tij_up_n;
        dI_dJ[conti0EqIdx+1][pressureIdx] -= tij_up_n;

        // partial derivative of the nonwetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dI[conti0EqIdx+1][saturationIdx] += rho_mu_flux_n*dKrn_dSn_inside*insideWeight_n;
        dI_dJ[conti0EqIdx+1][saturationIdx] += rho_mu_flux_n*dKrn_dSn_outside*outsideWeight_n;

        // partial derivative of the nonwetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dI[conti0EqIdx+1][saturationIdx] += tij_up_n*dpc_dSn_inside;
        dI_dJ[conti0EqIdx+1][saturationIdx] -= tij_up_n*dpc_dSn_outside;
    }

    /*!
     * \brief Adds flux derivatives for wetting and nonwetting phase for box method
     *
     * Compute derivatives for the wetting and the nonwetting phase flux with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param A The Jacobian Matrix
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
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
        static_assert(!FluidSystem::isCompressible(0),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(1),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(1),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(ModelTraits::numFluidPhases() == 2,
                      "2p/incompressiblelocalresidual.hh: Only two-phase models are allowed!");
        static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1,
                      "2p/incompressiblelocalresidual.hh: Analytic differentiation has to be checked for p0-s1 formulation!");

        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 1, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScv];
        const auto& outsideVolVars = curElemVolVars[outsideScv];

        const auto elemSol = elementSolution(element, curElemVolVars, fvGeometry);

        // old material law interface is deprecated: Replace this by
        // const auto& insidefluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element,
        //                                                                                           insideScv,
        //                                                                                           elemSol);
        // const auto& outsidefluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element,
        //                                                                                            outsideScv,
        //                                                                                            elemSol);
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto& insidefluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                          element,
                                                                          insideScv,
                                                                          elemSol);
        const auto& outsidefluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                           element,
                                                                           outsideScv,
                                                                           elemSol);

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(0);
        static const auto rho_n = insideVolVars.density(1);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(0);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(1);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(0);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(1);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(0);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(1);

        // let the Law for the advective fluxes calculate the transmissibilities
        const auto ti = AdvectionType::calculateTransmissibilities(problem,
                                                                   element,
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   scvf,
                                                                   elemFluxVarsCache[scvf]);

        // get the rows of the jacobian matrix for the inside/outside scv
        auto& dI_dJ_inside = A[insideScv.dofIndex()];
        auto& dI_dJ_outside = A[outsideScv.dofIndex()];

        // precalculate values
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        const auto rho_mu_flux_w = rhow_muw*flux_w;
        const auto rho_mu_flux_n = rhon_mun*flux_n;

        // add the partial derivatives w.r.t all scvs in the element
        for (const auto& scvJ : scvs(fvGeometry))
        {
            const auto globalJ = scvJ.dofIndex();
            const auto localJ = scvJ.indexInElement();

            // the transmissibily associated with the scvJ
            const auto tj = ti[localJ];

            // partial derivative of the wetting phase flux w.r.t. p_w
            const auto tj_up_w = tj*up_w;
            dI_dJ_inside[globalJ][conti0EqIdx+0][pressureIdx] += tj_up_w;
            dI_dJ_outside[globalJ][conti0EqIdx+0][pressureIdx] -= tj_up_w;

            // partial derivative of the nonwetting phase flux w.r.t. p_w
            const auto tj_up_n = tj*up_n;
            dI_dJ_inside[globalJ][conti0EqIdx+1][pressureIdx] += tj_up_n;
            dI_dJ_outside[globalJ][conti0EqIdx+1][pressureIdx] -= tj_up_n;

            // partial derivatives w.r.t. S_n (are the negative of those w.r.t sw)
            // relative permeability contributions only for inside/outside
            if (localJ == insideScvIdx)
            {
                // partial derivative of the wetting phase flux w.r.t. S_n
                const auto insideSw = insideVolVars.saturation(0);
                const auto dKrw_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrw_dsw(insideSw);
                const auto dFluxW_dSnJ = rho_mu_flux_w*dKrw_dSn_inside*insideWeight_w;
                dI_dJ_inside[globalJ][conti0EqIdx+0][saturationIdx] += dFluxW_dSnJ;
                dI_dJ_outside[globalJ][conti0EqIdx+0][saturationIdx] -= dFluxW_dSnJ;

                // partial derivative of the nonwetting phase flux w.r.t. S_n (k_rn contribution)
                const auto dKrn_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrn_dsw(insideSw);
                const auto dFluxN_dSnJ_krn = rho_mu_flux_n*dKrn_dSn_inside*insideWeight_n;
                dI_dJ_inside[globalJ][conti0EqIdx+1][saturationIdx] += dFluxN_dSnJ_krn;
                dI_dJ_outside[globalJ][conti0EqIdx+1][saturationIdx] -= dFluxN_dSnJ_krn;

                // partial derivative of the nonwetting phase flux w.r.t. S_n (p_c contribution)
                const auto dFluxN_dSnJ_pc = -1.0*tj_up_n*insidefluidMatrixInteraction.dpc_dsw(insideSw);
                dI_dJ_inside[globalJ][conti0EqIdx+1][saturationIdx] += dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][conti0EqIdx+1][saturationIdx] -= dFluxN_dSnJ_pc;
            }
            else if (localJ == outsideScvIdx)
            {
                // see comments for (globalJ == insideScvIdx)
                const auto outsideSw = outsideVolVars.saturation(0);
                const auto dKrw_dSn_outside = -1.0*outsidefluidMatrixInteraction.dkrw_dsw(outsideSw);
                const auto dFluxW_dSnJ = rho_mu_flux_w*dKrw_dSn_outside*outsideWeight_w;
                dI_dJ_inside[globalJ][conti0EqIdx+0][saturationIdx] += dFluxW_dSnJ;
                dI_dJ_outside[globalJ][conti0EqIdx+0][saturationIdx] -= dFluxW_dSnJ;

                const auto dKrn_dSn_outside = -1.0*outsidefluidMatrixInteraction.dkrn_dsw(outsideSw);
                const auto dFluxN_dSnJ_krn = rho_mu_flux_n*dKrn_dSn_outside*outsideWeight_n;
                dI_dJ_inside[globalJ][conti0EqIdx+1][saturationIdx] += dFluxN_dSnJ_krn;
                dI_dJ_outside[globalJ][conti0EqIdx+1][saturationIdx] -= dFluxN_dSnJ_krn;

                const auto dFluxN_dSnJ_pc = -1.0*tj_up_n*outsidefluidMatrixInteraction.dpc_dsw(outsideSw);
                dI_dJ_inside[globalJ][conti0EqIdx+1][saturationIdx] += dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][conti0EqIdx+1][saturationIdx] -= dFluxN_dSnJ_pc;
            }
            else
            {
                // old material law interface is deprecated: Replace this by
                // const auto& fluidMatrixInteractionJ = problem.spatialParams().fluidMatrixInteraction(element,
                //                                                                                      scvJ,
                //                                                                                      elemSol);

                // after the release of 3.3, when the deprecated interface is no longer supported
                const auto& fluidMatrixInteractionJ = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                             element,
                                                                             scvJ,
                                                                             elemSol);
                const auto swJ = curElemVolVars[scvJ].saturation(0);
                const auto dFluxN_dSnJ_pc = tj_up_n*fluidMatrixInteractionJ.dpc_dsw(swJ);
                dI_dJ_inside[globalJ][conti0EqIdx+1][saturationIdx] -= dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][conti0EqIdx+1][saturationIdx] += dFluxN_dSnJ_pc;
            }
        }
    }

    /*!
     * \brief Adds cell-centered Dirichlet flux derivatives for wetting and nonwetting phase
     *
     * Compute derivatives for the wetting and the nonwetting phase flux with respect to \f$p_w\f$
     * and \f$S_n\f$.
     *
     * \param derivativeMatrices The matrices containing the derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices>
    void addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                       const Problem& problem,
                                       const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& curElemVolVars,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const SubControlVolumeFace& scvf) const
    {
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 0, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars, scvf, 1, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];

        // old material law interface is deprecated: Replace this by
        // const auto& insidefluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element,
        //                                                                                           insideScv,
        //                                                                                           elementSolution<FVElementGeometry>(insideVolVars.priVars()));
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto& insidefluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem.spatialParams(),
                                                                          element,
                                                                          insideScv,
                                                                          elementSolution<FVElementGeometry>(insideVolVars.priVars()));

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(0);
        static const auto rho_n = insideVolVars.density(1);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(0);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(1);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(0);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(1);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(0);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(1);

        // get reference to the inside derivative matrix
        auto& dI_dI = derivativeMatrices[insideScvIdx];

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(0);
        const auto dKrw_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrw_dsw(insideSw);
        const auto dKrn_dSn_inside = -1.0*insidefluidMatrixInteraction.dkrn_dsw(insideSw);
        const auto dpc_dSn_inside = -1.0*insidefluidMatrixInteraction.dpc_dsw(insideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();
        // partial derivative of the wetting phase flux w.r.t. p_w
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        dI_dI[conti0EqIdx+0][pressureIdx] += tij*up_w;

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dI[conti0EqIdx+0][saturationIdx] += rhow_muw*flux_w*dKrw_dSn_inside*insideWeight_w;

        // partial derivative of the nonwetting phase flux w.r.t. p_w
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        dI_dI[conti0EqIdx+1][pressureIdx] += tij*up_n;

        // partial derivative of the nonwetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dI[conti0EqIdx+1][saturationIdx] += rhon_mun*flux_n*dKrn_dSn_inside*insideWeight_n;

        // partial derivative of the nonwetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dI[conti0EqIdx+1][saturationIdx] += tij*dpc_dSn_inside*up_n;
    }

    /*!
     * \brief Adds Robin flux derivatives for wetting and nonwetting phase
     *
     * \param derivativeMatrices The matrices containing the derivatives
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    { /* TODO maybe forward to problem for the user to implement the Robin derivatives?*/ }
};

} // end namespace Dumux

#endif
