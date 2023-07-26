// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowPorousMediumCouplingManagerStaggeredCCTpfa
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_STAGGERED_CCTPFA_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_STAGGERED_CCTPFA_HH

#include "couplingmanager_base.hh"
#include "couplingconditions_staggered_cctpfa.hh"

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for coupling freeflow and porous medium flow models
 *        specialization for staggered-cctpfa coupling.
 */
template<class MDTraits>
class FreeFlowPorousMediumCouplingManagerStaggeredCCTpfa
: public FreeFlowPorousMediumCouplingManagerBase<MDTraits>
{
    using ParentType = FreeFlowPorousMediumCouplingManagerBase<MDTraits>;

    using Scalar = typename MDTraits::Scalar;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    using SolutionVector = typename MDTraits::SolutionVector;

    using CouplingConditions = FFPMCouplingConditionsStaggeredCCTpfa<MDTraits, FreeFlowPorousMediumCouplingManagerStaggeredCCTpfa<MDTraits>>;

public:
    static constexpr auto freeFlowMomentumIndex = ParentType::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = ParentType::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = ParentType::porousMediumIndex;

public:
    /*!
     * \brief Returns the mass flux across the coupling boundary.
     */
    template<std::size_t i, std::size_t j>
    auto massCouplingCondition(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ,
                               const FVElementGeometry<i>& fvGeometry,
                               const typename FVElementGeometry<i>::SubControlVolumeFace& scvf,
                               const ElementVolumeVariables<i>& elemVolVars) const
    {
        static_assert(domainI != freeFlowMomentumIndex && domainJ != freeFlowMomentumIndex);

        const auto& couplingContext = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingContext(ii, fvGeometry, scvf);
        });

        const auto& freeFlowElement = [&]
        {
            if constexpr (domainI == freeFlowMassIndex)
                return fvGeometry.element();
            else
                return couplingContext.fvGeometry.element();
        }();

        const auto& freeFlowScvf = [&]
        {
            if constexpr (domainI == freeFlowMassIndex)
                return scvf;
            else
                return couplingContext.fvGeometry.scvf(couplingContext.freeFlowMassScvfIdx);

        }();

        // todo revise velocity (see ff mom pm mgr)

        couplingContext.velocity = this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(freeFlowElement, freeFlowScvf);
        return CouplingConditions::massCouplingCondition(domainI, domainJ, fvGeometry, scvf, elemVolVars, couplingContext);
    }


    //////////////////////// Conditions for FreeFlowMomentum - PorousMedium coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    NumEqVector<freeFlowMomentumIndex> momentumCouplingCondition(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                                                 Dune::index_constant<porousMediumIndex> domainJ,
                                                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                                 const typename FVElementGeometry<freeFlowMomentumIndex>::SubControlVolumeFace& scvf,
                                                                 const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars) const
    {
        if (scvf.isLateral())
            return NumEqVector<freeFlowMomentumIndex>(0.0);

        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
            domainI, fvGeometry, scvf
        );

        return CouplingConditions::momentumCouplingCondition(fvGeometry, scvf, elemVolVars, context);
    }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    auto darcyPermeability(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                           const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        if (scvf.isFrontal())
        {
            const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
                Dune::index_constant<freeFlowMomentumIndex>(), fvGeometry, scvf
            );

            return CouplingConditions::darcyPermeability(fvGeometry, scvf, context);
        }
        else
        {
            const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
            const auto& orthogonalScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());
            const auto& frontalScvfOnBoundary = fvGeometry.frontalScvfOnBoundary(orthogonalScv);
            const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
                Dune::index_constant<freeFlowMomentumIndex>(), fvGeometry, frontalScvfOnBoundary
            );

            return CouplingConditions::darcyPermeability(fvGeometry, frontalScvfOnBoundary, context);
        }
    }

    //////////////////////// Conditions for FreeFlowMomentum - FreeFlowMass coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).pressure(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the pressure at the center of a sub control volume corresponding to a given sub control volume face.
     *        This is used for setting a Dirichlet pressure for the mass model when a fixed pressure for the momentum balance is set at another
     *        boundary. Since the the pressure at the given scvf is solution-dependent and thus unknown a priori, we just use the value
     *        of the interior cell here.
     */
    Scalar cellPressure(const Element<freeFlowMomentumIndex>& element,
                        const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).cellPressure(
            element, scvf
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    auto insideAndOutsideDensity(const Element<freeFlowMomentumIndex>& element,
                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                 const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                 const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).insideAndOutsideDensity(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, scv, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).effectiveViscosity(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<freeFlowMassIndex>& element,
                      const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(
            element, scvf
        );
    }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupledLateralScvf(Dune::index_constant<freeFlowMomentumIndex> domainI,
                              Dune::index_constant<porousMediumIndex> domainJ,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupledLateralScvf(ii, scvf);
        });
    }
};

} // end namespace Dumux

#endif
