// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_CVFE_CVFE_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_CVFE_CVFE_HH

#include <deque>

#include <dune/common/exceptions.hh>

#include <dumux/common/concepts/variables_.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include "couplingmanager_base.hh"
#include "couplingconditions_cvfe_cvfe.hh"

namespace Dumux {


/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for coupling freeflow and porous medium flow models
 */
template<class MDTraits>
class FreeFlowPorousMediumCouplingManagerCvfe
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
    template<std::size_t id> using GridVariables =  GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables<id>>;
    template<std::size_t id> using ElementVariables = typename GridVariablesCache<id>::LocalView;

    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    using SolutionVector = typename MDTraits::SolutionVector;

    using MomentumDiscretizationMethod = typename GridGeometry<ParentType::freeFlowMomentumIndex>::DiscretizationMethod;

    using CouplingConditions = FFPMCouplingConditionsCvfe<MDTraits, FreeFlowPorousMediumCouplingManagerCvfe<MDTraits>>;

public:
    static constexpr auto freeFlowMomentumIndex = ParentType::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = ParentType::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = ParentType::porousMediumIndex;

    /*!
     * \brief Returns the mass flux across the coupling boundary.
     */
    template<std::size_t i, std::size_t j>
    auto massCouplingCondition(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ,
                               const FVElementGeometry<i>& fvGeometry,
                               const SubControlVolumeFace<i>& scvf,
                               const ElementVariables<i>& elemVars) const
    {
        static_assert(domainI != freeFlowMomentumIndex && domainJ != freeFlowMomentumIndex);

        const auto& couplingContext = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) {
            return cm.couplingContext(ii, fvGeometry, ipData(fvGeometry, scvf));
        });

        const auto& outsideFvGeometry = couplingContext.fvGeometry;
        const auto& ffMassFvGeometry = [&]() -> const auto& {
            if constexpr (i == porousMediumIndex)
                return outsideFvGeometry;
            else
                return fvGeometry;
        }();

        using VelocityVector = typename GridGeometry<freeFlowMomentumIndex>::GlobalCoordinate;
        VelocityVector intVelocity(0.0);
        // We apply the velocity-related quadrature rule
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf, typename GridGeometry<freeFlowMomentumIndex>::ScvfQuadratureRule{}))
        {
            const auto& faceIpData = qpData.ipData();
            const auto ffMassFaceIpData = ipData(ffMassFvGeometry, faceIpData.global());
            intVelocity += qpData.weight()
                         * this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).velocity(ffMassFvGeometry, ffMassFaceIpData);
        }

        // For the densities we interpolate at the interpolation point and do not use a quadrature rule
        const auto& ip = ipData(fvGeometry, scvf);
        const auto& outsideIp = ipData(outsideFvGeometry, ip.global());
        Scalar insideDensity = CouplingConditions::density(fvGeometry, elemVars, ip);
        Scalar outsideDensity = CouplingConditions::density(outsideFvGeometry, couplingContext, outsideIp);

        const auto normalIntVelocity = intVelocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = normalIntVelocity > 0.0;
        return CouplingConditions::advectiveFlux(insideDensity, outsideDensity, normalIntVelocity, insideIsUpstream);
    }

    //////////////////////// Conditions for FreeFlowMomentum - PorousMedium coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    template<class FaceIpData>
    Scalar pmPressure(Dune::index_constant<freeFlowMomentumIndex> domainI,
                      Dune::index_constant<porousMediumIndex> domainJ,
                      const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                      const ElementVariables<freeFlowMomentumIndex>& elemVars,
                      const FaceIpData& faceIpData) const
    {
        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
            domainI, fvGeometry, faceIpData
        );

        return CouplingConditions::pressure(context.fvGeometry, context, ipData(context.fvGeometry, faceIpData.global()));
    }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    template<class FaceIpData>
    auto darcyPermeability(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                           const FaceIpData& faceIpData) const
    {
        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
            Dune::index_constant<freeFlowMomentumIndex>(), fvGeometry, faceIpData
        );

        return CouplingConditions::darcyPermeability(fvGeometry, faceIpData, context);
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
     * \brief Returns the pressure at a given interpolation point
     */
    template <class IpData>
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const IpData& ipData,
                    const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).pressure(
            element, fvGeometry, ipData, considerPreviousTimeStep
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

    /*!
     * \brief Returns the density at a given interpolation point
     */
    template <class IpData>
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const IpData& ipData,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, fvGeometry, ipData, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the inside and outside densities at a given sub control volume face.
     */
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
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, fvGeometry, scv, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the effective viscosity at a given sub control volume face.
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
     * \brief Returns the effective viscosity at a interpolation point.
     */
    template <class IpData>
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const IpData& ipData,
                              const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).effectiveViscosity(
            element, fvGeometry, ipData, considerPreviousTimeStep
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
     * \brief Returns the velocity at an interpolation point
     */
    template <class IpData>
    auto velocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry,
                  const IpData& ipData) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).velocity(
            fvGeometry, ipData
        );
    }

    /*!
     * \brief Returns the velocity at the element center.
     */
    auto elementVelocity(const FVElementGeometry<freeFlowMassIndex>& fvGeometry) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).elementVelocity(
            fvGeometry
        );
    }

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    {
        if constexpr (MomentumDiscretizationMethod{} == DiscretizationMethods::fcdiamond)
        {
            // use coloring of the mass discretization for both domains
            // the diamond coloring is a subset (minimum amount of colors) of cctpfa/box coloring
            elementSets_ = computeColoring(this->problem(freeFlowMassIndex).gridGeometry()).sets;
        }
        else
        {
            // use coloring of the momentum discretization for both domains
            elementSets_ = computeColoring(this->problem(freeFlowMomentumIndex).gridGeometry()).sets;
        }

        // color the porous medium domain independently
        elementSetsPM_ = computeColoring(this->problem(porousMediumIndex).gridGeometry()).sets;
    }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainId the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        if constexpr (i == porousMediumIndex)
        {
            const auto& grid = this->problem(porousMediumIndex).gridGeometry().gridView().grid();
            for (const auto& elements : elementSetsPM_)
            {
                Dumux::parallelFor(elements.size(), [&](const std::size_t eIdx)
                {
                    const auto element = grid.entity(elements[eIdx]);
                    assembleElement(element);
                });
            }
        }
        else
        {
            const auto& grid = this->problem(freeFlowMomentumIndex).gridGeometry().gridView().grid();
            for (const auto& elements : elementSets_)
            {
                Dumux::parallelFor(elements.size(), [&](const std::size_t eIdx)
                {
                    const auto element = grid.entity(elements[eIdx]);
                    assembleElement(element);
                });
            }
        }
    }

private:
    std::deque<std::vector<ElementSeed<freeFlowMomentumIndex>>> elementSets_;
    std::deque<std::vector<ElementSeed<porousMediumIndex>>> elementSetsPM_;
};

template<class T>
struct CouplingManagerSupportsMultithreadedAssembly<FreeFlowPorousMediumCouplingManagerCvfe<T>>
: public std::true_type
{};

} // end namespace Dumux

#endif
