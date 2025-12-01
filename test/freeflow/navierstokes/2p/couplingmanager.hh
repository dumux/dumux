// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_CVFE_TWOP_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_CVFE_TWOP_HH

#include <dumux/multidomain/freeflow/couplingmanager_cvfe.hh>
#include <dumux/discretization/evalgradients.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for free flow systems
 * \note coupling manager for control volume finite element schemes
 */
template<class Traits>
class CVFEFreeFlowCouplingManagerTwoP
: public CVFEFreeFlowCouplingManager<Traits>
{
    using ParentType = CVFEFreeFlowCouplingManager<Traits>;
    using Scalar = typename Traits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename Traits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using VolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::VolumeVariables>;

public:
    static constexpr auto freeFlowMomentumIndex = typename Traits::template SubDomain<0>::Index();
    static constexpr auto freeFlowMassIndex = typename Traits::template SubDomain<1>::Index();

    // this can be used if the coupling manager is used inside a meta-coupling manager (e.g. multi-binary)
    // to manager the solution vector storage outside this class
    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;

    /*!
     * \brief Returns the pressure at a given interpolation point
     */
    template <class IpData>
    auto gradPhaseField(const Element<freeFlowMomentumIndex>& element,
                        const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                        const IpData& ipData) const
    {
        const auto& gg = this->problem(freeFlowMassIndex).gridGeometry();
        const auto& sol = this->curSol(freeFlowMassIndex);
        const auto elemSol = elementSolution(element, sol, gg);
        return evalGradientsAtLocalPos(element, element.geometry(), gg, elemSol, ipData.local())[VolumeVariables<freeFlowMassIndex>::Indices::phaseFieldIdx];
    }
};

namespace Detail {

// declaration (specialize for different discretization types)
template<class Traits, class DiscretizationMethod = typename Detail::MomentumDiscretizationMethod<Traits>::type>
struct CouplingManagerSupportsMultithreadedAssemblySelectorTwoP;

// multi-threading is not supported because we have only one coupling context instance and a mutable cache
template<class Traits, class D>
struct CouplingManagerSupportsMultithreadedAssemblySelectorTwoP<Traits, DiscretizationMethods::CVFE<D>>
{ using type = std::false_type; };

} // end namespace Detail

//! whether we support multithreaded assembly
template<class T>
struct CouplingManagerSupportsMultithreadedAssembly<CVFEFreeFlowCouplingManagerTwoP<T>>
: public Detail::CouplingManagerSupportsMultithreadedAssemblySelectorTwoP<T>::type
{};

} // end namespace Dumux

#endif
