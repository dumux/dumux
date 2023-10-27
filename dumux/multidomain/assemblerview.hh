// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Subdomain-specific views on multidomain assemblers.
 */
#ifndef DUMUX_MULTIDOMAIN_ASSEMBLER_VIEW_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLER_VIEW_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/type_traits.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Subdomain-specific view on a multidomain assembler.
 *        Allows retrieval of sub-domain specific objects w/o passing a domain id.
 * \todo This is not necessarily fv-specifiv (could be in other header).
 * \todo Can we get rid of some of the interfaces?
 */
template<typename MDAssembler, std::size_t domainId>
class MultiDomainAssemblerSubDomainView
{
    static constexpr Dune::index_constant<domainId> myId{};

    template<class A>
    using HasStaticIsImplicitCheck = decltype(A::isImplicit());

    template<class A>
    static constexpr bool hasStaticIsImplicit = Dune::Std::is_detected<HasStaticIsImplicitCheck, A>::value;

public:
    using CouplingManager = typename MDAssembler::CouplingManager;
    using SolutionVector = typename MDAssembler::SolutionVector;

    MultiDomainAssemblerSubDomainView(MDAssembler& assembler, Dune::index_constant<domainId>)
    : assembler_{assembler}
    {}

    template<std::size_t i>
    auto localResidual(Dune::index_constant<i> id) const { return assembler_.localResidual(id); }
    auto localResidual() const { return assembler_.localResidual(myId); }

    template<std::size_t i>
    const auto& problem(Dune::index_constant<i> id) const { return assembler_.problem(id); }
    const auto& problem() const { return assembler_.problem(myId); }

    template<std::size_t i>
    const auto& gridGeometry(Dune::index_constant<i> id) const { return assembler_.gridGeometry(id); }
    const auto& gridGeometry() const { return assembler_.gridGeometry(myId); }

    template<std::size_t i>
    const auto& gridVariables(Dune::index_constant<i> id) const { return assembler_.gridVariables(id); }
    const auto& gridVariables() const { return assembler_.gridVariables(myId); }

    const auto& prevSol() const { return assembler_.prevSol(); }
    bool isStationaryProblem() const { return assembler_.isStationaryProblem(); }

    template<class A = MDAssembler, typename std::enable_if_t<hasStaticIsImplicit<A>, int> = 0>
    static constexpr bool isImplicit() { return MDAssembler::isImplicit(); }

    template<class A = MDAssembler, typename std::enable_if_t<!hasStaticIsImplicit<A>, int> = 0>
    bool isImplicit() const { return assembler_.isImplicit(); }

private:
    MDAssembler& assembler_;
};

} // end namespace Dumux

#endif
