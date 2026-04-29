// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Coupling manager for coupling two domains both discretized with CVFE schemes
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_COUPLINGMANAGER_BINARY_CVFE_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_COUPLINGMANAGER_BINARY_CVFE_HH

#include <utility>
#include <memory>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <cassert>

#include <dune/common/exceptions.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include <dumux/multidomain/boundary/couplingmapper_boundaryfaces.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Coupling manager for coupling two domains both discretized with CVFE schemes
 */
template<class MDTraits>
class CouplingManagerBinaryCvfe
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

public:
    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;
private:
    static_assert(MDTraits::numSubDomains == 2, "CouplingManagerBinaryCvfe requires exactly two subdomains");

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridView = typename GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;

    static_assert(DiscretizationMethods::isCVFE<std::decay_t<decltype(GridGeometry<0>::discMethod)>>,
        "CouplingManagerBinaryCvfe requires subdomain 0 to be discretized with a CVFE scheme");
    static_assert(DiscretizationMethods::isCVFE<std::decay_t<decltype(GridGeometry<1>::discMethod)>>,
        "CouplingManagerBinaryCvfe requires subdomain 1 to be discretized with a CVFE scheme");

    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables<id>>;
    template<std::size_t id> using Variables = Concept::Variables_t<GridVariables<id>>;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;

    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    using CouplingMapper = MapperCoupledMatchingBoundaryFaces<MDTraits, CouplingManagerBinaryCvfe<MDTraits>>;
    using CouplingStencil = typename CouplingMapper::CouplingStencil;

    template<class ElementSolution, std::size_t domainIdx, std::size_t coupledDomainIdx>
    struct CouplingContext
    {
        FVElementGeometry<coupledDomainIdx> fvGeometry;
        const Problem<coupledDomainIdx>& problem;
        ElementSolution elemSol;

        template<class ScvOrIpData>
        auto vars(const ScvOrIpData& scvOrIpData) const
        {
            Variables<coupledDomainIdx> variables;
            if constexpr (Concept::FVGridVariables<GridVariables<coupledDomainIdx>>)
                variables.update(elemSol, problem, fvGeometry.element(), scvOrIpData);
            else
                variables.update(elemSol, problem, fvGeometry, ipData(fvGeometry, scvOrIpData));
            return variables;
        }

        template<class ScvOrIpData>
        auto operator[](const ScvOrIpData& scvOrIpData) const
        {
            return vars(scvOrIpData);
        }
    };

public:

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<domainIdx<0>()>> domain0Problem,
              std::shared_ptr<Problem<domainIdx<1>()>> domain1Problem,
              GridVariablesTuple&& gridVariables,
              const SolutionVectorStorage& curSol)
    {
        this->setSubProblems(std::make_tuple(domain0Problem, domain1Problem));
        gridVariables_ = std::move(gridVariables);
        this->attachSolution(curSol);
        couplingMapper_.update(*this);
    }

    // \}


    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual (called from the local assembler)
     */
    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler) const
    {}

    /*!
     * \brief Update the coupling context
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
        // couplingContext() recomputes variables on-the-fly via makeCouplingContext_
        // using the correct element from the coupling mapper, so no cache update needed here.
    }

    // \}

    /*!
     * \brief Access the coupling context for domain i
     */
    template<std::size_t i, class FaceIpData>
    auto couplingContext(Dune::index_constant<i> domainI,
                         const FVElementGeometry<i>& fvGeometry,
                         const FaceIpData& faceIpData) const
    { return makeCouplingContext_(domainI, fvGeometry, faceIpData); }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief Coupling stencil
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    // \}

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param fvGeometry the finite volume element geometry
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const FVElementGeometry<i>& fvGeometry,
                   const typename FVElementGeometry<i>::SubControlVolumeFace& scvf) const
    {
        return isCoupled(domainI, fvGeometry, fvGeometry.intersectionIndex(scvf));
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const typename FVElementGeometry<i>::SubControlVolume& scv) const
    { return couplingMapper_.isCoupledDof(domainI, scv.dofIndex()); }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param fvGeometry the finite volume element geometry
     * \param iIdx the index of the intersection
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const FVElementGeometry<i>& fvGeometry,
                   std::size_t iIdx) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(fvGeometry.element());
        return couplingMapper_.isCoupled(domainI, eIdx, iIdx);
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param dofIdx the global dof index
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   std::size_t dofIdx) const
    { return couplingMapper_.isCoupledDof(domainI, dofIdx); }

private:
    template<class FVElementGeometry, class FaceIpData>
    std::size_t intersectionIndex_(const FVElementGeometry& fvGeometry,
                                   const FaceIpData& faceIpData) const
    {
        if constexpr (requires { faceIpData.scvfIndex(); })
        {
            const auto& scvf = fvGeometry.scvf(faceIpData.scvfIndex());
            return fvGeometry.intersectionIndex(scvf);
        }
        else if constexpr (requires { faceIpData.boundaryFaceIndex(); })
            return fvGeometry.boundaryFace(faceIpData.boundaryFaceIndex()).intersectionIndex();
        else
            DUNE_THROW(Dune::InvalidStateException, "FaceIpData must provide either scvfIndex() or intersectionIndex()");
    }

    template<std::size_t i, class FaceIpData>
    auto makeCouplingContext_(Dune::index_constant<i> domainI,
                              const FVElementGeometry<i>& fvGeometry,
                              const FaceIpData& faceIpData) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(fvGeometry.element());
        const auto iIdx = intersectionIndex_(fvGeometry, faceIpData);

        if (!couplingMapper_.isCoupled(domainI, eIdx, iIdx))
            DUNE_THROW(Dune::InvalidStateException, "No coupling context found at element << " << eIdx << " << with intersection index of " << iIdx);

        const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, eIdx, iIdx);
        constexpr auto domainJ = Dune::index_constant<1-i>();
        const auto& otherGridGeometry = this->problem(domainJ).gridGeometry();
        const auto& otherElement = otherGridGeometry.element(otherElementIdx);
        auto otherFvGeometry = localView(otherGridGeometry).bind(otherElement);
        const auto elemSol = elementSolution(otherElement, this->curSol(domainJ), otherGridGeometry);
        return CouplingContext<decltype(elemSol), i, 1-i>{ std::move(otherFvGeometry), this->problem(domainJ), std::move(elemSol) };
    }

    CouplingMapper couplingMapper_;
    GridVariablesTuple gridVariables_;
};

} // end namespace Dumux

#endif
