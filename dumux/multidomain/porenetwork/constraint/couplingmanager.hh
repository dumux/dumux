// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMConstraint
 * \brief Coupling manager for PNM with constraint at throat
 */

#ifndef DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_COUPLINGMANAGER_HH

#include <iostream>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/porenetwork/constraint/couplingmapper.hh>

namespace Dumux {

/*!
 * \ingroup PNMConstraint
 * \brief Coupling manager for PNM with constraint at throat
 */
template<class MDTraits>
class PNMConstraintCouplingManager
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

public:
    static constexpr auto poreNetworkIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto constraintIndex = typename MDTraits::template SubDomain<1>::Index();

private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using Problem = GetPropType<SubDomainTypeTag<i>, Properties::Problem>;
    template<std::size_t i> using PrimaryVariables = GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>;
    template<std::size_t i> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<i>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t i> using VolumeVariables = typename GetPropType<SubDomainTypeTag<i>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using FVElementGeometry = typename GridGeometry<i>::LocalView;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;
    template<std::size_t i> using Element = typename GridView<i>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    template<std::size_t i>
    static constexpr bool isBox()
    { return GridGeometry<i>::discMethod == DiscretizationMethods::box; }

    struct ConstraintCouplingContext
    {
        FVElementGeometry<poreNetworkIndex> fvGeometry;
        ElementVolumeVariables<poreNetworkIndex> elemVolVars;
    };

    using CouplingStencil = std::vector<std::size_t>;
public:
    //! export traits
    using MultiDomainTraits = MDTraits;
    //! export stencil types
    using CouplingStencils = std::unordered_map<std::size_t, CouplingStencil>;


    void init(std::shared_ptr<Problem<domainIdx<0>()>> problem0,
              std::shared_ptr<Problem<domainIdx<1>()>> problem1,
              const SolutionVector& curSol)
    {
        static_assert(isBox<0>(), "Only box PNM implemented!");

        this->setSubProblems(std::make_tuple(problem0, problem1));
        this->updateSolution(curSol);
        couplingMapper_.update(*this);
    }

    /*!
     * \brief set the pointers to the grid variables
     * \param gridVariables A tuple of shared pointers to the grid variables
     */
    void setGridVariables(GridVariablesTuple&& gridVariables)
    { gridVariables_ = gridVariables; }

    /*!
     * \brief set a pointer to one of the grid variables
     * \param gridVariables a pointer to the grid variables
     * \param domainIdx the domain index of the grid variables
     */
    template<class  GridVariables, std::size_t i>
    void setGridVariables(std::shared_ptr<GridVariables> gridVariables, Dune::index_constant<i> domainIdx)
    { std::get<i>(gridVariables_) = gridVariables; }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainIdx) const
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param element the coupled element of domain í
     * \param domainJ the domain index of domain j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual (called from the local assembler)
     */
    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler = 0) const
    {
        bindCouplingContext_(domainI, element);
    }

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
        // ToDo make this more efficient
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
        bindCouplingContext_(domainI, localAssemblerI.fvGeometry().element());
    }

    using ParentType::updateCoupledVariables;
    template<class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(Dune::index_constant<poreNetworkIndex> domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        elemFluxVarsCache.update(localAssemblerI.fvGeometry().element(), localAssemblerI.fvGeometry(), elemVolVars);
    }

    /*!
     * \brief Returns the entry pressure at a given throat, i.e. element
     */
    Scalar pcEntry(const Element<constraintIndex>& element) const
    {
        assert(couplingContext_.size() == 1);
        return this->problem(poreNetworkIndex).spatialParams().pcEntry(element, couplingContext_[0].elemVolVars);
    }

    /*!
     * \brief Returns the element volume variables for an element
     */
    const ElementVolumeVariables<poreNetworkIndex>& elemVolVars(const Element<constraintIndex>& element) const
    {
        assert(couplingContext_.size() == 1);
        return couplingContext_[0].elemVolVars;
    }

    /*!
     * \brief Returns the theta, i.e. the state indicator for a given element
     */
    Scalar theta(const Element<constraintIndex>& element) const
    {
        const auto eIdx = this->problem(constraintIndex).gridGeometry().elementMapper().index(element);
        return this->curSol(constraintIndex)[eIdx][0];
    }

private:
    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    const GridVariables<i>& gridVars_(Dune::index_constant<i> domainIdx) const
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    GridVariables<i>& gridVars_(Dune::index_constant<i> domainIdx)
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    void bindCouplingContext_(Dune::index_constant<poreNetworkIndex> domainI, const Element<poreNetworkIndex>& element) const
    {
        // Nothing to do for now
    }

    void bindCouplingContext_(Dune::index_constant<constraintIndex> domainI, const Element<constraintIndex>& element) const
    {
        couplingContext_.clear();
        auto fvGeometry = localView(this->problem(poreNetworkIndex).gridGeometry());
        fvGeometry.bind(element);

        auto elemVolVars = localView(gridVars_(poreNetworkIndex).curGridVolVars());
        elemVolVars.bind(fvGeometry.element(), fvGeometry, this->curSol(poreNetworkIndex));

        couplingContext_.push_back({fvGeometry, elemVolVars});
    }

    /*!
     * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
     */
    GridVariablesTuple gridVariables_;

    PNMConstraintCouplingMapper<MDTraits> couplingMapper_;
    // ToDo store in some container
    mutable std::vector<ConstraintCouplingContext> couplingContext_;
};

} // end namespace Dumux

#endif
