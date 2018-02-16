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
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
 * \brief The interface of the coupling manager for multi domain problems
 */

#ifndef DUMUX_STAGGERED_COUPLING_MANAGER_HH
#define DUMUX_STAGGERED_COUPLING_MANAGER_HH

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
 * \brief Base coupling manager for the staggered discretization.
 */
template<class MDTraits, class Implementation>
class StaggeredCouplingManagerBase: public CouplingManager<MDTraits, Implementation>
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template PrimaryVariables<id>;

    using StaggeredSubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<0>;
    using GridView = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using LocalResidual = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, LocalResidual);
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

public:

    using Traits = MDTraits;

    static constexpr auto cellCenterIdx = Dune::index_constant<0>();
    static constexpr auto faceIdx = Dune::index_constant<1>();

    template<class... Args>
    StaggeredCouplingManagerBase(Args&&... args) {}


    void init(std::shared_ptr<const Problem<0>> problem)
    {
        problemTuple_ = std::make_tuple(problem, problem);
    }

    template<class... Args>
    void init(Args&&... args)
    {
        problemTuple_ = std::make_tuple(args...);
    }

    void init(typename Traits::ProblemTuple&& problemTuple)
    {
        problemTuple_ = std::move(problemTuple);
    }

    template<class... Args>
    void updateSolution(Args&&... args) {}

    template<std::size_t i, std::size_t j, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               Dune::index_constant<j> domainJ,
                               const std::size_t dofIndex,
                               const PrimaryVariables<j>& priVars,
                               const Assembler& assembler)
    {
        DUNE_THROW(Dune::NotImplemented, "No update of coupling context implemented in base coupling manager!");
    }


    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainId,
                             const Element& element,
                             const Assembler& assembler) {}

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    const CouplingStencil& couplingStencil(const Element& element,
                                           Dune::index_constant<0> domainI,
                                           Dune::index_constant<1> domainJ) const
    {
        const auto& connectivityMap = problem(domainI).fvGridGeometry().connectivityMap();
        const auto eIdx = problem(domainI).fvGridGeometry().elementMapper().index(element);
        return connectivityMap(domainI, domainJ, eIdx);
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(const Element& element,
                                           Dune::index_constant<i> domainI,
                                           Dune::index_constant<j> domainJ) const
    {
        DUNE_THROW(Dune::NotImplemented, "No coupling stencils for additional domains implemented in base coupling manager!");
        return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    const CouplingStencil& couplingStencil(const SubControlVolumeFace& scvf,
                                           Dune::index_constant<1> domainI,
                                           Dune::index_constant<0> domainJ) const
    {
        const auto& connectivityMap = problem(domainI).fvGridGeometry().connectivityMap();
        return connectivityMap(domainI, domainJ, scvf.index());
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(const SubControlVolumeFace& scvf,
                                           Dune::index_constant<i> domainI,
                                           Dune::index_constant<j> domainJ) const
    {
        DUNE_THROW(Dune::NotImplemented, "No coupling stencils for additional domains implemented in base coupling manager!");
        return emptyStencil_;
    }

    template<class LocalAssembler, std::size_t i, std::size_t j>
    auto evalCouplingResidual(const LocalAssembler& assembler,
                              Dune::index_constant<i> domainI,
                              Dune::index_constant<j> domainJ)
    {
        DUNE_THROW(Dune::NotImplemented, "No default coupling residual implemented");
        return assembler.evalLocalResidual();
    }

    template<class LocalAssembler, std::size_t j>
    auto evalCouplingResidual(const LocalAssembler& assembler,
                              Dune::index_constant<0> domainI,
                              Dune::index_constant<j> domainJ)
    {
        return assembler.evalLocalResidualForCellCenter();
    }

    template<class LocalAssembler, std::size_t j>
    auto evalCouplingResidual(const LocalAssembler& assembler,
                              const SubControlVolumeFace& scvf,
                              Dune::index_constant<1> domainI,
                              Dune::index_constant<j> domainJ)
    {
        return assembler.evalLocalResidualForFace(scvf);
    }

    template<class LocalAssembler, std::size_t i, std::size_t j>
    auto evalCouplingResidual(const LocalAssembler& assembler,
                              const SubControlVolumeFace& scvf,
                              Dune::index_constant<i> domainI,
                              Dune::index_constant<j> domainJ)
    {
        DUNE_THROW(Dune::NotImplemented, "No default coupling residual implemented");
        return assembler.evalLocalResidualForFace(scvf);
    }

    //! Return a reference to the problem
    template<std::size_t id>
    const Problem<id>& problem(Dune::index_constant<id> domainIdx) const
    {
        assert(std::get<id>(problemTuple_) && "No problem set. Call init() first!");
        return *std::get<id>(problemTuple_);
    }

private:
    typename Traits::ProblemTuple problemTuple_;

    std::vector<std::size_t> emptyStencil_;
};

template<class MDTraits>
class StaggeredCouplingManager : public StaggeredCouplingManagerBase<MDTraits, StaggeredCouplingManager<MDTraits>>
{
    using ParentType = StaggeredCouplingManagerBase<MDTraits, StaggeredCouplingManager<MDTraits>>;

public:
    using ParentType::ParentType;

};

} //end namespace Dumux

#endif
