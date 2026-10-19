// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Coupling manager of the small-deformation Cryer benchmark.
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_COUPLING_MANAGER_HH
#define DUMUX_TEST_POROELASTIC_CRYER_COUPLING_MANAGER_HH

#include <algorithm>
#include <memory>
#include <vector>

#include <dune/common/indices.hh>

#include <dumux/common/properties.hh>
#include <dumux/poromechanics/poroelastic/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Coupling manager of the small-deformation Cryer benchmark.
 *
 * Permeability and fluid density do not depend on the deformation. The flow residual of an element
 * therefore depends on the displacement only through the porosity in its storage term, that is, on
 * the displacement of its own vertices, and the transmissibilities do not change with the displacement.
 */
template<class MDTraits>
class CryerCouplingManager : public PoroMechanicsCouplingManager<MDTraits>
{
    using ParentType = PoroMechanicsCouplingManager<MDTraits>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using GridDiscretization = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridDiscretization<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndexType = typename GridView<id>::IndexSet::IndexType;

    static constexpr std::size_t flowIdx_ = ParentType::pmFlowId;
    static constexpr std::size_t mechIdx_ = ParentType::poroMechId;

public:
    using ParentType::couplingStencil;
    using ParentType::updateCoupledVariables;
    using ParentType::evalCouplingResidual;

    void init(std::shared_ptr<Problem<flowIdx_>> pmFlowProblem,
              std::shared_ptr<Problem<mechIdx_>> poroMechanicalProblem,
              const typename MDTraits::SolutionVector& curSol)
    {
        ParentType::init(pmFlowProblem, poroMechanicalProblem, curSol);

        const auto& flowGridDisc = this->problem(this->pmFlowId).gridGeometry();
        const auto& mechGridDisc = this->problem(this->poroMechId).gridGeometry();
        static constexpr int dim = GridView<flowIdx_>::dimension;
        stencils_.assign(flowGridDisc.gridView().size(0), {});
        for (const auto& element : elements(flowGridDisc.gridView()))
        {
            auto& stencil = stencils_[flowGridDisc.elementMapper().index(element)];
            for (int i = 0; i < element.subEntities(dim); ++i)
                stencil.push_back(mechGridDisc.vertexMapper().subIndex(element, i, dim));
            std::sort(stencil.begin(), stencil.end());
        }
    }

    const typename ParentType::template CouplingStencilType<flowIdx_>&
    couplingStencil(Dune::index_constant<flowIdx_>,
                    const Element<flowIdx_>& element,
                    Dune::index_constant<mechIdx_>) const
    { return stencils_[this->problem(this->pmFlowId).gridGeometry().elementMapper().index(element)]; }

    template<class PMFlowLocalAssembler, class UpdatableFluxVarCache>
    void updateCoupledVariables(Dune::index_constant<flowIdx_> pmFlowDomainId,
                                const PMFlowLocalAssembler& pmFlowLocalAssembler,
                                ElementVolumeVariables<flowIdx_>& elemVolVars,
                                UpdatableFluxVarCache&)
    {
        elemVolVars.bind(pmFlowLocalAssembler.element(),
                         pmFlowLocalAssembler.fvGeometry(),
                         this->curSol(pmFlowDomainId));
    }

    template<class PMFlowLocalAssembler>
    auto evalCouplingResidual(Dune::index_constant<flowIdx_>,
                              const PMFlowLocalAssembler& pmFlowLocalAssembler,
                              Dune::index_constant<mechIdx_>,
                              GridIndexType<mechIdx_>)
    {
        return pmFlowLocalAssembler.localResidual().evalStorage(pmFlowLocalAssembler.fvGeometry(),
                                                                pmFlowLocalAssembler.curElemVolVars());
    }

private:
    std::vector<typename ParentType::template CouplingStencilType<flowIdx_>> stencils_;
};

} // end namespace Dumux

#endif
