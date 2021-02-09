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
 * \ingroup PoroElastic
 * \brief Coupling manager for porous medium flow problems coupled to a poro-mechanical problem
 */

#ifndef DUMUX_POROMECHANICS_COUPLING_MANAGER_HH
#define DUMUX_POROMECHANICS_COUPLING_MANAGER_HH

#include <algorithm>
#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup PoroElastic
 * \brief Coupling manager for porous medium flow problems coupled to a poro-mechanical problem
 *
 *        Coupling manager for porous medium flow problems coupled to a poro-mechanical
 *        problem both defined on the same grid. Coupling occurs via the change of porosity
 *        and permeability due to mechanical deformations and the influence of the pore
 *        pressure on the effecive stresses acting on the porous medium.
 *
 * \tparam PMFlowId The porous medium flow domain id
 * \tparam PoroMechId The poro-mechanical domain id
 */
template< class MDTraits,
          std::size_t PMFlowId = 0,
          std::size_t PoroMechId = PMFlowId+1 >
class PoroMechanicsCouplingManager : public virtual CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using PrimaryVariables = typename GridVariables<id>::PrimaryVariables;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GridVolumeVariables<id>::VolumeVariables;
    template<std::size_t id> using GridGeometry = typename GridVariables<id>::GridGeometry;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;
    template<std::size_t id> using SubSolutionVector = typename GridVariables<id>::SolutionVector;
    template<std::size_t id>
    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element<id>>(),
                                                                  std::declval<SubSolutionVector<id>>(),
                                                                  std::declval<GridGeometry<id>>()))>;

    //! we assume that the two sub-problem operate on the same grid
    static_assert(std::is_same< GridView<PMFlowId>, GridView<PoroMechId> >::value,
                  "The grid types of the two sub-problems have to be equal!");

    //! this coupling manager is for cc - box only
    static_assert(GridGeometry<PoroMechId>::discMethod == DiscretizationMethod::box,
                  "Poro-mechanical problem must be discretized with the box scheme for this coupling manager!");

    static_assert(GridGeometry<PMFlowId>::discMethod == DiscretizationMethod::cctpfa ||
                  GridGeometry<PMFlowId>::discMethod == DiscretizationMethod::ccmpfa,
                  "Porous medium flow problem must be discretized with a cell-centered scheme for this coupling manager!");

    //! this does not work for enabled grid volume variables caching (update of local view in context has no effect)
    static_assert(!getPropValue<SubDomainTypeTag<PMFlowId>, Properties::EnableGridVolumeVariablesCache>(),
                  "Poromechanics framework does not yet work for enabled grid volume variables caching");

    //! Types used for coupling stencils
    template<std::size_t id>
    using CouplingIndexType = typename std::conditional< id == PMFlowId,
                                                         GridIndexType<PoroMechId>,
                                                         GridIndexType<PMFlowId> >::type;

    /*!
     * \brief Porous medium flow domain data required for the residual calculation of an
     *        element of the poro-mechanidal domain. We store the data required to do an
     *        update of all primary/secondary variables at the dof of the element.
     */
    struct PoroMechanicsCouplingContext
    {
        // We need optionals because the local views have no default constructor
        Element<PMFlowId> pmFlowElement;
        std::optional< FVElementGeometry<PMFlowId> > pmFlowFvGeometry;
        std::optional< ElementVolumeVariables<PMFlowId> > pmFlowElemVolVars;
    };

    /*!
     * \brief Poromechanical domain data required for the residual calculation of an
     *        element of the porous medium flow domain.
     */
    struct PorousMediumFlowCouplingContext
    {
        using Index = GridIndexType<PoroMechId>;
        using ElemSol = ElementSolution<PoroMechId>;
        using ElemDofs = std::vector<Index>;

        std::vector<Index> poroMechElemIndices;
        std::vector<ElemDofs> poroMechElemDofs;
        std::vector<ElemSol> poroMechElemSols;
    };

public:

    // export coupling context types
    template<std::size_t i>
    using CouplingContext = std::conditional_t<i == PMFlowId,
                                               PorousMediumFlowCouplingContext,
                                               PoroMechanicsCouplingContext>;

    // export the domain ids
    static constexpr auto pmFlowId = Dune::index_constant<PMFlowId>();
    static constexpr auto poroMechId = Dune::index_constant<PoroMechId>();

    /*!
     * \brief The types used for coupling stencils. An element of the poro-mechanical
     *        domain always only couples to the single dof (because we use a cell-centered
     *        scheme in the porous medium flow domain) of this same element.
     */
    template<std::size_t i, std::size_t j = (i == PMFlowId) ? PoroMechId : PMFlowId>
    using CouplingStencilType = typename std::conditional< i == PMFlowId,
                                                           std::vector< CouplingIndexType<i> >,
                                                           std::array< CouplingIndexType<i>, 1> >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param pmFlowProblem The porous medium flow problem
     * \param poroMechanicalProblem The poro-mechanical problem
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<PMFlowId> > pmFlowProblem,
              std::shared_ptr< Problem<PoroMechId> > poroMechanicalProblem,
              const SolutionVector& curSol)
    {
        // set the sub problems
        this->setSubProblem(pmFlowProblem, pmFlowId);
        this->setSubProblem(poroMechanicalProblem, poroMechId);

        // copy the solution vector
        ParentType::updateSolution(curSol);
        // set up the coupling map pmfow -> poromechanics
        initializeCouplingMap_();
    }

    /*!
     * \brief Return the coupling stencil for a given porous medium flow domain element
     */
    const CouplingStencilType<PMFlowId>& couplingStencil(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                                         const Element<PMFlowId>& element,
                                                         Dune::index_constant<PoroMechId> poroMechDomainId) const
    {
        return pmFlowCouplingMap_[ this->problem(pmFlowId).gridGeometry().elementMapper().index(element) ];
    }

    /*!
     * \brief Return the coupling element stencil for a given poro-mechanical domain element
     */
    const CouplingStencilType<PoroMechId> couplingStencil(Dune::index_constant<PoroMechId> poroMechDomainId,
                                                          const Element<PoroMechId>& element,
                                                          Dune::index_constant<PMFlowId> pmFlowDomainId) const
    {
        const auto eIdx = this->problem(pmFlowId).gridGeometry().elementMapper().index(element);
        return CouplingStencilType<PoroMechId>{ {eIdx} };
    }

    //! Construct the coupling context for an element of the porous medium flow domain
    template<class GridVariables>
    std::shared_ptr<PorousMediumFlowCouplingContext>
    makeCouplingContext(Dune::index_constant<PMFlowId> pmFlowDomainId,
                        const Element<PMFlowId>& element,
                        const GridVariables& gridVars) const
    {
        const auto& poroMechGridVars = gridVars[poroMechId];
        const auto& poroMechGridGeom = poroMechGridVars.gridGeometry();

        const auto globalI = gridVars[pmFlowId].gridGeometry().elementMapper().index(element);
        const auto& connectivityMapI = gridVars[pmFlowId].gridGeometry().connectivityMap()[globalI];

        std::vector< GridIndexType<PoroMechId> > elemIndices;
        std::vector< ElementSolution<PoroMechId> > elemSols;
        std::vector< std::vector<GridIndexType<PoroMechId>> > elemDofs;

        elemIndices.reserve(connectivityMapI.size()+1);
        elemSols.reserve(connectivityMapI.size()+1);
        elemDofs.reserve(connectivityMapI.size()+1);

        static constexpr int dim = GridView<PoroMechId>::dimension;
        {
            const auto& elementI = poroMechGridGeom.element(globalI);
            elemIndices.push_back(globalI);
            elemSols.push_back(elementSolution(elementI, poroMechGridVars.dofs(), poroMechGridGeom));

            elemDofs.push_back({});
            elemDofs.back().reserve(elementI.subEntities(dim));
            for (int i = 0; i < elementI.subEntities(dim); ++i)
                elemDofs.back().push_back(poroMechGridGeom.vertexMapper().subIndex(elementI, i, dim));
        }

        for (const auto& dataJ : connectivityMapI)
        {
            const auto& elementJ = poroMechGridGeom.element(dataJ.globalJ);
            elemIndices.push_back(dataJ.globalJ);
            elemSols.push_back(elementSolution(elementJ, poroMechGridVars.dofs(), poroMechGridGeom));

            elemDofs.push_back({});
            elemDofs.back().reserve(elementJ.subEntities(dim));
            for (int i = 0; i < elementJ.subEntities(dim); ++i)
                elemDofs.back().push_back(poroMechGridGeom.vertexMapper().subIndex(elementJ, i, dim));
        }

        using C = PorousMediumFlowCouplingContext;
        return std::make_shared<C>(C{std::move(elemIndices), std::move(elemDofs), std::move(elemSols)});
    }

    //! Construct the coupling context for an element of the poromechanical domain
    template<class GridVariables>
    std::shared_ptr<PoroMechanicsCouplingContext>
    makeCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainId,
                        const Element<PoroMechId>& element,
                        const GridVariables& gridVars) const
    {
        PoroMechanicsCouplingContext result;

        const auto& pmFlowGridVars = gridVars[Dune::index_constant<PMFlowId>()];
        const auto& pmFlowGridGeom = pmFlowGridVars.gridGeometry();
        const auto& pmFlowGridVolVars = pmFlowGridVars.gridVolVars();

        auto fvGeometry = localView( pmFlowGridGeom );
        auto elemVolVars = localView( pmFlowGridVolVars );
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, pmFlowGridVars.dofs());

        result.pmFlowElement = element;
        result.pmFlowFvGeometry.emplace(std::move(fvGeometry));
        result.pmFlowElemVolVars.emplace(std::move(elemVolVars));

        return std::make_shared<PoroMechanicsCouplingContext>(result);
    }

    void setCouplingContext(std::shared_ptr<PorousMediumFlowCouplingContext> context)
    { pmFlowCouplingContext_ = context; }

    void setCouplingContext(std::shared_ptr<PoroMechanicsCouplingContext> context)
    { poroMechCouplingContext_ = context; }

    // TODO: This would be called when assembling both domains, thus, in the case
    //       of assembly of a poromech element, we don't have the correct context!
    ElementSolution<PoroMechId> getPoroMechElemSol(const Element<pmFlowId>& element) const
    {
        const auto eIdx = this->problem(pmFlowId).gridGeometry().elementMapper().index(element);

        if (pmFlowCouplingContext_)
        {
            auto begin = pmFlowCouplingContext_->poroMechElemIndices.begin();
            auto end = pmFlowCouplingContext_->poroMechElemIndices.end();
            const auto it = std::find(begin, end, eIdx);

            if (it != end)
                return pmFlowCouplingContext_->poroMechElemSols[std::distance(begin, it)];
        }

        // make elem sol from stored solution.
        // TODO: How to resolve such circular dependencies?
        const auto poroMechElement = this->problem(poroMechId).gridGeometry().element(eIdx);
        return elementSolution(poroMechElement, this->curSol()[poroMechId], this->problem(poroMechId).gridGeometry());
    }

    /*!
     * \brief After deflection of the solution in the porous medium flow domain during
     *        element residual assembly in the poromechanics domain, we have to update
     *        the porous medium flow element variables of the coupling context
     */
    template< class PoroMechLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainId,
                               const PoroMechLocalAssembler& poroMechLocalAssembler,
                               Dune::index_constant<PMFlowId> pmFlowDomainId,
                               GridIndexType<PMFlowId> dofIdxGlobalJ,
                               const PrimaryVariables<PMFlowId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate the deflected pm flow domain primary variable
        // TODO: We still have to have and update the solution vector because it is used in updateCoupledVariables
        ParentType::updateCouplingContext(poroMechDomainId, poroMechLocalAssembler, pmFlowDomainId, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // now, update the coupling context (i.e. elemVolVars)
        const auto& element = poroMechCouplingContext_->pmFlowElement;
        const auto& fvGeometry = *(poroMechCouplingContext_->pmFlowFvGeometry);

        // TODO: Here the cm-internal solution is still used.
        //       This still works because the deflected solution is always the one of the current
        //       stage, and this function is only called to deflect the variables of the current stage.
        const auto& solution = this->curSol()[pmFlowId];
        poroMechCouplingContext_->pmFlowElemVolVars->bindElement(element, fvGeometry, solution);
    }

    /*!
     * \brief After deflection of the solution in the poromechanics domain during
     *        element residual assembly in that same domain, we have to update the
     *        porous medium flow element variables of the coupling context because
     *        the porosity/permeability might depend on the mechanical deformation
     */
    template< class PoroMechLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainIdI,
                               const PoroMechLocalAssembler& poroMechLocalAssembler,
                               Dune::index_constant<PoroMechId> poroMechDomainIdJ,
                               GridIndexType<PoroMechId> dofIdxGlobalJ,
                               const PrimaryVariables<PoroMechId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // // communicate the deflected displacement
        // TODO: We still have to have and update the solution vector because it is used in updateCoupledVariables
        ParentType::updateCouplingContext(poroMechDomainIdI, poroMechLocalAssembler, poroMechDomainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // now, update the coupling context (i.e. elemVolVars)
        // TODO: Here the cm-internal solution is still used.
        //       This still works because the deflected solution is always the one of the current
        //       stage, and this function is only called to deflect the variables of the current stage.
        const auto& solution = this->curSol()[pmFlowId];
        poroMechCouplingContext_->pmFlowElemVolVars->bindElement(poroMechCouplingContext_->pmFlowElement,
                                                                 *(poroMechCouplingContext_->pmFlowFvGeometry),
                                                                 solution);
    }

    /*!
     * \brief We need this overload to avoid ambiguity. However, for the porous medium flow
     *        domain weonly have to update the solution, which is done in the parent class.
     */
    template< std::size_t j, class PMFlowLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PMFlowId> pmFlowDomainId,
                               const PMFlowLocalAssembler& pmFlowLocalAssembler,
                               Dune::index_constant<j> domainIdJ,
                               GridIndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate the deflected displacement
        // TODO: We still have to have and update the solution vector because it is used in updateCoupledVariables
        ParentType::updateCouplingContext(pmFlowDomainId, pmFlowLocalAssembler, domainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // TODO: This is now much more expensive than the simple update before...
        //       Could an alternative be to store all stage solutions in the cm
        //       and not rely on a context here?
        if constexpr (j == PoroMechId)
        {
            for (int i = 0; i < pmFlowCouplingContext_->poroMechElemIndices.size(); ++i)
            {
                const auto& contextDofs = pmFlowCouplingContext_->poroMechElemDofs[i];
                const auto it = std::find(contextDofs.begin(), contextDofs.end(), dofIdxGlobalJ);

                if (it != contextDofs.end())
                {
                    const auto idxInContext = std::distance(contextDofs.begin(), it);
                    pmFlowCouplingContext_->poroMechElemSols[i][idxInContext][pvIdxJ] = priVarsJ[pvIdxJ];
                }
            }
        }
    }

    //! Pull up the base class' default implementation
    using ParentType::updateCoupledVariables;

    /*!
     * \brief Update the porous medium flow domain volume variables and flux variables cache
     *        after the coupling context has been updated. This has to be done because the
     *        mechanical deformation enters the porosity/permeability relationships.
     */
    template< class PMFlowLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                const PMFlowLocalAssembler& pmFlowLocalAssembler,
                                ElementVolumeVariables<PMFlowId>& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        // update the element volume variables to obtain the updated porosity/permeability
        elemVolVars.bind(pmFlowLocalAssembler.element(),
                         pmFlowLocalAssembler.fvGeometry(),
                         this->curSol()[pmFlowDomainId]);

        // update the transmissibilities subject to the new permeabilities
        elemFluxVarsCache.update(pmFlowLocalAssembler.element(),
                                 pmFlowLocalAssembler.fvGeometry(),
                                 elemVolVars);
    }

    /*!
     * \brief Update the poro-mechanics volume variables after the coupling context has been updated.
     *        This is necessary because the fluid density is stored in them and which potentially is
     *        solution-dependent.
     */
    template< class PoroMechLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(Dune::index_constant<PoroMechId> poroMechDomainId,
                                const PoroMechLocalAssembler& poroMechLocalAssembler,
                                ElementVolumeVariables<PoroMechId>& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        elemVolVars.bind(poroMechLocalAssembler.element(),
                         poroMechLocalAssembler.fvGeometry(),
                         this->curSol()[poroMechDomainId]);
    }

    /*!
     * \brief Evaluates the coupling element residual of the porous medium flow domain with
     *        respect to the poro-mechanical domain. The deformation might has an effect on
     *        both the permeability as well as the porosity. Thus, we need to compute fluxes
     *        and the storage term here.
     */
    template< class PMFlowLocalAssembler >
    auto
    evalCouplingResidual(Dune::index_constant<PMFlowId> pmFlowDomainId,
                         const PMFlowLocalAssembler& pmFlowLocalAssembler,
                         Dune::index_constant<PoroMechId> poroMechDomainId,
                         GridIndexType<PoroMechId> dofIdxGlobalJ)
    {
        // TODO: efficiency!
        return pmFlowLocalAssembler.evalLocalResidual();
    }

    /*!
     * \brief Evaluates the coupling element residual of the poromechanical domain with
     *        respect to the porous medium flow domain. The pressure has an effect on the
     *        mechanical stresses as well as the body forces. Thus, we have to compute
     *        the fluxes as well as the source term here.
     */
    template< class PoroMechLocalAssembler >
    auto
    evalCouplingResidual(Dune::index_constant<PoroMechId> poroMechDomainId,
                         const PoroMechLocalAssembler& poroMechLocalAssembler,
                         Dune::index_constant<PMFlowId> pmFlowDomainId,
                         GridIndexType<PMFlowId> dofIdxGlobalJ)
    {
        // TODO: efficiency!
        return poroMechLocalAssembler.evalLocalResidual();
    }

    //! Return the porous medium flow variables an element/scv of the poromech domain couples to
    const VolumeVariables<PMFlowId>& getPMFlowVolVars(const Element<PoroMechId>& element) const
    {
        //! If we do not yet have the queried object, build it first
        const auto eIdx = this->problem(poroMechId).gridGeometry().elementMapper().index(element);
        return (*poroMechCouplingContext_->pmFlowElemVolVars)[eIdx];
    }

    /*!
     * \brief the solution vector of the coupled problem
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    const SolutionVector& curSol() const
    { return ParentType::curSol(); }

private:
    /*!
     * \brief Initializes the pm flow domain coupling map. Since the elements
     *        of the poro-mechanical domain only couple to the same elements, we
     *        don't set up the maps here but return copy of the "stencil" always.
     */
    void initializeCouplingMap_()
    {
        // some references for convenience
        const auto& pmFlowGridGeom = this->problem(pmFlowId).gridGeometry();
        const auto& poroMechGridGeom = this->problem(poroMechId).gridGeometry();

        // make sure the two grids are really the same. Note that if the two grids
        // happen to have equal number of elements by chance, we don't detect this source of error.
        if (pmFlowGridGeom.gridView().size(0) != poroMechGridGeom.gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "The two sub-problems are assumed to operate on the same mesh!");

        pmFlowCouplingMap_.resize(pmFlowGridGeom.gridView().size(0));
        static constexpr int dim = GridView<PMFlowId>::dimension;
        for (const auto& element : elements(pmFlowGridGeom.gridView()))
        {
            const auto eIdx = pmFlowGridGeom.elementMapper().index(element);

            // firstly, the element couples to the nodal dofs in itself
            for (int i = 0; i < element.geometry().corners(); ++i)
                pmFlowCouplingMap_[eIdx].push_back( poroMechGridGeom.vertexMapper().subIndex(element, i , dim) );

            // the pm flow problem couples to the same elements as in its own stencil
            // due to the dependency of the residual on all permeabilities in its stencil,
            // which in turn depend on the mechanical deformations.
            const auto& inverseConnectivity = pmFlowGridGeom.connectivityMap()[eIdx];
            for (const auto& dataJ : inverseConnectivity)
                for (int i = 0; i < element.geometry().corners(); ++i)
                    pmFlowCouplingMap_[dataJ.globalJ].push_back( poroMechGridGeom.vertexMapper().subIndex(element, i , dim) );
        }

        // make stencils unique
        for (auto& stencil : pmFlowCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    // Container for storing the coupling element stencils for the pm flow domain
    std::vector< CouplingStencilType<PMFlowId> > pmFlowCouplingMap_;

    // the coupling contexts
    std::shared_ptr<PoroMechanicsCouplingContext> poroMechCouplingContext_ = nullptr;
    std::shared_ptr<PorousMediumFlowCouplingContext> pmFlowCouplingContext_ = nullptr;
};

} //end namespace Dumux

#endif
