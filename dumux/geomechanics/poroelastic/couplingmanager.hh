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

#include <dune/common/indices.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/identityindexmap.hh>

namespace Dumux {

/*!
 * \brief Coupling manager for porous medium flow problems coupled to a poro-mechanical problem
 *
 *        Coupling manager for porous medium flow problems coupled to a poro-mechanical
 *        problem both defined on the same grid. Coupling occurs via the change of porosity
 *        and permeability due to mechanical deformations and the influence of the pore
 *        pressure on the effecive stresses acting on the porous medium.
 *
 * \tparam MDTraits The multidomain traits of the coupled problem
 * \tparam PMFlowId The porous medium flow domain id
 * \tparam PoroMechId The poro-mechanical domain id
 * \tparam ElementIndexMap Element index map between the two sub-domains. Allows
 *                         obtaining the element index the overlap element within
 *                         the other sub-domain. Must fulfill the interface of
 *                         Dumux::IdentityIndexMap.
 */
template< class MDTraits,
          std::size_t PMFlowId = 0,
          std::size_t PoroMechId = PMFlowId+1,
          class ElementIndexMap = IdentityIndexMap<Dune::index_constant<PMFlowId>, std::size_t,
                                                   Dune::index_constant<PoroMechId>, std::size_t> >
class PoroMechanicsCouplingManager : public virtual CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // evaluate if the default, static identity index map is used
    using DefaultIndexMap = IdentityIndexMap<Dune::index_constant<PMFlowId>, std::size_t,
                                             Dune::index_constant<PoroMechId>, std::size_t>;
    static constexpr bool isDefaultIndexMap = std::is_same<ElementIndexMap, DefaultIndexMap>::value;

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using PrimaryVariables = typename GridVariables<id>::PrimaryVariables;

    template<std::size_t id> using GridGeometry = typename GridVariables<id>::GridGeometry;
    template<std::size_t id> using ElementGeometry = typename GridGeometry<id>::LocalView;

    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;


    // static_assert(GridGeometry<PMFlowId>::discMethod == DiscretizationMethod::cctpfa ||
    //               GridGeometry<PMFlowId>::discMethod == DiscretizationMethod::ccmpfa,
    //               "Porous medium flow problem must be discretized with a cell-centered scheme for this coupling manager!");

    static constexpr bool flowUsesBox = GridGeometry<PMFlowId>::discMethod == DiscretizationMethod::box;
    static constexpr bool poroMechUsesBox = GridGeometry<PoroMechId>::discMethod == DiscretizationMethod::box;
    static constexpr bool isBoxFemCombo = flowUsesBox && !poroMechUsesBox;

    static_assert(poroMechUsesBox || GridGeometry<PoroMechId>::discMethod == DiscretizationMethod::fem,
        "Poro-mechanical sub-domain must be discretized with box or finite element scheme!");

    // Porous-medium flow sub-domain specific types
    using PMFlowVolumeVariables = typename GridVariables<PMFlowId>::GridVolumeVariables::VolumeVariables;
    using PMFlowElementVolumeVariables = typename GridVariables<PMFlowId>::GridVolumeVariables::LocalView;

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
        // We need unique ptrs because the local views have no default constructor
        Element<PMFlowId> pmFlowElement;
        std::unique_ptr< ElementGeometry<PMFlowId> > pmFlowFvGeometry;
        std::unique_ptr< PMFlowElementVolumeVariables > pmFlowElemVolVars;
    };

public:

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
                                                           std::vector< CouplingIndexType<i> > >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param pmFlowProblem The porous medium flow problem
     * \param poroMechanicalProblem The poro-mechanical problem
     * \param curSol The current solution
     * \note this overload is only valid if the default identity
     *       index map is used which does not require any state.
     */
    template< bool isDefault = isDefaultIndexMap, std::enable_if_t<isDefault, int> = 0 >
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
     * \brief Initialize the coupling manager.
     *
     * \param pmFlowProblem The porous medium flow problem
     * \param poroMechanicalProblem The poro-mechanical problem
     * \param curSol The current solution
     * \param indexMap The element index map between the sub-domains
     */
    void init(std::shared_ptr< Problem<PMFlowId> > pmFlowProblem,
              std::shared_ptr< Problem<PoroMechId> > poroMechanicalProblem,
              const SolutionVector& curSol,
              const ElementIndexMap& indexMap)
    {
        // make sure the two grids are really the same. Note that if the two grids
        // happen to have equal number of elements by chance, we don't detect this source of error.
        if (pmFlowProblem->gridGeometry().gridView().size(0) != poroMechanicalProblem->gridGeometry().gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "The two sub-problems are assumed to operate on the same mesh!");

        indexMap_ = indexMap;

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
        const auto eIdx = this->problem(poroMechId).gridGeometry().elementMapper().index(element);

        if (!flowUsesBox)
            return { {static_cast<GridIndexType<poroMechId>>(indexMap_.map(poroMechId, eIdx))} };
        return poroMechCouplingMap_[eIdx];
    }

    //! Pull up the base class' default implementation
    using ParentType::bindCouplingContext;

    /*!
     * \brief For the assembly of the element residual of an element of the poro-mechanics domain,
     *        we have to prepare the element variables of the porous medium flow domain.
     */
    template< class Assembler >
    void bindCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainId,
                             const Element<PoroMechId>& element,
                             const Assembler& assembler)
    {
        // TODO: How to realize this in an explicit way?
        static_assert(Assembler::isImplicit(), "Poro-mechanics currently require implicit time discretization");

        // if box-fem is used, the vol vars are built on the fly
        if (isBoxFemCombo)
            return;

        // otherwise, reset and build coupling vol vars
        poroMechCouplingContext_.pmFlowFvGeometry.reset(nullptr);
        poroMechCouplingContext_.pmFlowElemVolVars.reset(nullptr);

        // prepare the fvGeometry and the element volume variables
        // these quantities will be used later to obtain the effective pressure
        auto fvGeometry = localView( this->problem(pmFlowId).gridGeometry() );
        auto elemVolVars = localView( assembler.gridVariables(Dune::index_constant<PMFlowId>()).curGridVolVars() );

        const auto eIdx = this->problem(poroMechId).gridGeometry().elementMapper().index(element);
        const auto pmFlowElementIndex = indexMap_.map(poroMechId, eIdx);
        const auto pmFlowElement = this->problem(pmFlowId).gridGeometry().element(pmFlowElementIndex);

        fvGeometry.bindElement(pmFlowElement);
        elemVolVars.bindElement(pmFlowElement, fvGeometry, this->curSol()[pmFlowId]);

        poroMechCouplingContext_.pmFlowElement = pmFlowElement;
        poroMechCouplingContext_.pmFlowFvGeometry = std::make_unique< ElementGeometry<PMFlowId> >(fvGeometry);
        poroMechCouplingContext_.pmFlowElemVolVars = std::make_unique< PMFlowElementVolumeVariables >(elemVolVars);
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
        ParentType::updateCouplingContext(poroMechDomainId, poroMechLocalAssembler, pmFlowDomainId, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // maybe update the coupling context (i.e. elemVolVars)
        if (!isBoxFemCombo)
        {
            const auto& element = poroMechCouplingContext_.pmFlowElement;
            const auto& fvGeometry = *poroMechCouplingContext_.pmFlowFvGeometry;
            poroMechCouplingContext_.pmFlowElemVolVars->bindElement(element, fvGeometry, this->curSol()[pmFlowDomainId]);
        }
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
        // communicate the deflected displacement
        ParentType::updateCouplingContext(poroMechDomainIdI, poroMechLocalAssembler, poroMechDomainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // maybe update the coupling context (i.e. elemVolVars)
        if (!isBoxFemCombo)
            (*poroMechCouplingContext_.pmFlowElemVolVars).bindElement(poroMechCouplingContext_.pmFlowElement,
                                                                      *poroMechCouplingContext_.pmFlowFvGeometry,
                                                                      this->curSol()[Dune::index_constant<PMFlowId>()]);
    }

    /*!
     * \brief We need this overload to avoid ambiguity. However, for the porous medium flow
     *        domain we only have to update the solution, which is done in the parent class.
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
        ParentType::updateCouplingContext(pmFlowDomainId, pmFlowLocalAssembler, domainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    //! Pull up the base class' default implementation
    using ParentType::updateCoupledVariables;

    /*!
     * \brief Update the porous medium flow domain volume variables and flux variables cache
     *        after the coupling context has been updated. This has to be done because the
     *        mechanical deformation enters the porosity/permeability relationships.
     */
    template< class PMFlowLocalAssembler, class UpdatableFluxVarCache,
              bool isBox = flowUsesBox, std::enable_if_t<!isBox, int> = 0>
    void updateCoupledVariables(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                const PMFlowLocalAssembler& pmFlowLocalAssembler,
                                PMFlowElementVolumeVariables& elemVolVars,
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
     * \brief Update the porous medium flow domain volume variables and flux variables cache
     *        after the coupling context has been updated. This has to be done because the
     *        mechanical deformation enters the porosity/permeability relationships.
     */
    template< class PMFlowLocalAssembler, class UpdatableFluxVarCache,
              bool isBox = flowUsesBox, std::enable_if_t<isBox, int> = 0>
    void updateCoupledVariables(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                const PMFlowLocalAssembler& pmFlowLocalAssembler,
                                PMFlowElementVolumeVariables& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        // update the element volume variables to obtain the updated porosity/permeability
        elemVolVars.bind(pmFlowLocalAssembler.element(),
                         pmFlowLocalAssembler.fvGeometry(),
                         this->curSol()[pmFlowDomainId]);
    }

    /*!
     * \brief Update the poro-mechanics volume variables after the coupling context has been updated.
     *        This is necessary because the fluid density is stored in them and which potentially is
     *        solution-dependent.
     * \note This function is only valid and used in the case of the box scheme used in the mechanical domain
     */
    // template< class PoroMechLocalAssembler, class UpdatableFluxVarCache, bool isB = isBox, std::enable_if_t<isB, int> = 0 >
    // void updateCoupledVariables(Dune::index_constant<PoroMechId> poroMechDomainId,
    //                             const PoroMechLocalAssembler& poroMechLocalAssembler,
    //                             typename GridVariables<PoroMechId>::GridVolumeVariables::LocalView& elemVolVars,
    //                             UpdatableFluxVarCache& elemFluxVarsCache)
    // {
    //     elemVolVars.bind(poroMechLocalAssembler.element(),
    //                      poroMechLocalAssembler.fvGeometry(),
    //                      this->curSol()[poroMechDomainId]);
    // }

    /*!
     * \brief Evaluates the coupling element residual of the porous medium flow domain with
     *        respect to the poro-mechanical domain. The deformation might has an effect on
     *        both the permeability as well as the porosity. Thus, we need to compute fluxes
     *        and the storage term here.
     */
    template< class PMFlowLocalAssembler >
    typename LocalResidual<PMFlowId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<PMFlowId> pmFlowDomainId,
                         const PMFlowLocalAssembler& pmFlowLocalAssembler,
                         Dune::index_constant<PoroMechId> poroMechDomainId,
                         GridIndexType<PoroMechId> dofIdxGlobalJ)
    {
        // return pmFlowLocalAssembler.evalLocalResidual();
        auto res = pmFlowLocalAssembler.localResidual().evalFluxAndSource(pmFlowLocalAssembler.element(),
                                                                          pmFlowLocalAssembler.fvGeometry(),
                                                                          pmFlowLocalAssembler.curElemVolVars(),
                                                                          pmFlowLocalAssembler.elemFluxVarsCache(),
                                                                          pmFlowLocalAssembler.elemBcTypes());

        // If the residual instationary, evaluate storage
        if (!pmFlowLocalAssembler.localResidual().isStationary())
            res += pmFlowLocalAssembler.localResidual().evalStorage(pmFlowLocalAssembler.element(),
                                                                    pmFlowLocalAssembler.fvGeometry(),
                                                                    pmFlowLocalAssembler.prevElemVolVars(),
                                                                    pmFlowLocalAssembler.curElemVolVars());

        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of the poromechanical domain with
     *        respect to the porous medium flow domain. The pressure has an effect on the
     *        mechanical stresses as well as the body forces. Thus, we have to compute
     *        the fluxes as well as the source term here.
     * \note This overload is for the case of the box scheme being used in the mechanical domain
     */
    template< class PoroMechLocalAssembler, bool isB = poroMechUsesBox, std::enable_if_t<isB, int> = 0 >
    typename LocalResidual<PoroMechId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<PoroMechId> poroMechDomainId,
                         const PoroMechLocalAssembler& poroMechLocalAssembler,
                         Dune::index_constant<PMFlowId> pmFlowDomainId,
                         GridIndexType<PMFlowId> dofIdxGlobalJ)
    {
        return poroMechLocalAssembler.localResidual().evalFluxAndSource(poroMechLocalAssembler.element(),
                                                                        poroMechLocalAssembler.fvGeometry(),
                                                                        poroMechLocalAssembler.curElemVolVars(),
                                                                        poroMechLocalAssembler.elemFluxVarsCache(),
                                                                        poroMechLocalAssembler.elemBcTypes());
    }

    /*!
     * \brief Evaluates the coupling element residual of the poromechanical domain with
     *        respect to the porous medium flow domain. The pressure has an effect on the
     *        mechanical stresses as well as the body forces. Thus, we have to compute
     *        the fluxes as well as the source term here.
     * \note This overload is for the case of a finite element scheme being used in the mechanical domain
     */
    template< class PoroMechLocalAssembler, bool isB = poroMechUsesBox, std::enable_if_t<!isB, int> = 0 >
    typename LocalResidual<PoroMechId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<PoroMechId> poroMechDomainId,
                         const PoroMechLocalAssembler& poroMechLocalAssembler,
                         Dune::index_constant<PMFlowId> pmFlowDomainId,
                         GridIndexType<PMFlowId> dofIdxGlobalJ)
    {
        return poroMechLocalAssembler.localResidual().eval(poroMechLocalAssembler.element(),
                                                           poroMechLocalAssembler.feGeometry(),
                                                           poroMechLocalAssembler.curElemSol());
    }

    //! Return the coupling context (used in mechanical sub-problem to compute effective pressure)
    [[deprecated("Obtain the volume variables directly calling getPMFlowVolVars(element). Will be removed after 3.1!")]]
    const PoroMechanicsCouplingContext& poroMechanicsCouplingContext() const
    { return poroMechCouplingContext_; }


    //! Return the porous medium flow variables an element/scv of the poromech domain couples to
    template<bool isC = isBoxFemCombo, std::enable_if_t<!isC, int> = 0>
    const PMFlowVolumeVariables& getPMFlowVolVars(const Element<PoroMechId>& element,
                                                  const GlobalPosition<PoroMechId>& ip) const
    {
        //! If we do not yet have the queried object, build it first
        const auto eIdx = this->problem(poroMechId).gridGeometry().elementMapper().index(element);
        return (*poroMechCouplingContext_.pmFlowElemVolVars)[eIdx];
    }

    //! Return the porous medium flow variables an element/scv of the poromech domain couples to
    template<bool isC = isBoxFemCombo, std::enable_if_t<isC, int> = 0>
    const PMFlowVolumeVariables getPMFlowVolVars(const Element<PoroMechId>& element,
                                                 const GlobalPosition<PoroMechId>& ip) const
    {
        PMFlowVolumeVariables volVars;

        const auto eIdx = this->problem(poroMechId).gridGeometry().elementMapper().index(element);
        const auto pmFlowElementIndex = indexMap_.map(poroMechId, eIdx);
        const auto pmFlowElement = this->problem(pmFlowId).gridGeometry().element(pmFlowElementIndex);

        const auto& fvGridGeometry = this->problem(pmFlowId).gridGeometry();
        auto fvGeometry = localView(this->problem(pmFlowId).gridGeometry());
        fvGeometry.bindElement(pmFlowElement);

        // interpolate solution and set it for each entry in element solution
        auto elemSol = elementSolution(pmFlowElement, curSol()[pmFlowId], fvGridGeometry);
        const auto ipSol = evalSolution(pmFlowElement, pmFlowElement.geometry(), fvGridGeometry, elemSol, ip);
        for (unsigned int i = 0; i < fvGeometry.numScv(); ++i)
            elemSol[i] = ipSol;

        // Update volume variables with the interpolated solution. Note that this standard
        // implementation only works for element-wise constant parameters as we simply use
        // the first element scv for the vol var update. For heterogeneities within the element
        // or more complex models (e.g. 2p with interface solver) this doesn't work!
        volVars.update(elemSol, this->problem(pmFlowId), pmFlowElement, *scvs(fvGeometry).begin());
        return volVars;
    }

    /*!
     * \brief the solution vector of the coupled problem
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    const SolutionVector& curSol() const
    { return ParentType::curSol(); }

    const ElementIndexMap& bulkIndexMap() const
    { return indexMap_; }

private:
    /*!
     * \brief Initializes the pm flow domain coupling map. Since the elements
     *        of the poro-mechanical domain only couple to the same elements, we
     *        don't set up the maps here but return copy of the "stencil" always.
     * \note This overload is active if a cell-centered schemes is used in the flow domain
     */
    template<bool usesBox = flowUsesBox, std::enable_if_t<!usesBox, int> = 0>
    void initializeCouplingMap_()
    {
        // some references for convenience
        const auto& pmFlowGridGeom = this->problem(pmFlowId).gridGeometry();
        const auto& poroMechGridGeom = this->problem(poroMechId).gridGeometry();

        pmFlowCouplingMap_.resize(pmFlowGridGeom.gridView().size(0));
        for (const auto& element : elements(pmFlowGridGeom.gridView()))
        {
            const auto eIdx = pmFlowGridGeom.elementMapper().index(element);
            const auto poroMechElemIdx = indexMap_.map(pmFlowId, eIdx);
            const auto poroMechElement = poroMechGridGeom.element(poroMechElemIdx);
            const auto elemDofsPoroMech = getPoroMechElementDofs_(poroMechElement);

            // firstly, the element couples to the dofs in this element
            pmFlowCouplingMap_[eIdx].insert( pmFlowCouplingMap_[eIdx].end(),
                                             elemDofsPoroMech.begin(),
                                             elemDofsPoroMech.end() );

            // the pm flow problem couples to the same elements as in its own stencil
            // due to the dependency of the residual on all permeabilities in its stencil,
            // which in turn depend on the mechanical deformations.
            const auto& inverseConnectivity = pmFlowGridGeom.connectivityMap()[eIdx];
            for (const auto& dataJ : inverseConnectivity)
                pmFlowCouplingMap_[dataJ.globalJ].insert( pmFlowCouplingMap_[dataJ.globalJ].end(),
                                                          elemDofsPoroMech.begin(),
                                                          elemDofsPoroMech.end() );
        }

        // make stencils unique
        for (auto& stencil : pmFlowCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    /*!
     * \brief Initializes the pm flow domain coupling map. Since the elements
     *        of the poro-mechanical domain only couple to the same elements, we
     *        don't set up the maps here but return copy of the "stencil" always.
     * \note This overload is active if a cell-centered schemes is used in the flow domain
     */
    template<bool usesBox = flowUsesBox, std::enable_if_t<usesBox, int> = 0>
    void initializeCouplingMap_()
    {
        // some references for convenience
        const auto& pmFlowGridGeom = this->problem(pmFlowId).gridGeometry();
        const auto& poroMechGridGeom = this->problem(poroMechId).gridGeometry();

        pmFlowCouplingMap_.resize(pmFlowGridGeom.gridView().size(0));
        poroMechCouplingMap_.resize(poroMechGridGeom.gridView().size(0));
        for (const auto& element : elements(pmFlowGridGeom.gridView()))
        {
            const auto eIdx = pmFlowGridGeom.elementMapper().index(element);
            const auto poroMechElemIdx = indexMap_.map(pmFlowId, eIdx);
            const auto poroMechElement = poroMechGridGeom.element(poroMechElemIdx);

            // the elements couple to all dofs in the other domain
            const auto elemDofsPoroMech = getPoroMechElementDofs_(poroMechElement);
            const auto elemDofsFlow = getFlowElementDofs_(element);

            auto& poroMechStencil = poroMechCouplingMap_[poroMechElemIdx];
            auto& pmFlowStencil = pmFlowCouplingMap_[eIdx];

            poroMechStencil.insert(poroMechStencil.end(), elemDofsFlow.begin(), elemDofsFlow.end());
            pmFlowStencil.insert(pmFlowStencil.end(), elemDofsPoroMech.begin(), elemDofsPoroMech.end());
        }

        // make stencils unique
        for (auto& stencil : pmFlowCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
        for (auto& stencil : poroMechCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    //! extract dofs inside an element (box scheme)
    template<bool useBox = poroMechUsesBox, std::enable_if_t<useBox, int> = 0>
    std::vector<GridIndexType<poroMechId>> getPoroMechElementDofs_(const Element<poroMechId>& element) const
    {
        std::vector<GridIndexType<poroMechId>> result;

        auto fvGeometry = localView(this->problem(poroMechId).gridGeometry());
        fvGeometry.bindElement(element);

        result.reserve(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            result.push_back(scv.dofIndex());

        return result;
    }

    //! extract dofs inside an element (fem scheme)
    template<bool useBox = poroMechUsesBox, std::enable_if_t<!useBox, int> = 0>
    std::vector<GridIndexType<poroMechId>> getPoroMechElementDofs_(const Element<poroMechId>& element) const
    {
        std::vector<GridIndexType<poroMechId>> result;

        auto feGeometry = localView(this->problem(poroMechId).gridGeometry());
        feGeometry.bind(element);

        const auto& feBasisLocalView = feGeometry.feBasisLocalView();
        result.reserve(feBasisLocalView.size());
        for (unsigned int i = 0; i < feBasisLocalView.size(); ++i)
            result.push_back(feBasisLocalView.index(i));

        return result;
    }

    //! extract dofs inside an element (box scheme)
    std::vector<GridIndexType<poroMechId>> getFlowElementDofs_(const Element<pmFlowId>& element) const
    {
        std::vector<GridIndexType<pmFlowId>> result;

        auto fvGeometry = localView(this->problem(pmFlowId).gridGeometry());
        fvGeometry.bindElement(element);

        result.reserve(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            result.push_back(scv.dofIndex());

        return result;
    }

    // data members
    ElementIndexMap indexMap_;  //!< maps element indices between the domains
    std::vector< CouplingStencilType<PMFlowId> > pmFlowCouplingMap_; //!< coupling stencils
    std::vector< CouplingStencilType<PoroMechId> > poroMechCouplingMap_; //!< coupling stencils
    PoroMechanicsCouplingContext poroMechCouplingContext_; //!< coupling context
};

} //end namespace Dumux

#endif
