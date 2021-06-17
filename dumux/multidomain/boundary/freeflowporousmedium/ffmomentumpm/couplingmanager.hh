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
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMomentumPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include "couplingmapper.hh"

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class FreeFlowMomentumPorousMediumCouplingManager
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMomentumIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto porousMediumIndex = typename MDTraits::template SubDomain<1>::Index();

private:
    // obtain the type tags of the sub problems
    using FreeFlowMomentumTypeTag = typename MDTraits::template SubDomain<freeFlowMomentumIndex>::TypeTag;
    using PorousMediumTypeTag = typename MDTraits::template SubDomain<porousMediumIndex>::TypeTag;

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridView = typename GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using PrimaryVariables = typename MDTraits::template SubDomain<id>::PrimaryVariables;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;

    using VelocityVector = typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate;

    struct FreeFlowMomentumCouplingContext
    {
        VolumeVariables<porousMediumIndex> volVars;
        FVElementGeometry<porousMediumIndex> fvGeometry;
        std::size_t freeFlowMomentumScvfIdx;
        std::size_t porousMediumScvfIdx;
        VelocityVector gravity;
    };

    struct PorousMediumCouplingContext
    {
        std::size_t porousMediumScvfIdx;
        std::size_t freeFlowMomentumScvfIdx;
        std::size_t freeFlowMomentumDofIdx;
        VelocityVector faceVelocity;
    };

    using CouplingMapper = FreeFlowMomentumPorousMediumCouplingMapper<MDTraits, FreeFlowMomentumPorousMediumCouplingManager<MDTraits>>;

public:

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<porousMediumIndex>> porousMediumProblem,
              typename ParentType::SolutionVectorTuple& curSol)
    {
        this->attachSolution(curSol);

        freeFlowMomentumProblemPtr_ = &(*freeFlowMomentumProblem);
        porousMediumProblemPtr_ = &(*porousMediumProblem);
        // if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
        //     DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        couplingMapper_.update(*this);
    }

    // \}


    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    using ParentType::evalCouplingResidual;

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual (called from the local assembler)
     */
    template<std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler) const
    { bindCouplingContext(domainI, element); }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<std::size_t i>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element) const
    {
        auto fvGeometry = localView(problem(domainI).gridGeometry());
        fvGeometry.bindElement(element);
        bindCouplingContext(domainI, fvGeometry);
    }

    /*!
     * \brief prepares all data and variables that are necessary to evaluate the residual of an Darcy element (i.e. Stokes information)
     */
    template<std::size_t i>
    void bindCouplingContext(Dune::index_constant<i> domainI, const FVElementGeometry<i>& fvGeometry) const
    {
        auto& context = std::get<domainI>(couplingContext_);
        context.clear();

        const auto eIdx = problem(domainI).gridGeometry().elementMapper().index(fvGeometry.element());

        // do nothing if the element is not coupled to the other domain
        if (!isCoupledElement(domainI, eIdx))
            return;

        couplingContextBoundForElement_[domainI] = eIdx;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (couplingMapper_.isCoupled(domainI, scvf))
            {
                if constexpr (domainI == freeFlowMomentumIndex)
                {
                    const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                    constexpr auto domainJ = Dune::index_constant<1-domainI>();
                    const auto& otherGridGeometry = problem(domainJ).gridGeometry();
                    const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                    auto otherFvGeometry = localView(otherGridGeometry);
                    otherFvGeometry.bindElement(otherElement);

                    // there is only one scv for TPFA
                    context.push_back({volVars(domainJ, otherElement, *std::begin(scvs(otherFvGeometry))),
                                       otherFvGeometry,
                                       scvf.index(),
                                       couplingMapper_.flipScvfIndex(domainI, scvf),
                                       problem(domainJ).spatialParams().gravity(scvf.center())});
                }
                else
                {
                    const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                    constexpr auto domainJ = Dune::index_constant<1-domainI>();
                    const auto& otherGridGeometry = problem(domainJ).gridGeometry();
                    const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                    auto otherFvGeometry = localView(otherGridGeometry);
                    otherFvGeometry.bindElement(otherElement);

                    const auto otherScvfIdx = couplingMapper_.flipScvfIndex(domainI, scvf);
                    const auto& otherScvf = otherFvGeometry.scvf(otherScvfIdx);
                    const auto& otherScv = otherFvGeometry.scv(otherScvf.insideScvIdx());

                    context.push_back({scvf.index(),
                                       otherScvfIdx,
                                       otherScv.dofIndex(),
                                       faceVelocity(fvGeometry.element(), scvf)});
                }
            }
        }
    }

    /*!
     * \brief Update the coupling context for the Darcy residual w.r.t. Darcy DOFs
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
        auto& context = std::get<i>(couplingContext_);

        for (auto& c : context)
        {
            const auto& otherElement = c.fvGeometry.element();
            const auto otherElementIdx = problem(domainJ).gridGeometry().elementMapper().index(otherElement);

            if (otherElementIdx != dofIdxGlobalJ)
                continue;

            const auto elemSol = elementSolution(otherElement, this->curSol(domainJ), problem(domainJ).gridGeometry());
            for (const auto& scv : scvs(c.fvGeometry))
                c.volVars.update(elemSol, problem(domainJ), otherElement, scv);
        }
    }

    // \}

    /*!
     * \brief Access the coupling context needed for the Stokes domain
     */
    template<std::size_t i>
    const auto& couplingContext(Dune::index_constant<i> domainI,
                                const FVElementGeometry<i>& fvGeometry,
                                const SubControlVolumeFace<i> scvf) const
    {
        auto& contexts = std::get<i>(couplingContext_);

        if (contexts.empty() || couplingContextBoundForElement_[i] != scvf.insideScvIdx())
            bindCouplingContext(domainI, fvGeometry);

        for (const auto& context : contexts)
        {
            const auto expectedScvfIdx = domainI == freeFlowMomentumIndex ? context.freeFlowMomentumScvfIdx : context.porousMediumScvfIdx;
            if (scvf.index() == expectedScvfIdx)
                return context;
        }

        DUNE_THROW(Dune::InvalidStateException, "No coupling context found at scvf " << scvf.center());
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief The Stokes cell center coupling stencil w.r.t. Darcy DOFs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<porousMediumIndex> domainI,
                                           const Element<porousMediumIndex>& element,
                                           Dune::index_constant<freeFlowMomentumIndex> domainJ) const
    {
        const auto eIdx = problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given sub-control volume of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                           const Element<freeFlowMomentumIndex>& elementI,
                                           const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                           Dune::index_constant<porousMediumIndex> domainJ) const
    {
        return couplingMapper_.couplingStencil(domainI, elementI, scvI, domainJ);
    }

    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        return localAssemblerI.evalLocalResidual();
    }

    // \}

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupledLateralScvf(Dune::index_constant<freeFlowMomentumIndex> domainI, const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    { return  couplingMapper_.isCoupledLateralScvf(domainI, scvf); }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI, const SubControlVolumeFace<i>& scvf) const
    {
        if constexpr (i == freeFlowMomentumIndex)
            return couplingMapper_.isCoupled(domainI, scvf) || couplingMapper_.isCoupledLateralScvf(domainI, scvf);
        else
            return couplingMapper_.isCoupled(domainI, scvf);
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scv the sub control volume
     */
    bool isCoupled(Dune::index_constant<freeFlowMomentumIndex> domainI,
                   const SubControlVolume<freeFlowMomentumIndex>& scv) const
    { return couplingMapper_.isCoupled(domainI, scv); }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement(Dune::index_constant<i> domainI, std::size_t eIdx) const
    { return couplingMapper_.isCoupledElement(domainI, eIdx); }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<porousMediumIndex>& element,
                      const SubControlVolumeFace<porousMediumIndex>& scvf) const
    {
        // create a unit normal vector oriented in positive coordinate direction
        auto velocity = scvf.unitOuterNormal();
        using std::abs;
        std::for_each(velocity.begin(), velocity.end(), [](auto& v){ v = abs(v); });

        // create the actual velocity vector
        velocity *= this->curSol(freeFlowMomentumIndex)[couplingMapper_.outsideDofIndex(Dune::index_constant<porousMediumIndex>(), scvf)];

        return velocity;
    }


    //! Return the volume variables of domain i for a given element and scv
    template<std::size_t i>
    VolumeVariables<i> volVars(Dune::index_constant<i> domainI,
                               const Element<i>& element,
                               const SubControlVolume<i>& scv) const
    {
        VolumeVariables<i> volVars;
        const auto elemSol = elementSolution(element, this->curSol(domainI), problem(domainI).gridGeometry());
        volVars.update(elemSol, problem(domainI), element, scv);
        return volVars;
    }

    template<std::size_t i>
    Problem<i>& problem(Dune::index_constant<i> domainI)
    {
        if constexpr (i == freeFlowMomentumIndex)
        {
            assert(freeFlowMomentumProblemPtr_);
            return *freeFlowMomentumProblemPtr_;
        }
        else
        {
            assert(porousMediumProblemPtr_);
            return *porousMediumProblemPtr_;
        }
    }

    template<std::size_t i>
    const Problem<i>& problem(Dune::index_constant<i> domainI) const
    {
        if constexpr (i == freeFlowMomentumIndex)
        {
            assert(freeFlowMomentumProblemPtr_);
            return *freeFlowMomentumProblemPtr_;
        }
        else
        {
            assert(porousMediumProblemPtr_);
            return *porousMediumProblemPtr_;
        }
    }

private:

    mutable std::tuple<std::vector<FreeFlowMomentumCouplingContext>, std::vector<PorousMediumCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    Problem<freeFlowMomentumIndex>* freeFlowMomentumProblemPtr_ = nullptr;
    Problem<porousMediumIndex>* porousMediumProblemPtr_ = nullptr;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
