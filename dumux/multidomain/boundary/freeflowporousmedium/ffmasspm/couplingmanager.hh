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
 * \copydoc Dumux::FreeFlowMassPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/boundary/darcydarcy/couplingmapper.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class FreeFlowMassPorousMediumCouplingManager
: public CouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ParentType = CouplingManager<MDTraits>;

public:
    static constexpr auto freeFlowMassIndex = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto porousMediumIndex = typename MDTraits::template SubDomain<1>::Index();

private:

    // obtain the type tags of the sub problems
    using FreeFlowMassTypeTag = typename MDTraits::template SubDomain<freeFlowMassIndex>::TypeTag;
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

    using VelocityVector = typename Element<freeFlowMassIndex>::Geometry::GlobalCoordinate;

    struct FreeFlowMassCouplingContext
    {
        VolumeVariables<porousMediumIndex> volVars;
        FVElementGeometry<porousMediumIndex> fvGeometry;
        std::size_t freeFlowMassScvfIdx;
        std::size_t porousMediumScvfIdx;
        mutable VelocityVector velocity; // velocity needs to be set externally, not availabe in this class
    };

    struct PorousMediumCouplingContext
    {
        VolumeVariables<freeFlowMassIndex> volVars;
        FVElementGeometry<freeFlowMassIndex> fvGeometry;
        std::size_t porousMediumScvfIdx;
        std::size_t freeFlowMassScvfIdx;
        mutable VelocityVector velocity; // velocity needs to be set externally, not availabe in this class
    };

    using CouplingMapper = DarcyDarcyBoundaryCouplingMapper<MDTraits>; // TODO rename/generalize class

public:

    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<porousMediumIndex>> darcyProblem,
              typename ParentType::SolutionVectorTuple& curSol)
    {
        this->attachSolution(curSol);


        freeFlowMassProblemPtr_ = &(*freeFlowMassProblem);
        porousMediumProblemPtr_ = &(*darcyProblem);
        // if (Dune::FloatCmp::ne(stokesProblem->gravity(), darcyProblem->spatialParams().gravity({})))
        //     DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        std::cout << problem(freeFlowMassIndex).name() << std::endl;
        std::cout << problem(porousMediumIndex).name() << std::endl;

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
            if (isCoupled(domainI, scvf))
            {
                const auto otherElementIdx = couplingMapper_.outsideElementIndex(domainI, scvf);
                constexpr auto domainJ = Dune::index_constant<1-domainI>();
                const auto& otherGridGeometry = problem(domainJ).gridGeometry();
                const auto& otherElement = otherGridGeometry.element(otherElementIdx);
                auto otherFvGeometry = localView(otherGridGeometry);
                otherFvGeometry.bindElement(otherElement);

                // there is only one scv for TPFA
                const auto& otherScv = [&]
                {
                    for (const auto scv : scvs(otherFvGeometry))
                        return scv;
                    DUNE_THROW(Dune::InvalidStateException, "Error");
                }();

                context.push_back({volVars(domainJ, otherElement, otherScv),
                                   otherFvGeometry,
                                   scvf.index(),
                                   couplingMapper_.flipScvfIndex(domainI, scvf),
                                   VelocityVector{}});
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

        for(const auto& context : contexts)
        {
            const auto expectedScvfIdx = domainI == freeFlowMassIndex ? context.freeFlowMassScvfIdx : context.porousMediumScvfIdx;
            if(scvf.index() == expectedScvfIdx)
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
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto eIdx = problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    // \}

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI, const SubControlVolumeFace<i>& scvf) const
    { return couplingMapper_.isCoupled(domainI, scvf); }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledElement(Dune::index_constant<i> domainI, std::size_t eIdx) const
    { return couplingMapper_.isCoupledElement(domainI, eIdx); }


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
        if constexpr (i == freeFlowMassIndex)
        {
            assert(freeFlowMassProblemPtr_);
            return *freeFlowMassProblemPtr_;
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
        if constexpr (i == freeFlowMassIndex)
        {
            assert(freeFlowMassProblemPtr_);
            return *freeFlowMassProblemPtr_;
        }
        else
        {
            assert(porousMediumProblemPtr_);
            return *porousMediumProblemPtr_;
        }
    }

private:

    mutable std::tuple<std::vector<FreeFlowMassCouplingContext>, std::vector<PorousMediumCouplingContext>> couplingContext_;
    mutable std::array<std::size_t, 2> couplingContextBoundForElement_;

    Problem<freeFlowMassIndex>* freeFlowMassProblemPtr_ = nullptr;
    Problem<porousMediumIndex>* porousMediumProblemPtr_ = nullptr;

    CouplingMapper couplingMapper_;
};

} // end namespace Dumux

#endif
