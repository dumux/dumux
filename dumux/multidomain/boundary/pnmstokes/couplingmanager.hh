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
 * \ingroup BoubdaryCoupling
 * \ingroup BoxModel
 * \copydoc Dumux::PNMStokesCouplingManager
 */

#ifndef DUMUX_PNM_STOKES_COUPLINGMANAGER_HH
#define DUMUX_PNM_STOKES_COUPLINGMANAGER_HH

#include <iostream>
#include <string>
#include <utility>
#include <memory>

#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>

#include "couplingdata.hh"
#include "couplingmapper.hh"

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionBoundary
 * \brief Coupling manager for low-dimensional domains coupled at the boundary to the bulk
 *        domain.
 */
template<class MDTraits, bool useNeumannNeumannCoupling = true>
class PNMStokesCouplingManager
: public StaggeredCouplingManager<MDTraits>
{
    using Scalar = typename MDTraits::Scalar;
    using ThisType = PNMStokesCouplingManager<MDTraits, useNeumannNeumannCoupling>;
    using ParentType = StaggeredCouplingManager<MDTraits>;

public:
    static constexpr auto bulkFaceIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto bulkCellCenterIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto cellCenterIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto faceIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto bulkIdx = bulkCellCenterIdx;
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<2>::Index();

private:

    using SolutionVector = typename MDTraits::SolutionVector;
    using CouplingData = PNMStokesCouplingData<MDTraits, ThisType>;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridVolumeVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume  = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace  = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;
    template<std::size_t id> using ElementResidualVector = typename LocalResidual<id>::ElementResidualVector;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

    template<std::size_t id> using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;

    static constexpr auto lowDimDim = GridView<lowDimIdx>::dimension;
    static constexpr auto dimWorld = GridView<bulkCellCenterIdx>::dimensionworld;

    static_assert(lowDimDim == 1, "The bounding box coupling manager only works with one-dimensional low-dim grids");

    struct BulkCouplingContext
    {
        Element<lowDimIdx> element;
        FVElementGeometry<lowDimIdx> fvGeometry;
        ElementVolumeVariables<lowDimIdx> elemVolVars;
        ElementFluxVariablesCache<lowDimIdx> elemFluxVarsCache;
    };

    struct LowDimCouplingContext
    {
        Element<bulkIdx> element;
        FVElementGeometry<bulkIdx> fvGeometry;
        std::size_t scvfIdx;
        VolumeVariables<bulkIdx> volVars;
        Scalar velocity;

        const auto& getBulkScvf() const
        { return fvGeometry.scvf(scvfIdx); }
    };

public:

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;
    using ParentType::couplingStencil;
    using ParentType::updateCouplingContext;
    using ParentType::evalCouplingResidual;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    //! Initialize the coupling manager
    void init(std::shared_ptr<const Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<const Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        this->setSubProblem(bulkProblem, bulkCellCenterIdx);
        this->setSubProblem(bulkProblem, bulkFaceIdx);
        this->setSubProblem(lowDimProblem, lowDimIdx);

        this->curSol() = curSol;
        computeStencils();
        couplingData_ = std::make_unique<CouplingData>(*this);
    }

    //! Update after the grid has changed
    void update()
    { }

    // \}

    //! Update the solution vector before assembly
    void updateSolution(const SolutionVector& curSol)
    { this->curSol() = curSol; }

    //! Prepare the coupling stencils
    void computeStencils()
    {
        Dune::Timer timer;

        auto data = PNMStokesCouplingMapper::computeCouplingMapsAndStencils(*this);

        lowDimToBulkCellCenterCouplingStencils_ = std::move(data.lowDimToBulkCellCenterStencils);
        lowDimToBulkFaceCouplingStencils_ = std::move(data.lowDimToBulkFaceStencils);
        bulkCellCenterCouplingStencils_ = std::move(data.bulkCellCenterToLowDimStencils);
        bulkFaceCouplingStencils_ = std::move(data.bulkFaceToLowDimStencils);

        isCoupledLowDimDof_ = std::move(data.isCoupledLowDimDof);
        isCoupledBulkFaceDof_ = std::move(data.isCoupledBulkFaceDof);
        isCoupledBulkFrontalFaceDof_ = std::move(data.isCoupledBulkFrontalFaceDof);

        lowDimElementToBulkElementsMap_ = std::move(data.lowDimElementToBulkElementsMap);
        bulkElementToLowDimElementMap_ = std::move(data.bulkElementToLowDimElementMap);

        for(auto&& stencil : lowDimToBulkCellCenterCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : lowDimToBulkFaceCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : bulkCellCenterCouplingStencils_)
            removeDuplicates_(stencil.second);
        for(auto&& stencil : bulkFaceCouplingStencils_)
            removeDuplicates_(stencil.second);

        std::cout << "Computing coupling maps and stencils toook " << timer.elapsed() << " seconds." << std::endl;
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    template<class LocalAssemblerI, bool enable = !useNeumannNeumannCoupling, std::enable_if_t<enable, int> = 0>
    decltype(auto) evalCouplingResidual(Dune::index_constant<lowDimIdx> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<bulkCellCenterIdx> domainJ,
                                        std::size_t dofIdxGlobalJ)
    {
        ElementResidualVector<lowDimIdx> partialDerivs(localAssemblerI.element().subEntities(GridView<lowDimIdx>::dimension));
        partialDerivs = 0.0;

        for(const auto& scv : scvs(localAssemblerI.fvGeometry()))
        {
            if(isCoupledDof(lowDimIdx, scv.dofIndex()))
                partialDerivs[scv.indexInElement()][0] = this->curSol()[lowDimIdx][scv.dofIndex()] - couplingData().bulkPrivar(localAssemblerI.element(), scv);
        }

        return partialDerivs;
    }

    // ! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    // ! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    template<class LocalAssemblerI, bool enable = useNeumannNeumannCoupling, std::enable_if_t<enable, int> = 0>
    decltype(auto) evalCouplingResidual(Dune::index_constant<faceIdx> domainI,
                                        const SubControlVolumeFace<bulkIdx>& scvfI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<lowDimIdx> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        if (!isCoupledBulkFrontalFace(scvfI))
            return localAssemblerI.evalLocalResidualForFace(scvfI);

        using FaceResidualValue = typename LocalResidual<bulkIdx>::FaceResidualValue;
        FaceResidualValue partialDeriv(0.0);
        const auto& element = localAssemblerI.element();
        const auto& fvGeometry = localAssemblerI.fvGeometry();
        const auto& elemVolVars = localAssemblerI.curElemVolVars();
        const auto& elemFaceVars = localAssemblerI.curElemFaceVars();
        partialDeriv = couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvfI)
                       * scvfI.area() * elemVolVars[scvfI.insideScvIdx()].extrusionFactor();
        return partialDeriv;
    }


    //! Bind the bulk coupling context (i.e. lowDim information) TODO remove Assembler
    template<std::size_t i, class Assembler = int>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<bulkIdx>& element, const Assembler& assembler = 0) const
    {
        bulkCouplingContext_.clear();
        const auto ownElementIdx  = this->problem(bulkIdx).gridGeometry().elementMapper().index(element);
        boundBulkElemIdx_ = ownElementIdx;
        bool bindElement = false;
        std::size_t bulkElementIdx;

        // if the scvf belongs to an element that is directly coupled to a lowDim dof,
        // bind the element itself
        if (bulkCellCenterCouplingStencils_.count(ownElementIdx))
        {
            bindElement = true;
            bulkElementIdx = ownElementIdx;
        }

        // if we assemble another element that is not directly coupled to the lowDim dof
        // but shares a scvf with a neighbor element that does, bind that neighbor
        else
        {
            auto bulkFvGeometry = localView(this->problem(bulkIdx).gridGeometry());
            bulkFvGeometry.bind(element);
            for (const auto& scvf : scvfs(bulkFvGeometry))
            {
                if (bulkFaceCouplingStencils_.count(scvf.dofIndex()))
                {
                    bindElement = true;
                    bulkElementIdx = scvf.outsideScvIdx();
                }
            }
        }

        // do nothing if the element is not coupled to the other domain
        if (!bindElement)
            return;

        auto fvGeometry = localView(this->problem(lowDimIdx).gridGeometry());
        auto elemVolVars = localView(gridVariables(lowDimIdx).curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables(lowDimIdx).gridFluxVarsCache());

        const auto lowDimElemIdx = bulkElementToLowDimElementMap().at(bulkElementIdx);
        const auto& adjacentElement = this->problem(lowDimIdx).gridGeometry().boundingBoxTree().entitySet().entity(lowDimElemIdx);

        fvGeometry.bindElement(adjacentElement);
        elemVolVars.bind(adjacentElement, fvGeometry, this->curSol()[lowDimIdx]);
        elemFluxVarsCache.bind(adjacentElement, fvGeometry, elemVolVars);

        // add the coupling context
        bulkCouplingContext_.push_back({adjacentElement, fvGeometry, elemVolVars, elemFluxVarsCache});
    }

    //! Bind the coupling context for a low dim element TODO remove Assembler
    template<std::size_t i, class Assembler = int>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<lowDimIdx>& element, const Assembler& assembler = 0) const
    {
        lowDimCouplingContext_.clear();

        const auto lowDimElementIdx = this->problem(lowDimIdx).gridGeometry().elementMapper().index(element);
        boundLowDimElemIdx_ = lowDimElementIdx;

        // do nothing if the element is not coupled to the other domain
        if (!lowDimElementToBulkElementsMap().count(lowDimElementIdx))
            return;

        // prepare the coupling context
        auto fvGeometry = localView(this->problem(bulkIdx).gridGeometry());

        const auto& stencil = lowDimToBulkFaceCouplingStencils_.at(lowDimElementIdx);

        for (const auto& bulkElemIdx : lowDimElementToBulkElementsMap().at(lowDimElementIdx))
        {
            const auto& bulkElement = this->problem(bulkIdx).gridGeometry().boundingBoxTree().entitySet().entity(bulkElemIdx);
            fvGeometry.bindElement(bulkElement);

            const auto& cellCenterPriVars = this->curSol()[bulkCellCenterIdx][bulkElemIdx];
            using PrimaryVariablesType = typename VolumeVariables<bulkCellCenterIdx>::PrimaryVariables;
            PrimaryVariablesType priVars = makePriVarsFromCellCenterPriVars<PrimaryVariablesType>(cellCenterPriVars);
            const auto elemSol = elementSolution<FVElementGeometry<bulkCellCenterIdx>>(std::move(priVars));

            VolumeVariables<bulkIdx> bulkVolVars;
            for (const auto& scv : scvs(fvGeometry))
                bulkVolVars.update(elemSol, this->problem(bulkIdx), bulkElement, scv);

            // add the context
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (std::any_of(stencil.begin(), stencil.end(), [&](const auto x){ return scvf.dofIndex() == x; } ))
                {
                    const Scalar faceVelocity = this->curSol()[bulkFaceIdx][scvf.dofIndex()];
                    lowDimCouplingContext_.push_back({bulkElement, fvGeometry, scvf.index(), bulkVolVars, faceVelocity});
                    break; // only one scvf per element can be coupled to the lowDim dof
                }
            }
        }
    }

    /*!
     * \brief Update the coupling context for the lowDim residual w.r.t to the bulk cell center dofs
     */
    template<class LocalAssemblerI, class BulkCellCenterPriVars>
    void updateCouplingContext(Dune::index_constant<lowDimIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<bulkCellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const BulkCellCenterPriVars& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVars[pvIdxJ];

        for (auto& data : lowDimCouplingContext_)
        {
            const auto bulkElemIdx = this->problem(bulkIdx).gridGeometry().elementMapper().index(data.element);

            const auto& cellCenterPriVars = this->curSol()[bulkCellCenterIdx][bulkElemIdx];
            using PrimaryVariablesType = typename VolumeVariables<bulkCellCenterIdx>::PrimaryVariables;
            PrimaryVariablesType priVars = makePriVarsFromCellCenterPriVars<PrimaryVariablesType>(cellCenterPriVars);
            const auto elemSol = elementSolution<FVElementGeometry<bulkCellCenterIdx>>(std::move(priVars));

            for(const auto& scv : scvs(data.fvGeometry))
                data.volVars.update(elemSol, this->problem(bulkIdx), data.element, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the lowDim residual w.r.t to the bulk cell center dofs
     */
    template<class LocalAssemblerI, class BulkFacePriVars, bool enable = useNeumannNeumannCoupling, std::enable_if_t<enable, int> = 0>
    void updateCouplingContext(Dune::index_constant<lowDimIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<bulkFaceIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const BulkFacePriVars& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVars[pvIdxJ];

        for (auto& data : lowDimCouplingContext_)
        {
            const auto bulkScvf = data.getBulkScvf();
            if (bulkScvf.dofIndex() == dofIdxGlobalJ)
            {
                data.velocity = priVars;
                return; // only one context needs to be updated per bulk dof
            }
        }
    }

    /*!
     * \brief Update the coupling context for the bulk face residual w.r.t to the lowDim dofs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<bulkFaceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<lowDimIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<lowDimIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVars[pvIdxJ];

        for (auto& data : bulkCouplingContext_)
        {
            const auto elemSol = elementSolution(data.element, this->curSol()[lowDimIdx], this->problem(lowDimIdx).gridGeometry());

            for (const auto& scv : scvs(data.fvGeometry))
                getVolVarAccess_(domainJ, gridVars_(lowDimIdx).curGridVolVars(), data.elemVolVars, scv).update(elemSol, this->problem(lowDimIdx), data.element, scv);
        }
    }

    /*!
     * \brief Update the coupling context for the bulk face residual w.r.t to the lowDim dofs
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<bulkCellCenterIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<lowDimIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables<lowDimIdx>& priVars,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVars[pvIdxJ];

        for (auto& data : bulkCouplingContext_)
        {
            const auto elemSol = elementSolution(data.element, this->curSol()[lowDimIdx], this->problem(lowDimIdx).gridGeometry());

            for (const auto& scv : scvs(data.fvGeometry))
                getVolVarAccess_(domainJ, gridVars_(lowDimIdx).curGridVolVars(), data.elemVolVars, scv).update(elemSol, this->problem(lowDimIdx), data.element, scv);
        }
    }

    // \}

    //! Access the coupling data
    const CouplingData& couplingData() const
    { return *couplingData_; }

    //! Access the coupling context needed for the bulk domain
    const auto& bulkCouplingContext(const Element<bulkIdx>& element, const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        if (bulkCouplingContext_.empty() || boundBulkElemIdx_ != scvf.insideScvIdx())
            bindCouplingContext(bulkIdx, element);

        return bulkCouplingContext_;
    }

    //! Access the coupling context needed for the lowDim domain
    const auto& lowDimCouplingContext(const Element<lowDimIdx>& element, const SubControlVolume<lowDimIdx>& scv) const
    {
        if (lowDimCouplingContext_.empty() || boundLowDimElemIdx_ != scv.elementIndex())
            bindCouplingContext(lowDimIdx, element);

        return lowDimCouplingContext_;
    }

    /*!
     * \brief The coupling stencils
     */
    // \{

    /*!
     * \brief Return the coupling element stencil for a given bulk domain element
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<bulkCellCenterIdx> domainI,
                                           const Element<bulkIdx>& element,
                                           Dune::index_constant<lowDimIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (bulkCellCenterCouplingStencils_.count(eIdx))
            return bulkCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return the coupling element stencil for a given bulk domain element
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief Return the coupling element stencil for a given bulk domain element
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<lowDimIdx> domainI,
                                           const Element<lowDimIdx>& element,
                                           Dune::index_constant<bulkFaceIdx> domainJ) const
    {
        if (!useNeumannNeumannCoupling)
            return emptyStencil_;

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (lowDimToBulkFaceCouplingStencils_.count(eIdx))
            return lowDimToBulkFaceCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The coupling stencil of domain I, i.e. which domain J dofs
     *        the given domain I element's residual depends on.
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const SubControlVolumeFace<bulkIdx>& scvf,
                                           Dune::index_constant<j> domainJ) const
    { return emptyStencil_; }

    /*!
     * \brief The coupling stencil of a bulk face w.r.t to lowDim dofs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<bulkFaceIdx> domainI,
                                           const SubControlVolumeFace<bulkIdx>& scvf,
                                           Dune::index_constant<lowDimIdx> domainJ) const
    {
        const auto faceDofIdx = scvf.dofIndex();
        if(bulkFaceCouplingStencils_.count(faceDofIdx))
            return bulkFaceCouplingStencils_.at(faceDofIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return the coupling element stencil for a given bulk domain element
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<lowDimIdx> domainI,
                                           const Element<lowDimIdx>& element,
                                           Dune::index_constant<bulkCellCenterIdx> domainJ) const
    {
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (lowDimToBulkCellCenterCouplingStencils_.count(eIdx))
            return lowDimToBulkCellCenterCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    // \}

    //! Clear all internal data members
    void clear()
    { DUNE_THROW(Dune::NotImplemented, "clear() not implemented!"); }


    template<class IdType>
    const std::vector<std::size_t>& getAdditionalDofDependencies(IdType id, std::size_t bulkElementIdx) const
    { return emptyStencil_; }


    //! Returns whether a given bulk scvf is coupled to another domain
    template<std::size_t i>
    bool isCoupledEntity(Dune::index_constant<i> id, const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        return isCoupledBulkFaceDof_[scvf.dofIndex()];
    }

    //! Returns whether a given bulk scvf is directly coupled to another domain
    bool isCoupledBulkFrontalFace(const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        return isCoupledBulkFrontalFaceDof_[scvf.dofIndex()];
    }

    bool isCoupledDof(Dune::index_constant<2>, std::size_t lowDimDofIdx) const
    {
        return isCoupledLowDimDof_[lowDimDofIdx];
    }

    bool isCoupledDof(Dune::index_constant<lowDimIdx>, Dune::index_constant<bulkIdx>, std::size_t lowDimDofIdx) const
    {
        return isCoupledLowDimDof_[lowDimDofIdx];
    }

    const auto& lowDimElementToBulkElementsMap() const
    {
        return lowDimElementToBulkElementsMap_;
    }

    const auto& bulkElementToLowDimElementMap() const
    {
        return bulkElementToLowDimElementMap_;
    }

    /*!
     * \brief set the pointers to the grid variables
     * \param problems A tuple of shared pointers to the grid variables
     */
    void setGridVariables(GridVariablesTuple&& gridVariables)
    { gridVariables_ = gridVariables; }

    /*!
     * \brief set a pointer to one of the grid variables
     * \param problem a pointer to the grid variables
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

protected:

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

    //! Return a reference to an empty stencil
    std::vector<std::size_t>& emptyStencil()
    { return emptyStencil_; }

    void removeDuplicates_(std::vector<std::size_t>& stencil)
    {
        std::sort(stencil.begin(), stencil.end());
        stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
    }

private:

    template<std::size_t i>
    VolumeVariables<i>& getVolVarAccess_(Dune::index_constant<i> domainIdx, GridVolumeVariables<i>& gridVolVars, ElementVolumeVariables<i>& elemVolVars, const SubControlVolume<i>& scv)
    {
        if constexpr (getPropValue<SubDomainTypeTag<i>, Properties::EnableGridVolumeVariablesCache>())
            return gridVolVars.volVars(scv);
        else
            return elemVolVars[scv];
    }

    std::vector<bool> isCoupledLowDimDof_;
    std::vector<bool> isCoupledBulkFaceDof_;
    std::vector<bool> isCoupledBulkFrontalFaceDof_;

    std::unique_ptr<CouplingData> couplingData_;

    std::unordered_map<std::size_t, std::vector<std::size_t> > bulkCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > bulkFaceCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > lowDimToBulkCellCenterCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > lowDimToBulkFaceCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

    std::unordered_map<std::size_t, std::vector<std::size_t>> lowDimElementToBulkElementsMap_;
    std::unordered_map<std::size_t, std::size_t> bulkElementToLowDimElementMap_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    mutable std::vector<BulkCouplingContext> bulkCouplingContext_;
    mutable std::vector<LowDimCouplingContext> lowDimCouplingContext_;

    mutable std::size_t boundBulkElemIdx_;
    mutable std::size_t boundLowDimElemIdx_;

    /*!
     * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
     */
    GridVariablesTuple gridVariables_;
};

} // end namespace Dumux

#endif
