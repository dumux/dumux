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
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief A multidomain assembler for Jacobian and residual contribution per element (staggered method)
 */
#ifndef DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/simpleassemblystructs.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all multidomain local assemblers (staggered)
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual assembler implementation
 * \tparam implicit whether the assembly is explicit or implicit in time
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, bool isImplicit = true>
class SubDomainStaggeredLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler,Implementation, isImplicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler,Implementation, isImplicit>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename Assembler::SolutionVector;
    using SubSolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;

    static constexpr auto numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();

    using SimpleMassBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMassBalanceSummands);
    using SimpleMomentumBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummands);
    using SimpleMomentumBalanceSummandsVector = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummandsVector);

public:
    static constexpr auto domainId = typename Dune::index_constant<id>();
    static constexpr auto cellCenterId = typename Dune::index_constant<0>();
    static constexpr auto faceId = typename Dune::index_constant<1>();

    using ParentType::ParentType;

    explicit SubDomainStaggeredLocalAssemblerBase(const Assembler& assembler,
                                                  const Element& element,
                                                  const SolutionVector& curSol,
                                                  CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.fvGridGeometry(domainId)),
                 localView(assembler.gridVariables(domainId).curGridVolVars()),
                 localView(assembler.gridVariables(domainId).prevGridVolVars()),
                 localView(assembler.gridVariables(domainId).gridFluxVarsCache()),
                 assembler.localResidual(domainId),
                 (element.partitionType() == Dune::GhostEntity))
    , curElemFaceVars_(localView(assembler.gridVariables(domainId).curGridFaceVars()))
    , prevElemFaceVars_(localView(assembler.gridVariables(domainId).prevGridFaceVars()))
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleCoefficientMatrixAndRHS(JacobianMatrixRow& coefficentMatrixRow, SubSol& RHS, GridVariablesTuple& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        assembleCoefficientMatrixAndRHSImpl_(domainId, coefficentMatrixRow, RHS, gridVariables);
    }

    /*!
     * \brief Assemble the residual only (for convergence information)
     */
    template<class SubSol>
    void assembleResidual(SubSol& res)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        assembleResidualImpl_(domainId, res);
    }

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    CellCenterResidualValue evalLocalResidualForCellCenter(SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            simpleMassBalanceSummands.setToZero(this->element(), this->fvGeometry());
            return CellCenterResidualValue(0.0);
        }

        return isImplicit ? evalLocalResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars(), simpleMassBalanceSummands)
                   : evalLocalResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars(), simpleMassBalanceSummands);
    }

    /*!
     * \brief Evaluates the complete local residual for the current cell center.
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    CellCenterResidualValue evalLocalResidualForCellCenter(const ElementVolumeVariables& elemVolVars,
                                                           const ElementFaceVariables& elemFaceVars,
                                                           SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        evalLocalFluxAndSourceResidualForCellCenter(elemVolVars, elemFaceVars, simpleMassBalanceSummands);

        if (!this->assembler().isStationaryProblem())
            evalLocalStorageResidualForCellCenter(simpleMassBalanceSummands);

        this->localResidual().evalBoundaryForCellCenter(this->problem(), this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), simpleMassBalanceSummands);

        CellCenterResidualValue res(0.0);

        for (auto&& scvf : scvfs(this->fvGeometry()))
        {
            res += simpleMassBalanceSummands.coefficients[scvf.localFaceIdx()] * elemFaceVars[scvf].velocitySelf();
        }
        res -= simpleMassBalanceSummands.RHS;

        return res;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    void evalLocalFluxAndSourceResidualForCellCenter(SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        isImplicit ? evalLocalFluxAndSourceResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars(), simpleMassBalanceSummands)
                          : evalLocalFluxAndSourceResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars(), simpleMassBalanceSummands);
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVolVars The element volume variables
     * \param elemVolVars The element face variables
     */
    void evalLocalFluxAndSourceResidualForCellCenter(const ElementVolumeVariables& elemVolVars, const ElementFaceVariables& elemFaceVars, SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        this->localResidual().evalFluxAndSourceForCellCenter(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), simpleMassBalanceSummands);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    void evalLocalStorageResidualForCellCenter(SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        this->localResidual().evalStorageForCellCenter(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars(), simpleMassBalanceSummands);
    }

    /*!
     * \brief Convenience function to evaluate the  local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf, SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            simpleMomentumBalanceSummands.setToZero(scvf);
            return FaceResidualValue(0.0);
        }



        return isImplicit ? evalLocalResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars(), simpleMomentumBalanceSummands)
                          : evalLocalResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars(), simpleMomentumBalanceSummands);
    }

    /*!
     * \brief Evaluates the complete local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        evalLocalFluxAndSourceResidualForFace(scvf, elemVolVars, elemFaceVars, simpleMomentumBalanceSummands);

        if (!this->assembler().isStationaryProblem())
            evalLocalStorageResidualForFace(scvf, simpleMomentumBalanceSummands);

        this->localResidual().evalBoundaryForFace(this->problem(), this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf, simpleMomentumBalanceSummands);

        FaceResidualValue res(0.0);

        res +=
        simpleMomentumBalanceSummands.selfCoefficient * elemFaceVars[scvf].velocitySelf()
        + simpleMomentumBalanceSummands.oppositeCoefficient * elemFaceVars[scvf].velocityOpposite();

        const int numSubFaces = scvf.pairData().size();
        for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx){
            res += simpleMomentumBalanceSummands.parallelCoefficients[localSubFaceIdx] * elemFaceVars[scvf].velocityParallel(localSubFaceIdx);
            res += simpleMomentumBalanceSummands.innerNormalCoefficients[localSubFaceIdx] * elemFaceVars[scvf].velocityNormalInside(localSubFaceIdx);
            res += simpleMomentumBalanceSummands.outerNormalCoefficients[localSubFaceIdx] * elemFaceVars[scvf].velocityNormalOutside(localSubFaceIdx);
        }

        res +=
        simpleMomentumBalanceSummands.pressureCoefficient * elemVolVars[scvf.insideScvIdx()].pressure()
        - simpleMomentumBalanceSummands.RHS;

        return res;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    void evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf, SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        isImplicit ? evalLocalFluxAndSourceResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars(), simpleMomentumBalanceSummands)
                          : evalLocalFluxAndSourceResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars(), simpleMomentumBalanceSummands);
    }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    void evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf,
                                                            const ElementVolumeVariables& elemVolVars,
                                                            const ElementFaceVariables& elemFaceVars,
                                                            SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        this->localResidual().evalFluxAndSourceForFace(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf, simpleMomentumBalanceSummands);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    void evalLocalStorageResidualForFace(const SubControlVolumeFace& scvf, SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        this->localResidual().evalStorageForFace(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars(), this->prevElemFaceVars(), this->curElemFaceVars(), scvf, simpleMomentumBalanceSummands);
    }

    const Problem& problem() const
    { return this->assembler().problem(domainId); }

    //! The current element volume variables
    ElementFaceVariables& curElemFaceVars()
    { return curElemFaceVars_; }

    //! The element volume variables of the provious time step
    ElementFaceVariables& prevElemFaceVars()
    { return prevElemFaceVars_; }

    //! The current element volume variables
    const ElementFaceVariables& curElemFaceVars() const
    { return curElemFaceVars_; }

    //! The element volume variables of the provious time step
    const ElementFaceVariables& prevElemFaceVars() const
    { return prevElemFaceVars_; }

    CouplingManager& couplingManager()
    { return couplingManager_; }

private:

    //! Assembles the residuals for the cell center dofs.
    template<class SubSol>
    void assembleResidualImpl_(Dune::index_constant<0>, SubSol& res)
    {
        SimpleMassBalanceSummands simpleMassBalanceSummands(this->element(), this->fvGeometry());

        const auto cellCenterGlobalI = problem().fvGridGeometry().elementMapper().index(this->element());
        res[cellCenterGlobalI] = this->asImp_().assembleCellCenterResidualImpl(simpleMassBalanceSummands);
    }

    //! Assembles the residuals for the face dofs.
    template<class SubSol>
    void assembleResidualImpl_(Dune::index_constant<1>, SubSol& res)
    {
        for(auto&& scvf : scvfs(this->fvGeometry()))
        {
            SimpleMomentumBalanceSummands simpleMomentumBalanceSummands(scvf);
            res[scvf.dofIndex()] +=  this->asImp_().assembleFaceResidualImpl(scvf, simpleMomentumBalanceSummands);
        }
    }

    //! Assembles the residuals and derivatives for the cell center dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleCoefficientMatrixAndRHSImpl_(Dune::index_constant<0>, JacobianMatrixRow& coefficentMatrixRow, SubSol& RHS, GridVariablesTuple& gridVariables)
    {
        SimpleMassBalanceSummands simpleMassBalanceSummands(this->element(), this->fvGeometry());

        const auto cellCenterGlobalI = problem().fvGridGeometry().elementMapper().index(this->element());

        this->asImp_().assembleCellCenterCoefficientMatrixAndRHSImpl(coefficentMatrixRow[domainId], *std::get<domainId>(gridVariables), simpleMassBalanceSummands);
        RHS[cellCenterGlobalI] = simpleMassBalanceSummands.RHS;

        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<Dune::Hybrid::size(coefficentMatrixRow), domainId>{};
        forEach(otherDomainIds, [&, domainId = domainId](auto&& domainJ)
        {
            this->asImp_().assembleCoefficientMatrixCellCenterCoupling(domainJ, coefficentMatrixRow[domainJ], *std::get<domainJ>(gridVariables), simpleMassBalanceSummands);
        });
    }

    //! Assembles the residuals and derivatives for the face dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleCoefficientMatrixAndRHSImpl_(Dune::index_constant<1>, JacobianMatrixRow& coefficentMatrixRow, SubSol& RHS, GridVariablesTuple& gridVariables)
    {
        SimpleMomentumBalanceSummandsVector simpleMomentumBalanceSummandsVector;
        for (auto&& scvf : scvfs(this->fvGeometry())){
            SimpleMomentumBalanceSummands zeroValue(scvf);
            zeroValue.setToZero(scvf);
            simpleMomentumBalanceSummandsVector.push_back(zeroValue);
        }

        this->asImp_().assembleFaceCoefficientMatrixAndRHSImpl(coefficentMatrixRow[domainId], *std::get<domainId>(gridVariables), simpleMomentumBalanceSummandsVector);

        for(auto&& scvf : scvfs(this->fvGeometry()))
        {
            RHS[scvf.dofIndex()] += simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].RHS;
        }

        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<Dune::Hybrid::size(coefficentMatrixRow), domainId>{};
        forEach(otherDomainIds, [&, domainId = domainId](auto&& domainJ)
        {
            this->asImp_().assembleCoefficientMatrixFaceCoupling(domainJ, coefficentMatrixRow[domainJ], *std::get<domainJ>(gridVariables), simpleMomentumBalanceSummandsVector);
        });
    }

    ElementFaceVariables curElemFaceVars_;
    ElementFaceVariables prevElemFaceVars_;
    CouplingManager& couplingManager_; //!< the coupling manager
};

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all implicit multidomain local assemblers (staggered)
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual assembler implementation
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainStaggeredLocalAssemblerImplicitBase : public SubDomainStaggeredLocalAssemblerBase<id, TypeTag, Assembler, Implementation>
{
    using ParentType = SubDomainStaggeredLocalAssemblerBase<id, TypeTag, Assembler, Implementation>;
    static constexpr auto domainId = Dune::index_constant<id>();
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        auto& couplingManager = this->couplingManager();
        const auto& element = this->element();
        const auto& curSol = this->curSol();
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& curElemFaceVars = this->curElemFaceVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);
        curElemVolVars.bind(element, fvGeometry, curSol);
        curElemFaceVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
        if (!this->assembler().isStationaryProblem())
        {
            this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol());
            this->prevElemFaceVars().bindElement(element, fvGeometry, this->assembler().prevSol());
        }
    }

};

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief The staggered multidomain local assembler
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 * \tparam implicit whether the assembler is explicit or implicit in time
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainStaggeredLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief Staggered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using FaceVariables = typename ElementFaceVariables::FaceVariables;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*ModelTraits::dim());
    static constexpr auto domainI = Dune::index_constant<id>();
    static constexpr auto cellCenterId = typename Dune::index_constant<0>();
    static constexpr auto faceId = typename Dune::index_constant<1>();

    static constexpr auto numEq = ModelTraits::numEq();
    static constexpr auto numEqCellCenter = CellCenterPrimaryVariables::dimension;
    static constexpr auto numEqFace = FacePrimaryVariables::dimension;

    using SimpleMassBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMassBalanceSummands);
    using SimpleMomentumBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummands);
    using SimpleMomentumBalanceSummandsVector = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummandsVector);

public:
    using ParentType::ParentType;

    CellCenterResidualValue assembleCellCenterResidualImpl(SimpleMassBalanceSummands& simpleMassBalanceSummands)
    {
        return this->evalLocalResidualForCellCenter(simpleMassBalanceSummands);
    }

    FaceResidualValue assembleFaceResidualImpl(const SubControlVolumeFace& scvf, SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
        return this->evalLocalResidualForFace(scvf, simpleMomentumBalanceSummands);
    }

    /*!
     * \brief Computes the coefficient matrix entries for the given element and adds them
     *        to the global matrix.
     *
     * \return The right hand side entries for the element at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    void assembleCellCenterCoefficientMatrixAndRHSImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables, SimpleMassBalanceSummands& simpleMassBalanceSummands)
    {
        assert(domainI == cellCenterId);

        this->evalLocalResidualForCellCenter(simpleMassBalanceSummands);

        // CCToCC Block is empty for SIMPLE algorithm
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    void assembleFaceCoefficientMatrixAndRHSImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables, SimpleMomentumBalanceSummandsVector& simpleMomentumBalanceSummandsVector)
    {
        assert(domainI == faceId);

        // get an alias for convenience
        const auto& fvGeometry = this->fvGeometry();

        //fill the coefficientMatrix and RHS
        for(auto&& scvf : scvfs(fvGeometry))
        {
            this->evalLocalResidualForFace(scvf, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()]);

            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            for(int pvIdx = 0; pvIdx < numEqFace; ++pvIdx)
            {
                //coefficient in the momentum balance equation of the selfVelocity
                updateGlobalJacobian_(A, faceGlobalI, scvf.dofIndex(), pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].selfCoefficient);

                //coefficient in the momentum balance equation of the oppositeVelocity
                updateGlobalJacobian_(A, faceGlobalI, scvf.dofIndexOpposingFace(), pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].oppositeCoefficient);

                //coefficients in the momentum balance equation of the parallelVelocities
                const int numSubFaces = scvf.pairData().size();
                for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                {
                    const auto& data = scvf.pairData(localSubFaceIdx);

                    updateGlobalJacobian_(A, faceGlobalI, data.normalPair.first, pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].innerNormalCoefficients[localSubFaceIdx]);
                    if(data.outerParallelFaceDofIdx >= 0)
                    {
                        updateGlobalJacobian_(A, faceGlobalI, data.outerParallelFaceDofIdx, pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].parallelCoefficients[localSubFaceIdx]);
                        //note on comparing this with StaggeredFreeFlowConnectivityMap computeFaceToCellCenterStencil_: SIMPLE does not consider the normalPair contributions, such terms are included in the right-hand side
                    }
                    if(!scvf.boundary()){
                        updateGlobalJacobian_(A, faceGlobalI, data.normalPair.second, pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].outerNormalCoefficients[localSubFaceIdx]);
                    }
                }
            }
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class GridVariables>
    void assembleCoefficientMatrixCellCenterCoupling(Dune::index_constant<1> domainJ, JacobianBlock& A, GridVariables& gridVariables, const SimpleMassBalanceSummands& simpleMassBalanceSummands)
    {
        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->problem().fvGridGeometry();
        // const auto& curSol = this->curSol()[domainI];
        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto cellCenterGlobalI = fvGridGeometry.elementMapper().index(element);

        for(const auto& scvfJ : scvfs(fvGeometry))
        {
             const auto globalJ = scvfJ.dofIndex();

            for(int pvIdx = 0; pvIdx < numEqFace; ++pvIdx)
            {
                //fill the coefficientMatrix
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, simpleMassBalanceSummands.coefficients[scvfJ.localFaceIdx()]);
            }
        }
    }

    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleCoefficientMatrixCellCenterCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                            const CellCenterResidualValue& res, GridVariables& gridVariables, SimpleMassBalanceSummands& simpleMassBalanceSummands)
    {
        //TODO adapt for SIMPLE
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class GridVariables>
    void assembleCoefficientMatrixFaceCoupling(Dune::index_constant<0> domainJ, JacobianBlock& A, GridVariables& gridVariables, const SimpleMomentumBalanceSummandsVector& simpleMomentumBalanceSummandsVector)
    {
        //TODO adapt for SIMPLE
    }

    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleCoefficientMatrixFaceCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A, GridVariables& gridVariables, SimpleMomentumBalanceSummandsVector& simpleMomentumBalanceSummandsVector)
    {
        // get an alias for convenience
        const auto& fvGeometry = this->fvGeometry();

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        for(auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            for(int pvIdx = 0; pvIdx < numEqCellCenter; ++pvIdx)
            {
                //note on comparing this with StaggeredFreeFlowConnectivityMap computeFaceToCellCenterStencil_: SIMPLE does not consider outerParallelElements, such terms are included in velocity coefficients
                updateGlobalJacobian_(A, faceGlobalI, scvf.insideScvIdx(), pvIdx, simpleMomentumBalanceSummandsVector[scvf.localFaceIdx()].pressureCoefficient);
            }
        }
    }

    template<class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDerivatives(const std::vector<std::size_t>& additionalDofDependencies,
                                   JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {

    }

    /*!
     * \brief Updates the current global Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col'. Specialization for cc methods.
     */
    template<class SubMatrix, class CCOrFacePrimaryVariables>
    static void updateGlobalJacobian_(SubMatrix& matrix,
                                      const int globalI,
                                      const int globalJ,
                                      const int pvIdx,
                                      const CCOrFacePrimaryVariables &partialDeriv)
    {
        for (int eqIdx = 0; eqIdx < partialDeriv.size(); eqIdx++)
        {
            // A[i][col][eqIdx][pvIdx] is the rate of change of
            // the residual of equation 'eqIdx' at dof 'i'
            // depending on the primary variable 'pvIdx' at dof
            // 'col'.

            assert(pvIdx >= 0);
            assert(eqIdx < matrix[globalI][globalJ].size());
            assert(pvIdx < matrix[globalI][globalJ][eqIdx].size());
            matrix[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
        }
    }
};

} // end namespace Dumux

#endif
