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
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \ingroup MultiDomain
 * \brief A multidomain assembler for Jacobian and residual contribution per element (staggered method)
 */
#ifndef DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

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

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using SolutionVector = typename Assembler::SolutionVector;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using ElementFaceVariables = typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    static constexpr auto domainId = typename Dune::index_constant<id>();
    static constexpr auto cellCenterId = GridGeometry::cellCenterIdx();
    static constexpr auto faceId = GridGeometry::faceIdx();

    static constexpr auto numEqCellCenter = CellCenterResidualValue::dimension;
    static constexpr auto faceOffset = numEqCellCenter;

    using ParentType::ParentType;

    explicit SubDomainStaggeredLocalAssemblerBase(const Assembler& assembler,
                                                  const Element& element,
                                                  const SolutionVector& curSol,
                                                  CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.gridGeometry(domainId)),
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
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        assembleJacobianAndResidualImpl_(domainId, jacRow, res, gridVariables);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class SubSol>
    void assembleResidual(SubSol& res)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        if constexpr (domainId == cellCenterId)
        {
            const auto cellCenterGlobalI = problem().gridGeometry().elementMapper().index(this->element());
            res[cellCenterGlobalI] = this->asImp_().assembleCellCenterResidualImpl();
        }
        else
        {
            for (auto&& scvf : scvfs(this->fvGeometry()))
                res[scvf.dofIndex()] +=  this->asImp_().assembleFaceResidualImpl(scvf);
        }
    }

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    CellCenterResidualValue evalLocalResidualForCellCenter() const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            return CellCenterResidualValue(0.0);
        }

        return isImplicit ? evalLocalResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current cell center.
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    CellCenterResidualValue evalLocalResidualForCellCenter(const ElementVolumeVariables& elemVolVars,
                                                           const ElementFaceVariables& elemFaceVars) const
    {
        auto residual = evalLocalFluxAndSourceResidualForCellCenter(elemVolVars, elemFaceVars);

        if (!this->assembler().isStationaryProblem())
            residual += evalLocalStorageResidualForCellCenter();

        // handle cells with a fixed Dirichlet value
        const auto cellCenterGlobalI = problem().gridGeometry().elementMapper().index(this->element());
        const auto& scvI = this->fvGeometry().scv(cellCenterGlobalI);
        for (int pvIdx = 0; pvIdx < numEqCellCenter; ++pvIdx)
        {
            static constexpr auto offset = numEq - numEqCellCenter;
            if (this->problem().isDirichletCell(this->element(), this->fvGeometry(), scvI, pvIdx + offset))
                residual[pvIdx] = this->curSol()[cellCenterId][cellCenterGlobalI][pvIdx] - this->problem().dirichlet(this->element(), scvI)[pvIdx + offset];
        }

        return residual;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    CellCenterResidualValue evalLocalFluxAndSourceResidualForCellCenter() const
    {
        return isImplicit ? evalLocalFluxAndSourceResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalFluxAndSourceResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars());
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    CellCenterResidualValue evalLocalFluxAndSourceResidualForCellCenter(const ElementVolumeVariables& elemVolVars, const ElementFaceVariables& elemFaceVars) const
    {
        return this->localResidual().evalFluxAndSourceForCellCenter(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache());
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    CellCenterResidualValue evalLocalStorageResidualForCellCenter() const
    {
        return this->localResidual().evalStorageForCellCenter(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars());
    }

    /*!
     * \brief Convenience function to evaluate the  local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf) const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            return FaceResidualValue(0.0);
        }

        return isImplicit ? evalLocalResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars) const
    {
        auto residual = evalLocalFluxAndSourceResidualForFace(scvf, elemVolVars, elemFaceVars);

        if (!this->assembler().isStationaryProblem())
            residual += evalLocalStorageResidualForFace(scvf);

        this->localResidual().evalDirichletBoundariesForFace(residual, this->problem(), this->element(),
                                                             this->fvGeometry(), scvf, elemVolVars, elemFaceVars,
                                                             this->elemBcTypes(), this->elemFluxVarsCache());

        return residual;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf) const
    {
        return isImplicit ? evalLocalFluxAndSourceResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalFluxAndSourceResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    FaceResidualValue evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf,
                                                            const ElementVolumeVariables& elemVolVars,
                                                            const ElementFaceVariables& elemFaceVars) const
    {
        return this->localResidual().evalFluxAndSourceForFace(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalStorageResidualForFace(const SubControlVolumeFace& scvf) const
    {
        return this->localResidual().evalStorageForFace(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars(), this->prevElemFaceVars(), this->curElemFaceVars(), scvf);
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

    //! Assembles the residuals and derivatives for the cell center dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    auto assembleJacobianAndResidualImpl_(Dune::index_constant<cellCenterId>, JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        auto& gridVariablesI = *std::get<domainId>(gridVariables);
        const auto cellCenterGlobalI = problem().gridGeometry().elementMapper().index(this->element());
        const auto residual = this->asImp_().assembleCellCenterJacobianAndResidualImpl(jacRow[domainId], gridVariablesI);
        res[cellCenterGlobalI] = residual;


        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<JacobianMatrixRow::size(), domainId>{};
        forEach(otherDomainIds, [&](auto&& domainJ)
        {
            this->asImp_().assembleJacobianCellCenterCoupling(domainJ, jacRow[domainJ], residual, gridVariablesI);
        });

        // handle cells with a fixed Dirichlet value
        incorporateDirichletCells_(jacRow);
    }

    //! Assembles the residuals and derivatives for the face dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleJacobianAndResidualImpl_(Dune::index_constant<faceId>, JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        auto& gridVariablesI = *std::get<domainId>(gridVariables);
        const auto residual = this->asImp_().assembleFaceJacobianAndResidualImpl(jacRow[domainId], gridVariablesI);

        for(auto&& scvf : scvfs(this->fvGeometry()))
            res[scvf.dofIndex()] += residual[scvf.localFaceIdx()];

        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<JacobianMatrixRow::size(), domainId>{};
        forEach(otherDomainIds, [&](auto&& domainJ)
        {
            this->asImp_().assembleJacobianFaceCoupling(domainJ, jacRow[domainJ], residual, gridVariablesI);
        });
    }

    //! If specified in the problem, a fixed Dirichlet value can be assigned to cell centered unknows such as pressure
    template<class JacobianMatrixRow>
    void incorporateDirichletCells_(JacobianMatrixRow& jacRow)
    {
        const auto cellCenterGlobalI = problem().gridGeometry().elementMapper().index(this->element());

        // overwrite the partial derivative with zero in case a fixed Dirichlet BC is used
        static constexpr auto offset = numEq - numEqCellCenter;
        for (int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
        {
            if (problem().isDirichletCell(this->element(), this->fvGeometry(), this->fvGeometry().scv(cellCenterGlobalI), eqIdx + offset))
            {
                using namespace Dune::Hybrid;
                forEach(integralRange(Dune::Hybrid::size(jacRow)), [&, domainId = domainId](auto&& i)
                {
                    auto& ccRowI = jacRow[i][cellCenterGlobalI];
                    for (auto col = ccRowI.begin(); col != ccRowI.end(); ++col)
                    {
                        ccRowI[col.index()][eqIdx] = 0.0;
                        // set the diagonal entry to 1.0
                        if ((i == domainId) && (col.index() == cellCenterGlobalI))
                            ccRowI[col.index()][eqIdx][eqIdx] = 1.0;
                    }
                });
            }
        }
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
 *
 * The assembly of the cellCenterResidual is done element-wise, the assembly of the face residual is done half-element-wise.
 *
 * A sketch of what this means can be found in the following image:
 *
 * \image html staggered_halfelementwise.png
 *
 * Half-element wise assembly means, that integrals are
split into contributions from the left and right part of the staggered control volume. For an example term \f$\int_{\Omega}\varrho u\text{d}\Omega\f$ this reads
\f$\int_{\Omega}\varrho u\text{d}\Omega = \frac{1}{2}\Omega_\text{left}\varrho_\text{left}u+\frac{1}{2}\Omega_\text{right}\varrho_\text{right}u\f$.

During assembly, \f$\frac{1}{2}\Omega_\text{left}\varrho_\text{left}u\f$ is added to the residual of the
staggered control volume \f$\Omega\f$, when the loops reach the scvf within element \f$\Omega_\text{left}\f$. \f$\frac{1}{2}\Omega_\text{right}\varrho_\text{right}u\f$ is added to the residual of \f$\Omega\f$, when the loops reach the scvf within in element \f$\Omega_\text{right}\f$.
Other terms are split analogously.
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridFaceVariables = GetPropType<TypeTag, Properties::GridFaceVariables>;
    using ElementFaceVariables = typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;
    using FaceVariables = typename ElementFaceVariables::FaceVariables;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr bool enableGridFluxVarsCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    static constexpr int maxNeighbors = 4*(2*ModelTraits::dim());
    static constexpr auto domainI = Dune::index_constant<id>();
    static constexpr auto cellCenterId = GridGeometry::cellCenterIdx();
    static constexpr auto faceId = GridGeometry::faceIdx();

    static constexpr auto numEq = ModelTraits::numEq();
    static constexpr auto numEqCellCenter = CellCenterPrimaryVariables::dimension;
    static constexpr auto numEqFace = FacePrimaryVariables::dimension;

public:
    using ParentType::ParentType;

    CellCenterResidualValue assembleCellCenterResidualImpl()
    {
        return this->evalLocalResidualForCellCenter();
    }

    FaceResidualValue assembleFaceResidualImpl(const SubControlVolumeFace& scvf)
    {
        return this->evalLocalResidualForFace(scvf);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    CellCenterResidualValue assembleCellCenterJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        assert(domainI == cellCenterId);

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        const auto& gridGeometry = this->problem().gridGeometry();
        const auto& curSol = this->curSol()[domainI];

        const auto cellCenterGlobalI = gridGeometry.elementMapper().index(element);
        const auto origResidual = this->evalLocalResidualForCellCenter();

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all cell center residuals in the element w.r.t. to other cell center dofs. //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        // lambda to evaluate the derivatives for cell center dofs with respect to neighbor cells
        auto evaluateCellCenterDerivatives = [&](const std::size_t globalJ)
        {
            // get the volVars of the element with respect to which we are going to build the derivative
            auto&& scvJ = fvGeometry.scv(globalJ);
            const auto elementJ = fvGeometry.gridGeometry().element(globalJ);
            auto& curVolVars =  this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
            const auto origVolVars(curVolVars);

            for (int pvIdx = 0; pvIdx < numEqCellCenter; ++pvIdx)
            {
                CellCenterPrimaryVariables cellCenterPriVars = curSol[globalJ];
                using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

                constexpr auto offset = numEq - numEqCellCenter;

                auto evalResidual = [&](Scalar priVar)
                {
                    // update the volume variables
                    priVars[pvIdx + offset] = priVar;
                    auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                    curVolVars.update(elemSol, this->problem(), elementJ, scvJ);

                    // update the coupling context
                    cellCenterPriVars[pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalJ, cellCenterPriVars, pvIdx);

                    // compute element residual
                    return this->evalLocalResidualForCellCenter();
                };

                // create the vector storing the partial derivatives
                CellCenterResidualValue partialDeriv(0.0);

                // derive the residuals numerically
                const auto& paramGroup = this->problem().paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto eps = this->couplingManager().numericEpsilon(domainI, paramGroup);
                NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, origResidual,
                                                          eps(priVars[pvIdx + offset], pvIdx), numDiffMethod);

                // update the global jacobian matrix with the current partial derivatives
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the original volVars
                curVolVars = origVolVars;

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalJ, curSol[globalJ], pvIdx);
            }
        };

        // get the list of cell center dofs that have an influence on the cell center resdiual of the current element
        const auto& connectivityMap = gridGeometry.connectivityMap();

        // evaluate derivatives w.r.t. own dof
        evaluateCellCenterDerivatives(cellCenterGlobalI);

        // evaluate derivatives w.r.t. all other related cell center dofs
        for (const auto& globalJ : connectivityMap(cellCenterId, cellCenterId, cellCenterGlobalI))
             evaluateCellCenterDerivatives(globalJ);

        return origResidual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    auto assembleFaceJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        assert(domainI == faceId);

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->problem().gridGeometry();
        const auto& curSol = this->curSol()[domainI];

        using FaceSolutionVector = GetPropType<TypeTag, Properties::FaceSolutionVector>; // TODO: use reserved vector
        FaceSolutionVector origResiduals;
        origResiduals.resize(fvGeometry.numScvf());
        origResiduals = 0.0;

        // treat the local residua of the face dofs:
        for (auto&& scvf : scvfs(fvGeometry))
            origResiduals[scvf.localFaceIdx()] = this->evalLocalResidualForFace(scvf);

        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all face residuals in the element w.r.t. to other face dofs.         //
        // Note that we do an element-wise assembly, therefore this is only the contribution of the      //
        // current element while the contribution of the element on the opposite side of the scvf will   //
        // be added separately.                                                                          //
        ///////////////////////////////////////////////////////////////////////////////////////////////////

        for (auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            using FaceSolution = GetPropType<TypeTag, Properties::StaggeredFaceSolution>;
            const auto origFaceSolution = FaceSolution(scvf, curSol, gridGeometry);

            // Lambda to evaluate the derivatives for faces
            auto evaluateFaceDerivatives = [&](const std::size_t globalJ)
            {
                // get the faceVars of the face with respect to which we are going to build the derivative
                auto& faceVars = getFaceVarAccess_(gridVariables.curGridFaceVars(), this->curElemFaceVars(), scvf);
                const auto origFaceVars = faceVars;

                for (int pvIdx = 0; pvIdx < numEqFace; ++pvIdx)
                {
                    auto faceSolution = origFaceSolution;

                    auto evalResidual = [&](Scalar priVar)
                    {
                        // update the volume variables
                        faceSolution[globalJ][pvIdx] = priVar;
                        faceVars.update(faceSolution, problem, element, fvGeometry, scvf);

                        // update the coupling context
                        this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalJ, faceSolution[globalJ], pvIdx);

                        // compute face residual
                        return this->evalLocalResidualForFace(scvf);
                    };

                    // derive the residuals numerically
                    FaceResidualValue partialDeriv(0.0);
                    const auto& paramGroup = problem.paramGroup();
                    static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                    static const auto eps = this->couplingManager().numericEpsilon(domainI, paramGroup);
                    NumericDifferentiation::partialDerivative(evalResidual, faceSolution[globalJ][pvIdx], partialDeriv, origResiduals[scvf.localFaceIdx()],
                                                              eps(faceSolution[globalJ][pvIdx], pvIdx), numDiffMethod);

                    // update the global jacobian matrix with the current partial derivatives
                    updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original faceVars
                    faceVars = origFaceVars;

                    // restore the undeflected state of the coupling context
                    this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalJ, origFaceSolution[globalJ], pvIdx);
                }
            };

            // evaluate derivatives w.r.t. own dof
            evaluateFaceDerivatives(scvf.dofIndex());

            // get the list of face dofs that have an influence on the resdiual of the current face
            const auto& connectivityMap = gridGeometry.connectivityMap();

            // evaluate derivatives w.r.t. all other related face dofs
            for (const auto& globalJ : connectivityMap(faceId, faceId, scvf.index()))
               evaluateFaceDerivatives(globalJ);
        }

        return origResiduals;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class GridVariables>
    void assembleJacobianCellCenterCoupling(Dune::index_constant<faceId> domainJ, JacobianBlock& A,
                                            const CellCenterResidualValue& origResidual, GridVariables& gridVariables)
    {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Calculate derivatives of all cell center residuals in the element w.r.t. to all coupled faces dofs //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->problem().gridGeometry();
        const auto& curSol = this->curSol()[domainJ];
        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto cellCenterGlobalI = gridGeometry.elementMapper().index(element);

        for (const auto& scvfJ : scvfs(fvGeometry))
        {
            const auto globalJ = scvfJ.dofIndex();

            // get the faceVars of the face with respect to which we are going to build the derivative
            auto& faceVars = getFaceVarAccess_(gridVariables.curGridFaceVars(), this->curElemFaceVars(), scvfJ);
            const auto origFaceVars(faceVars);

            for (int pvIdx = 0; pvIdx < numEqFace; ++pvIdx)
            {
                auto facePriVars = curSol[globalJ];

                auto evalResidual = [&](Scalar priVar)
                {
                    // update the face variables
                    facePriVars[pvIdx] = priVar;
                    faceVars.updateOwnFaceOnly(facePriVars);

                    // update the coupling context
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, facePriVars, pvIdx);

                    // compute element residual
                    return this->evalLocalResidualForCellCenter();
                };

                // create the vector storing the partial derivatives
                CellCenterResidualValue partialDeriv(0.0);

                // derive the residuals numerically
                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);
                NumericDifferentiation::partialDerivative(evalResidual, facePriVars[pvIdx], partialDeriv, origResidual,
                                                          epsCoupl(facePriVars[pvIdx], pvIdx), numDiffMethod);

                // update the global jacobian matrix with the current partial derivatives
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the original faceVars
                faceVars = origFaceVars;

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, curSol[globalJ], pvIdx);
            }
        }
    }

    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCellCenterCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                            const CellCenterResidualValue& res, GridVariables& gridVariables)
    {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Calculate derivatives of all cell center residuals in the element w.r.t. all dofs in the coupling stencil. //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();

        // get stencil informations
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);

        if (stencil.empty())
            return;

        for (const auto globalJ : stencil)
        {
            const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
            const auto& curSol = this->curSol()[domainJ];
            const auto origPriVarsJ = curSol[globalJ];

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    auto deflectedPriVarsJ = origPriVarsJ;
                    deflectedPriVarsJ[pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, deflectedPriVarsJ, pvIdx);
                    return this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
                };

                // create the vector storing the partial derivatives
                CellCenterResidualValue partialDeriv(0.0);

                // derive the residuals numerically
                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);
                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDeriv, origResidual,
                                                          epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                const auto cellCenterGlobalI = this->problem().gridGeometry().elementMapper().index(element);
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, origPriVarsJ, pvIdx);
            }
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class ElementResidualVector, class GridVariables>
    void assembleJacobianFaceCoupling(Dune::index_constant<cellCenterId> domainJ, JacobianBlock& A,
                                      const ElementResidualVector& origResiduals, GridVariables& gridVariables)
    {
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Calculate derivatives of all face residuals in the element w.r.t. all coupled cell center dofs //
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->problem().gridGeometry();
        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto& curSol = this->curSol()[domainJ];

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for (const auto& globalJ : connectivityMap(faceId, cellCenterId, scvf.index()))
            {
                // get the volVars of the element with respect to which we are going to build the derivative
                auto&& scvJ = fvGeometry.scv(globalJ);
                const auto elementJ = fvGeometry.gridGeometry().element(globalJ);
                auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), this->curElemVolVars(), scvJ);
                const auto origVolVars(curVolVars);
                const auto origCellCenterPriVars = curSol[globalJ];

                for (int pvIdx = 0; pvIdx < numEqCellCenter; ++pvIdx)
                {
                    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
                    PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(origCellCenterPriVars);

                    constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;

                    auto evalResidual = [&](Scalar priVar)
                    {
                        // update the volume variables
                        priVars[pvIdx + offset] = priVar;
                        auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                        curVolVars.update(elemSol, problem, elementJ, scvJ);

                        // update the coupling context
                        auto deflectedCellCenterPriVars = origCellCenterPriVars;
                        deflectedCellCenterPriVars[pvIdx] = priVar;
                        this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, deflectedCellCenterPriVars, pvIdx);

                        // compute face residual
                        return this->evalLocalResidualForFace(scvf);
                    };

                    // derive the residuals numerically
                    FaceResidualValue partialDeriv(0.0);
                    const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                    static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                    static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);
                    NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, origResiduals[scvf.localFaceIdx()],
                                                              epsCoupl(priVars[pvIdx + offset], pvIdx), numDiffMethod);

                    // update the global jacobian matrix with the current partial derivatives
                    updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original volVars
                    curVolVars = origVolVars;

                    // restore the undeflected state of the coupling context
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, origCellCenterPriVars, pvIdx);
                }
            }
        }
    }

    template<std::size_t otherId, class JacobianBlock, class ElementResidualVector, class GridVariables>
    void assembleJacobianFaceCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                      const ElementResidualVector& res, GridVariables& gridVariables)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Calculate derivatives of all face residuals in the element w.r.t. all dofs in the coupling stencil. //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol()[domainJ];

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // get stencil informations
            const auto& stencil = this->couplingManager().couplingStencil(domainI, scvf, domainJ);

            if (stencil.empty())
                continue;

            // build derivatives with for face dofs w.r.t. cell center dofs
            for (const auto& globalJ : stencil)
            {
                const auto origPriVarsJ = curSol[globalJ];
                const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, scvf, *this, domainJ, globalJ);

                for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
                {
                    auto evalCouplingResidual = [&](Scalar priVar)
                    {
                        auto deflectedPriVars = origPriVarsJ;
                        deflectedPriVars[pvIdx] = priVar;
                        this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, deflectedPriVars, pvIdx);
                        return this->couplingManager().evalCouplingResidual(domainI, scvf, *this, domainJ, globalJ);
                    };

                    // derive the residuals numerically
                    FaceResidualValue partialDeriv(0.0);
                    const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                    static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                    static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);
                    NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDeriv, origResidual,
                                                              epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                    // update the global stiffness matrix with the current partial derivatives
                    updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the undeflected state of the coupling context
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, origPriVarsJ, pvIdx);
                }
            }
        }
    }

    template<class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDerivatives(const std::vector<std::size_t>& additionalDofDependencies,
                                   JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    { }

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
                                      const CCOrFacePrimaryVariables& partialDeriv)
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

private:

    FaceVariables& getFaceVarAccess_(GridFaceVariables& gridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    {
        if constexpr (getPropValue<TypeTag, Properties::EnableGridFaceVariablesCache>())
            return gridFaceVariables.faceVars(scvf.index());
        else
            return elemFaceVars[scvf];
    }
};

} // end namespace Dumux

#endif
