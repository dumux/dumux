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
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/diffmethod.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<class TypeTag, class Assembler, class Implementation>
class CCLocalAssemblerBase
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Element = typename GridView::template Codim<0>::Entity;

public:

    explicit CCLocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
    , fvGeometry_(localView(assembler.fvGridGeometry()))
    , curElemVolVars_(localView(assembler.gridVariables().curGridVolVars()))
    , prevElemVolVars_(localView(assembler.gridVariables().prevGridVolVars()))
    , elemFluxVarsCache_(localView(assembler.gridVariables().gridFluxVarsCache()))
    , localResidual_(assembler.localResidual())
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables)
    {
        asImp_().bindLocalViews();
        const auto globalI = assembler().fvGridGeometry().elementMapper().index(this->element());
        res[globalI] = asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        asImp_().bindLocalViews();
        asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SolutionVector& res)
    {
        asImp_().bindLocalViews();
        const auto globalI = assembler().fvGridGeometry().elementMapper().index(this->element());
        res[globalI] = asImp_().assembleResidualImpl(); // forward to the internal implementation
    }

    /*!
     * \brief Computes the epsilon used for numeric differentiation
     *        for a given value of a primary variable.
     *
     * \param priVar The value of the primary variable
     */
    static Scalar numericEpsilon(const Scalar priVar)
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
    }

    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

    LocalResidualValues evalLocalResidual(const ElementVolumeVariables& elemVolVars) const
    {
        if (!assembler().isStationaryProblem())
        {
            auto residual = evalLocalFluxSourceResidual(elemVolVars);
            residual += evalLocalStorageResidual();
            return residual;
        }
        else
            return evalLocalFluxSourceResidual(elemVolVars);
    }

    LocalResidualValues evalLocalFluxSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return localResidual_.evalFluxSource(element_, fvGeometry_, elemVolVars, elemFluxVarsCache_, elemBcTypes_)[0];
    }

    LocalResidualValues evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(element_, fvGeometry_, prevElemVolVars_, curElemVolVars_)[0];
    }

    LocalResidualValues evalFluxResidual(const Element& neigbor,
                                         const SubControlVolumeFace& scvf) const
    {
        return localResidual_.evalFlux(problem(), neigbor, fvGeometry_, curElemVolVars_, elemFluxVarsCache_, scvf);
    }

    const Problem& problem() const
    { return assembler_.problem(); }

    const Assembler& assembler() const
    { return assembler_; }

    const Element& element() const
    { return element_; }

    bool elementIsGhost() const
    { return (element_.partitionType() == Dune::GhostEntity); }

    const SolutionVector& curSol() const
    { return curSol_; }

    FVElementGeometry& fvGeometry()
    { return fvGeometry_; }

    ElementVolumeVariables& curElemVolVars()
    { return curElemVolVars_; }

    ElementVolumeVariables& prevElemVolVars()
    { return prevElemVolVars_; }

    ElementFluxVariablesCache& elemFluxVarsCache()
    { return elemFluxVarsCache_; }

    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    const ElementVolumeVariables& curElemVolVars() const
    { return curElemVolVars_; }

    const ElementVolumeVariables& prevElemVolVars() const
    { return prevElemVolVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution

    FVElementGeometry fvGeometry_;
    ElementVolumeVariables curElemVolVars_;
    ElementVolumeVariables prevElemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all implicit local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<class TypeTag, class Assembler, class Implementation>
class CCLocalAssemblerImplicitBase : public CCLocalAssemblerBase<TypeTag, Assembler, Implementation>
{
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, Implementation>;
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol();
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        fvGeometry.bind(element);
        curElemVolVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
        if (!this->assembler().isStationaryProblem())
            this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol());
    }

    using ParentType::evalLocalResidual;
    LocalResidualValues evalLocalResidual() const
    { return this->evalLocalResidual(this->curElemVolVars()); }

    using ParentType::evalLocalFluxSourceResidual;
    LocalResidualValues evalLocalFluxSourceResidual() const
    { return this->evalLocalFluxSourceResidual(this->curElemVolVars()); }

    /*!
     * \brief Computes the residual
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleResidualImpl()
    { return this->elementIsGhost() ? LocalResidualValues(0.0) : evalLocalResidual(); }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all explicit local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<class TypeTag, class Assembler, class Implementation>
class CCLocalAssemblerExplicitBase : public CCLocalAssemblerBase<TypeTag, Assembler, Implementation>
{
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, Implementation>;
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol();
        const auto& prevSol = this->assembler().prevSol();
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& prevElemVolVars = this->prevElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        fvGeometry.bind(element);
        curElemVolVars.bindElement(element, fvGeometry, curSol);
        prevElemVolVars.bind(element, fvGeometry, prevSol);
        elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
    }

    using ParentType::evalLocalResidual;
    LocalResidualValues evalLocalResidual() const
    { return this->evalLocalResidual(this->prevElemVolVars()); }

    using ParentType::evalLocalFluxSourceResidual;
    LocalResidualValues evalLocalFluxSourceResidual() const
    { return this->evalLocalFluxSourceResidual(this->prevElemVolVars()); }

    /*!
     * \brief Computes the residual
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleResidualImpl()
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        return this->elementIsGhost() ? LocalResidualValues(0.0) : evalLocalResidual();
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class CCLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public CCLocalAssemblerImplicitBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = CCLocalAssemblerImplicitBase<TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        // assemble the undeflected residual
        const auto residual = this->assembleResidualImpl();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get the numeric differentiation type (-1 : backward differences, 0: central differences, 1: forward differences)
        static const std::string group = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        static const int numericDifferenceMethod = getParamFromGroup<int>(group, "Implicit.NumericDifferenceMethod");

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->assembler().fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        std::vector<Element> neighborElements;
        neighborElements.reserve(numNeighbors);

        // all undeflected fluxes influences by this dof's primary variables
        Dune::BlockVector<LocalResidualValues> origFlux(numNeighbors); origFlux = 0.0;

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        unsigned int j = 0;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements.emplace_back(fvGridGeometry.element(dataJ.globalJ));
            for (const auto scvfIdx : dataJ.scvfsJ)
                origFlux[j++] += this->evalFluxResidual(neighborElements.back(), fvGeometry.scvf(scvfIdx));
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        ElementSolutionVector elemSol(origPriVars);

        // derivatives in the neighbors with repect to the current elements
        Dune::BlockVector<LocalResidualValues> neighborDeriv(numNeighbors);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            // reset derivatives of element dof with respect to itself
            // as well as neighbor derivatives
            LocalResidualValues partialDeriv(0.0);
            neighborDeriv = 0.0;

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with repect to the initial solution this
            // results in a delta of zero for ghosts
            if (this->elementIsGhost()) partialDeriv[pvIdx] = 1.0;

            Scalar eps = ParentType::numericEpsilon(curVolVars.priVar(pvIdx));
            Scalar delta = 0;

            // we are using forward or central differences, i.e. we need to calculate f(x + \epsilon)
            if (numericDifferenceMethod >= 0)
            {
                // deflect primary variables
                elemSol[0][pvIdx] += eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                if (!this->elementIsGhost()) partialDeriv = this->evalLocalResidual();

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        neighborDeriv[k] += this->evalFluxResidual(neighborElements[k], fvGeometry.scvf(scvfIdx));
            }

            // we are using backward differences, i.e. we don't need to calculate f(x + \epsilon)
            // we can recycle the (already calculated) residual f(x) for non-ghosts
            else
            {
                if (!this->elementIsGhost()) partialDeriv = residual;
                neighborDeriv = origFlux;
            }

            // we are using backward or central differences, i.e. we need to calculate f(x - \epsilon)
            if (numericDifferenceMethod <= 0)
            {
                // deflect the primary variables
                elemSol[0][pvIdx] -= delta + eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables and subtract it
                if (!this->elementIsGhost()) partialDeriv -= this->evalLocalResidual();

                // calculate the fluxes into element with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        neighborDeriv[k] -= this->evalFluxResidual(neighborElements[k], fvGeometry.scvf(scvfIdx));
            }

            // we are using forward differences, i.e. we don't need to calculate f(x - \epsilon)
            // we can recycle the (already calculated) residual f(x) for non-ghosts
            else
            {
                if (!this->elementIsGhost()) partialDeriv -= residual;
                neighborDeriv -= origFlux;
            }

            // divide difference in residuals by the magnitude of the
            // deflections between the two function evaluation
            if (!this->elementIsGhost()) partialDeriv /= delta;
            neighborDeriv /= delta;

            // add the current partial derivatives to the global jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                // the diagonal entries
                A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];

                // off-diagonal entries
                j = 0;
                for (const auto& dataJ : connectivityMap[globalI])
                    A[dataJ.globalJ][globalI][eqIdx][pvIdx] += neighborDeriv[j++][eqIdx];
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // return the original residual
        return residual;
    }
};


/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and explicits time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public CCLocalAssemblerExplicitBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false> >
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>;
    using ParentType = CCLocalAssemblerExplicitBase<TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // assemble the undeflected residual
        auto residual = this->elementIsGhost() ? LocalResidualValues(0.0) : this->evalLocalFluxSourceResidual();
        const auto storageResidual = this->elementIsGhost() ? LocalResidualValues(0.0) : this->evalLocalStorageResidual();
        residual += storageResidual;

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get the numeric differentiation type (-1 : backward differences, 0: central differences, 1: forward differences)
        static const std::string group = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        static const int numericDifferenceMethod = getParamFromGroup<int>(group, "Implicit.NumericDifferenceMethod");

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->assembler().fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        ElementSolutionVector elemSol(origPriVars);

        // derivatives in the neighbors with repect to the current elements
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            // reset derivatives of element dof with respect to itself
            // as well as neighbor derivatives
            LocalResidualValues partialDeriv(0.0);
            if (this->elementIsGhost()) partialDeriv[pvIdx] = 1.0;

            Scalar eps = ParentType::numericEpsilon(curVolVars.priVar(pvIdx));
            Scalar delta = 0;

            // we are using forward or central differences, i.e. we need to calculate f(x + \epsilon)
            if (numericDifferenceMethod >= 0)
            {
                // deflect primary variables
                elemSol[0][pvIdx] += eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);

                // calculate the residual with the deflected primary variables
                if (!this->elementIsGhost()) partialDeriv = this->evalLocalStorageResidual();
            }

            // we are using backward differences, i.e. we don't need to calculate f(x + \epsilon)
            // we can recycle the (already calculated) residual f(x) for non-ghosts
            else
            {
                if (!this->elementIsGhost()) partialDeriv = storageResidual;
            }

            // we are using backward or central differences, i.e. we need to calculate f(x - \epsilon)
            if (numericDifferenceMethod <= 0)
            {
                // deflect the primary variables
                elemSol[0][pvIdx] -= delta + eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);

                // calculate the residual with the deflected primary variables and subtract it
                if (!this->elementIsGhost()) partialDeriv -= this->evalLocalStorageResidual();
            }

            // we are using forward differences, i.e. we don't need to calculate f(x - \epsilon)
            // we can recycle the (already calculated) residual f(x) for non-ghosts
            else
            {
                if (!this->elementIsGhost()) partialDeriv -= storageResidual;
            }

            // divide difference in residuals by the magnitude of the
            // deflections between the two function evaluation
            if (!this->elementIsGhost()) partialDeriv /= delta;

            // add the current partial derivatives to the global jacobian matrix
            // only diagonal entries for explicit jacobians
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];
        }

        // return the original residual
        return residual;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public CCLocalAssemblerExplicitBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true> >
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>;
    using ParentType = CCLocalAssemblerImplicitBase<TypeTag, Assembler, ThisType>;
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        // assemble the undeflected residual
        const auto residual = this->assembleResidualImpl();

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // if the problem is instationary, add derivative of storage term
        if (!this->assembler().isStationaryProblem())
            this->localResidual().addStorageDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add source term derivatives
        this->localResidual().addSourceDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add flux derivatives for each scvf
        for (const auto& scvf : scvfs(fvGeometry))
        {
            // inner faces
            if (!scvf.boundary())
                this->localResidual().addFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

            // boundary faces
            else
            {
                const auto& bcTypes = problem.boundaryTypes(element, scvf);

                // add Dirichlet boundary flux derivatives
                if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                    this->localResidual().addCCDirichletFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                // add Robin ("solution dependent Neumann") boundary flux derivatives
                else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
                    this->localResidual().addRobinFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                else
                    DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
            }
        }

        // return element residual
        return residual;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false>
: public CCLocalAssemblerExplicitBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false> >
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>;
    using ParentType = CCLocalAssemblerImplicitBase<TypeTag, Assembler, ThisType>;
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        // assemble the undeflected residual
        const auto residual = this->assembleResidualImpl();

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // add hand-code derivative of storage term
        this->localResidual().addStorageDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // return the original residual
        return residual;
    }
};

} // end namespace Dumux

#endif
