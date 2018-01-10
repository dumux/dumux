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
#ifndef DUMUX_FV_LOCAL_ASSEMBLER_BASE_HH
#define DUMUX_FV_LOCAL_ASSEMBLER_BASE_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericdifferentiation.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<class TypeTag, class Assembler, class Implementation>
class FVLocalAssemblerBase
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
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

    explicit FVLocalAssemblerBase(const Assembler& assembler,
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
    , elementIsGhost_((element.partitionType() == Dune::GhostEntity))
    {}

    auto evalLocalResidual() const
    {
        if(this->assembler().isStationaryProblem() && !isImplicit())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if(elementIsGhost())
        {
            using ResdiualType = decltype(evalLocalResidual(std::declval<ElementVolumeVariables>()));
            return ResdiualType(0.0);
        }

        return isImplicit() ? evalLocalResidual(curElemVolVars())
                            : evalLocalResidual(prevElemVolVars());
    }


    auto evalLocalResidual(const ElementVolumeVariables& elemVolVars) const
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

    auto evalLocalFluxSourceResidual() const
    {
        return isImplicit() ? evalLocalFluxSourceResidual(curElemVolVars())
                            : evalLocalFluxSourceResidual(prevElemVolVars());
     }

    auto evalLocalFluxSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return localResidual_.evalFluxSource(element_, fvGeometry_, elemVolVars, elemFluxVarsCache_, elemBcTypes_);
    }

    auto evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(element_, fvGeometry_, prevElemVolVars_, curElemVolVars_);
    }

    auto evalFluxResidual(const Element& neigbor,
                                         const SubControlVolumeFace& scvf) const
    {
        return localResidual_.evalFlux(problem(), neigbor, fvGeometry_, curElemVolVars_, elemFluxVarsCache_, scvf);
    }

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

        if(isImplicit())
        {
            curElemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
            if (!this->assembler().isStationaryProblem())
                prevElemVolVars.bindElement(element, fvGeometry, this->assembler().prevSol());
        }
        else
        {
            curElemVolVars.bindElement(element, fvGeometry, curSol);
            prevElemVolVars.bind(element, fvGeometry, prevSol);
            elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
        }
    }

    const Problem& problem() const
    { return assembler_.problem(); }

    const Assembler& assembler() const
    { return assembler_; }

    const Element& element() const
    { return element_; }

    bool elementIsGhost() const
    { return elementIsGhost_; }

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

    LocalResidual& localResidual()
    { return localResidual_; }

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

    const LocalResidual& localResidual() const
    { return localResidual_; }

    template<class T = TypeTag, typename std::enable_if_t<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag, typename std::enable_if_t<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

    constexpr bool isImplicit() const
    { return Implementation::isImplicit(); }
    // { return Implementation::IsImplicit::value; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution

    FVElementGeometry fvGeometry_;
    ElementVolumeVariables curElemVolVars_;
    ElementVolumeVariables prevElemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_; //!< whether the element's partitionType is ghost
};


} // end namespace Dumux

#endif
