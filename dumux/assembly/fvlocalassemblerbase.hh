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
 * \copydoc Dumux::FVLocalAssemblerBase
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

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A base class for all local assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The assembler implementation
 * \tparam useImplicitAssembly Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool useImplicitAssembly>
class FVLocalAssemblerBase
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = typename Assembler::ResidualType;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVolumeVariables = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    /*!
     * \brief The constructor. Delegates to the general constructor.
     */
    explicit FVLocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol)
    : FVLocalAssemblerBase(assembler,
                           element,
                           curSol,
                           localView(assembler.gridGeometry()),
                           localView(assembler.gridVariables().curGridVolVars()),
                           localView(assembler.gridVariables().prevGridVolVars()),
                           localView(assembler.gridVariables().gridFluxVarsCache()),
                           assembler.localResidual(),
                           element.partitionType() == Dune::GhostEntity)
    {}

    /*!
     * \brief The constructor. General version explicitly expecting each argument.
     */
    explicit FVLocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementVolumeVariables& prevElemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const LocalResidual& localResidual,
                                  const bool elementIsGhost)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
    , fvGeometry_(fvGeometry)
    , curElemVolVars_(curElemVolVars)
    , prevElemVolVars_(prevElemVolVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    , localResidual_(localResidual)
    , elementIsGhost_(elementIsGhost)
    {}

    /*!
     * \brief Returns true if the assembler considers implicit assembly.
     */
    static constexpr bool isImplicit()
    { return useImplicitAssembly; }

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalResidual() const
    {
        if (!isImplicit())
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (elementIsGhost())
            return ElementResidualVector(0.0);

        return isImplicit() ? evalLocalResidual(curElemVolVars())
                            : evalLocalResidual(prevElemVolVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current element.
     * \param elemVolVars The element volume variables
     */
    ElementResidualVector evalLocalResidual(const ElementVolumeVariables& elemVolVars) const
    {
        if (!assembler().isStationaryProblem())
        {
            ElementResidualVector residual = evalLocalFluxAndSourceResidual(elemVolVars);
            residual += evalLocalStorageResidual();
            return residual;
        }
        else
            return evalLocalFluxAndSourceResidual(elemVolVars);
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalFluxAndSourceResidual() const
    {
        return isImplicit() ? evalLocalFluxAndSourceResidual(curElemVolVars())
                            : evalLocalFluxAndSourceResidual(prevElemVolVars());
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVolVars The element volume variables
     */
    ElementResidualVector evalLocalFluxAndSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return localResidual_.evalFluxAndSource(element_, fvGeometry_, elemVolVars, elemFluxVarsCache_, elemBcTypes_);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(element_, fvGeometry_, prevElemVolVars_, curElemVolVars_);
    }

    /*!
     * \brief Convenience function bind and prepare all relevant variables required for the
     *        evaluation of the local residual.
     */
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

        if (isImplicit())
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

    /*!
     * \brief Enforces Dirichlet constraints if enabled in the problem
     */
    template<typename ApplyFunction, class P = Problem, typename std::enable_if_t<P::enableInternalDirichletConstraints(), int> = 0>
    void enforceInternalDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet constraints strongly by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        for (const auto& scvI : scvs(this->fvGeometry()))
        {
            const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scvI);
            if (internalDirichletConstraints.any())
            {
                const auto dirichletValues = this->problem().internalDirichlet(this->element(), scvI);
                // set the Dirichlet conditions in residual and jacobian
                for (int eqIdx = 0; eqIdx < internalDirichletConstraints.size(); ++eqIdx)
                    if (internalDirichletConstraints[eqIdx])
                        applyDirichlet(scvI, dirichletValues, eqIdx, eqIdx);
            }
        }
    }

    template<typename ApplyFunction, class P = Problem, typename std::enable_if_t<!P::enableInternalDirichletConstraints(), int> = 0>
    void enforceInternalDirichletConstraints(const ApplyFunction& applyDirichlet)
    {}

    //! The problem
    const Problem& problem() const
    { return assembler_.problem(); }

    //! The assembler
    const Assembler& assembler() const
    { return assembler_; }

    //! The current element
    const Element& element() const
    { return element_; }

    //! Returns if element is a ghost entity
    bool elementIsGhost() const
    { return elementIsGhost_; }

    //! The current solution
    const SolutionVector& curSol() const
    { return curSol_; }

    //! The global finite volume geometry
    FVElementGeometry& fvGeometry()
    { return fvGeometry_; }

    //! The current element volume variables
    ElementVolumeVariables& curElemVolVars()
    { return curElemVolVars_; }

    //! The element volume variables of the provious time step
    ElementVolumeVariables& prevElemVolVars()
    { return prevElemVolVars_; }

    //! The element flux variables cache
    ElementFluxVariablesCache& elemFluxVarsCache()
    { return elemFluxVarsCache_; }

    //! The local residual for the current element
    LocalResidual& localResidual()
    { return localResidual_; }

    //! The element's boundary types
    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    //! The finite volume geometry
    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    //! The current element volume variables
    const ElementVolumeVariables& curElemVolVars() const
    { return curElemVolVars_; }

    //! The element volume variables of the provious time step
    const ElementVolumeVariables& prevElemVolVars() const
    { return prevElemVolVars_; }

    //! The element flux variables cache
    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

    //! The element's boundary types
    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    //! The local residual for the current element
    const LocalResidual& localResidual() const
    { return localResidual_; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template<class T = TypeTag, typename std::enable_if_t<!getPropValue<T, Properties::EnableGridVolumeVariablesCache>(), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag, typename std::enable_if_t<getPropValue<T, Properties::EnableGridVolumeVariablesCache>(), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

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
