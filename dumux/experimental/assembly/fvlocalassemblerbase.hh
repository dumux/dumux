// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Assembly
 * \copydoc Dumux::FVLocalAssemblerBase
 */
#ifndef DUMUX_EXPERIMENTAL_FV_LOCAL_ASSEMBLER_BASE_HH
#define DUMUX_EXPERIMENTAL_FV_LOCAL_ASSEMBLER_BASE_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/assembly/diffmethod.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \brief A base class for all local assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The assembler implementation
 */
template<class TypeTag, class Assembler, class Implementation>
class FVLocalAssemblerBase
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = typename Assembler::SolutionVector;
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
    using LocalResidual = std::decay_t<decltype(std::declval<Assembler>().localResidual())>;
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
                           localView(assembler.gridVariables().gridFluxVarsCache()),
                           assembler.localResidual(),
                           element.partitionType() == Dune::GhostEntity,
                           assembler.isImplicit())
    {}

    /*!
     * \brief The constructor. General version explicitly expecting each argument.
     */
    explicit FVLocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const LocalResidual& localResidual,
                                  const bool elementIsGhost,
                                  const bool isImplicit)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
    , fvGeometry_(fvGeometry)
    , curElemVolVars_(curElemVolVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    , localOperator_(localResidual)
    , elementIsGhost_(elementIsGhost)
    , isImplicit_(isImplicit)
    {}

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalResidual() const
    {
        if (elementIsGhost())
            return ElementResidualVector(0.0);

        return evalLocalResidual(curElemVolVars());
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
            residual += evalStorage();
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
        return evalLocalFluxAndSourceResidual(curElemVolVars());
    }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVolVars The element volume variables
     */
    ElementResidualVector evalLocalFluxAndSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return localOperator_.evalFluxAndSource(element_, fvGeometry_, elemVolVars, elemFluxVarsCache_, elemBcTypes_);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalStorage() const
    {
        return localOperator_.evalStorage(fvGeometry_, curElemVolVars_);
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
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        fvGeometry.bind(element);
        if (std::abs(this->localResidual().spatialWeight()) < 1e-6)
            curElemVolVars.bindElement(element, fvGeometry, curSol);
        else
        {
            curElemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
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
            const auto internalDirichletConstraints = asImp_().problem().hasInternalDirichletConstraint(this->element(), scvI);
            if (internalDirichletConstraints.any())
            {
                const auto dirichletValues = asImp_().problem().internalDirichlet(this->element(), scvI);
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

    //! The element flux variables cache
    ElementFluxVariablesCache& elemFluxVarsCache()
    { return elemFluxVarsCache_; }

    //! The local residual for the current element
    LocalResidual& localResidual()
    { return localOperator_; }

    //! The element's boundary types
    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    //! The finite volume geometry
    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    //! The current element volume variables
    const ElementVolumeVariables& curElemVolVars() const
    { return curElemVolVars_; }

    //! The element flux variables cache
    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

    //! The element's boundary types
    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    //! The local residual for the current element
    const LocalResidual& localResidual() const
    { return localOperator_; }

    //! If the time stepping scheme is implicit
    bool isImplicit() const
    { return isImplicit_; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template<class T = TypeTag, typename std::enable_if_t<!GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag, typename std::enable_if_t<GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

private:

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution

    FVElementGeometry fvGeometry_;
    ElementVolumeVariables curElemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localOperator_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_; //!< whether the element's partitionType is ghost
    bool isImplicit_; //!< whether the time stepping scheme is implicit
};

} // end namespace Dumux

#endif
