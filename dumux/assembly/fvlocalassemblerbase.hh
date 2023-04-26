// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
    using ElementVariables = typename GridVariables::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;

    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
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
                                  GridVariables& gridVars)
    : FVLocalAssemblerBase(assembler,
                           element,
                           gridVars,
                           localView(gridVars.gridGeometry()),
                           localView(gridVars),
                           localView(assembler.prevGridVariables()),
                           assembler.localResidual(),
                           element.partitionType() == Dune::GhostEntity)
    {}

    /*!
     * \brief The constructor. General version explicitly expecting each argument.
     * \todo do we really need this beast in the public interface???
     */
    explicit FVLocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  GridVariables& gridVars,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVariables& curElemVars,
                                  const ElementVariables& prevElemVars,
                                  const LocalResidual& localResidual,
                                  const bool elementIsGhost)
    : assembler_(assembler)
    , element_(element)
    , gridVars_(gridVars)
    , fvGeometry_(fvGeometry)
    , curElemVars_(curElemVars)
    , prevElemVars_(prevElemVars)
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

        return isImplicit() ? evalLocalResidual(curElemVars_)
                            : evalLocalResidual(prevElemVars_);
    }

    /*!
     * \brief Evaluates the complete local residual for the current element.
     * \param elemVolVars The element volume variables
     */
    ElementResidualVector evalLocalResidual(const ElementVariables& elemVars) const
    {
        if (!assembler().isStationaryProblem())
        {
            ElementResidualVector residual = evalLocalFluxAndSourceResidual(elemVars);
            residual += evalLocalStorageResidual();
            return residual;
        }
        else
            return evalLocalFluxAndSourceResidual(elemVars);
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalFluxAndSourceResidual() const
    {
        return isImplicit() ? evalLocalFluxAndSourceResidual(curElemVars_)
                            : evalLocalFluxAndSourceResidual(prevElemVars_);
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVars The element variables
     */
    ElementResidualVector evalLocalFluxAndSourceResidual(const ElementVariables& elemVars) const
    {
        return localResidual_.evalFluxAndSource(fvGeometry_.element(), fvGeometry_, elemVars, elemBcTypes_);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    ElementResidualVector evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(fvGeometry_, prevElemVars_.elemVolVars(), curElemVars_.elemVolVars());
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
        const auto& prevSol = this->assembler().prevGridVariables().dofs();
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
                prevElemVolVars.bindElement(element, fvGeometry, prevSol);
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
    const auto& curSol() const
    { return gridVars_.dofs(); }

    //! The global finite volume geometry
    const auto& gridGeometry()
    { return fvGeometry_.gridGeometry(); }

    //! The local finite volume geometry
    FVElementGeometry& fvGeometry()
    { return fvGeometry_; }

    //! The current element variables
    auto& curElemVars()
    { return curElemVars_; }

    //! The current element volume variables
    auto& curElemVolVars()
    { return curElemVars_.elemVolVars(); }

    //! The element volume variables of the provious time step
    auto& prevElemVolVars()
    { return prevElemVars_.elemVolVars(); }

    //! The element flux variables cache
    auto& elemFluxVarsCache()
    { return curElemVars_.elemFluxVarsCache(); }

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
    const auto& curElemVolVars() const
    { return curElemVars_.elemVolVars(); }

    //! The element volume variables of the provious time step
    const auto& prevElemVolVars() const
    { return prevElemVars_.elemVolVars(); }

    //! The element flux variables cache
    const auto& elemFluxVarsCache() const
    { return curElemVars_.elemFluxVarsCache(); }

    //! The element's boundary types
    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    //! The local residual for the current element
    const LocalResidual& localResidual() const
    { return localResidual_; }

protected:
    const GridVariables& gridVariables() const { return gridVars_; }
    GridVariables& gridVariables() { return gridVars_; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template<class T = TypeTag, typename std::enable_if_t<!GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    VolumeVariables& getVolVarAccess(GridVariables& gridVars, ElementVariables& elemVars, const SubControlVolume& scv)
    { return elemVars.volumeVariables(scv); }

    template<class T = TypeTag, typename std::enable_if_t<GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    VolumeVariables& getVolVarAccess(GridVariables& gridVars, ElementVariables& elemVars, const SubControlVolume& scv)
    { return gridVars.gridVolVars.volVars(scv); }

private:

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    GridVariables& gridVars_; //!< the current state of the variables

    FVElementGeometry fvGeometry_;
    ElementVariables curElemVars_;
    ElementVariables prevElemVars_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_; //!< whether the element's partitionType is ghost
};


} // end namespace Dumux

#endif
