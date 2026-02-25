// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \copydoc Dumux::LocalAssemblerBase
 */
#ifndef DUMUX_LOCAL_ASSEMBLER_BASE_HH
#define DUMUX_LOCAL_ASSEMBLER_BASE_HH

#include <dune/grid/common/gridenums.hh> // for GhostEntity

#include <dumux/common/properties.hh>

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
class LocalAssemblerBase
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = typename Assembler::SolutionVector;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using LocalResidual = std::decay_t<decltype(std::declval<Assembler>().localResidual())>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    /*!
     * \brief The constructor. Delegates to the general constructor.
     */
    explicit LocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol)
    : LocalAssemblerBase(assembler,
                         element,
                         curSol,
                         localView(assembler.gridGeometry()),
                         localView(assembler.gridVariables().curGridVars()),
                         localView(assembler.gridVariables().prevGridVars()),
                         assembler.localResidual(),
                         element.partitionType() == Dune::GhostEntity)
    {}

    /*!
     * \brief The constructor. General version explicitly expecting each argument.
     */
    explicit LocalAssemblerBase(const Assembler& assembler,
                                const Element& element,
                                const SolutionVector& curSol,
                                const FVElementGeometry& fvGeometry,
                                const ElementVariables& curElemVars,
                                const ElementVariables& prevElemVars,
                                const LocalResidual& localResidual,
                                const bool elementIsGhost)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
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
     *        element variables.
     */
    ElementResidualVector evalLocalResidual() const
    {
        if (!isImplicit())
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (elementIsGhost())
            return ElementResidualVector(0.0);

        return isImplicit() ? evalLocalResidual(curElemVars())
                            : evalLocalResidual(prevElemVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current element.
     * \param elemVars The element variables
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
     *        element variables.
     */
    ElementResidualVector evalLocalFluxAndSourceResidual() const
    {
        return isImplicit() ? evalLocalFluxAndSourceResidual(curElemVars())
                            : evalLocalFluxAndSourceResidual(prevElemVars());
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVars The element variables
     */
    ElementResidualVector evalLocalFluxAndSourceResidual(const ElementVariables& elemVars) const
    {
        return localResidual_.evalFluxAndSource(element_, fvGeometry_, elemVars, elemBcTypes_);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element variables.
     */
    ElementResidualVector evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(element_, fvGeometry_, prevElemVars_, curElemVars_);
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
        auto&& curElemVars = this->curElemVars();
        auto&& prevElemVars = this->prevElemVars();

        // bind the caches
        fvGeometry.bind(element);

        if (isImplicit())
        {
            curElemVars.bind(element, fvGeometry, curSol);
            if (!this->assembler().isStationaryProblem())
                prevElemVars.bindElement(element, fvGeometry, this->assembler().prevSol());
        }
        else
        {
            curElemVars.bindElement(element, fvGeometry, curSol);
            prevElemVars.bind(element, fvGeometry, prevSol);
        }
    }

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

    //! The current element variables
    ElementVariables& curElemVars()
    { return curElemVars_; }

    //! The element variables of the previous time step
    ElementVariables& prevElemVars()
    { return prevElemVars_; }

    //! The local residual for the current element
    LocalResidual& localResidual()
    { return localResidual_; }

    //! The element's boundary types
    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    //! The finite volume geometry
    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    //! The current element variables
    const ElementVariables& curElemVars() const
    { return curElemVars_; }

    //! The element variables of the previous time step
    const ElementVariables& prevElemVars() const
    { return prevElemVars_; }

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

    // template<class T = TypeTag, typename std::enable_if_t<!GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    // VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVars, ElementVariables& elemVars, const SubControlVolume& scv)
    // { return elemVars[scv]; }

    // template<class T = TypeTag, typename std::enable_if_t<GetPropType<T, Properties::GridVariables>::GridVolumeVariables::cachingEnabled, int> = 0>
    // VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVars, ElementVariables& elemVars, const SubControlVolume& scv)
    // { return gridVars.volVars(scv); }

private:

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution

    FVElementGeometry fvGeometry_;
    ElementVariables curElemVars_;
    ElementVariables prevElemVars_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_; //!< whether the element's partitionType is ghost
};


} // end namespace Dumux

#endif
