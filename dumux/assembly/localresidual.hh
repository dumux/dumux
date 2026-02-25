// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief The element-wise residual for grid-based discretization schemes
 */
#ifndef DUMUX_BASE_LOCAL_RESIDUAL_HH
#define DUMUX_BASE_LOCAL_RESIDUAL_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise residual for grid-based discretization schemes
 * \note This class defines the interface used by the assembler using
 *       static polymorphism. Implementations are specialized for a certain discretization scheme
 */
template<class TypeTag>
class LocalResidual
{
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Extrusion = Extrusion_t<GridGeometry>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! the container storing all element residuals
    using ElementResidualVector = ReservedBlockVector<NumEqVector, Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>()>;

    //! the constructor
    LocalResidual(const Problem* problem,
                  const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
    {}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current  time level
     */
    ElementResidualVector evalStorage(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVariables& prevElemVars,
                                      const ElementVariables& curElemVars) const
    {
        assert(!this->isStationary() && "no time loop set for storage term evaluation");

        // initialize the residual vector for all local dofs in this element
        ElementResidualVector residual(Detail::LocalDofs::numLocalDofs(fvGeometry));

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(fvGeometry))
            this->asImp().evalStorage(residual, this->problem(), element, fvGeometry, prevElemVars, curElemVars, scv);

        // allow for additional contributions (e.g. hybrid CVFE / FE schemes)
        this->asImp().addToElementStorageResidual(residual, this->problem(), element, fvGeometry, prevElemVars, curElemVars);

        return residual;
    }

    /*!
     * \brief Compute the flux and source
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element at the current time level
     * \param elemFluxVarsCache The element flux variables cache
     * \param bcTypes The element boundary types
     */
    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVariables& elemVars,
                                            const ElementBoundaryTypes& bcTypes) const
    {
        // initialize the residual vector for all local dofs in this element
        ElementResidualVector residual(Detail::LocalDofs::numLocalDofs(fvGeometry));

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(fvGeometry))
            this->asImp().evalSource(residual, this->problem(), element, fvGeometry, elemVars, scv);

        // forward to the local residual specialized for the discretization methods
        for (auto&& scvf : scvfs(fvGeometry))
            this->asImp().evalFlux(residual, this->problem(), element, fvGeometry, elemVars, bcTypes, scvf);

        // allow for additional contributions (e.g. hybrid CVFE / FE schemes)
        this->asImp().addToElementFluxAndSourceResidual(residual, this->problem(), element, fvGeometry, elemVars, bcTypes);

        return residual;
    }

    //! add additional storage contributions (e.g. hybrid CVFE or FE schemes)
    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {}

    //! add additional flux and source contributions (e.g. hybrid CVFE or FE schemes)
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& curElemVars,
                                           const ElementBoundaryTypes& bcTypes) const
    {}

    // \}


    /*!
     * \name Model specific interface
     * \note The following method are the model specific implementations of the local residual
     */
    // \{

    /*!
     * \brief Calculate the source term integral of the equation
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note has to be implemented by the model specific residual class
     *
     */
    template<class SubControlVolume>
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a storageIntegral method!");
    }

    /*!
     * \brief Calculate the source term integral of the equation
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    template<class SubControlVolume>
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
            source += qpData.weight() * problem.source(fvGeometry, elemVars, qpData.ipData());

        // ToDo: point source data with ipData
        // add contribution from possible point sources
        if (!problem.pointSourceMap().empty())
            source += Extrusion::volume(fvGeometry, scv) * problem.scvPointSources(fvGeometry.element(), fvGeometry, elemVars, scv);

        source *= elemVars[scv].extrusionFactor();

        return source;
    }

    /*!
     * \brief Calculate the flux integral of the equation
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables associated with the element stencil
     * \param scvf The sub-control volume over which we integrate the flux
     *
     * \note has to be implemented by the model specific residual class
     *
     */
    template<class SubControlVolumeFace>
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a fluxIntegral method!");
    }

    // \}

    /*!
     * \name Discretization specific interface
     * \note The following method are the discretization specific wrapper methods
     */
    // \{

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    template<class SubControlVolume>
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVariables& prevElemVars,
                     const ElementVariables& curElemVars,
                     const SubControlVolume& scv) const
    {
        // mass balance within the element. this is the
        // \f$\frac{\partial S}{\partial t}\f$ term if using implicit or explicit
        // euler as time discretization.
        //
        // TODO: We might need a more explicit way for
        // doing the time discretization...

        //! Compute storage with the model specific storage residual
        NumEqVector prevStorage = this->asImp().storageIntegral(fvGeometry, prevElemVars, scv, /*previous time level?*/true);
        NumEqVector storage = this->asImp().storageIntegral(fvGeometry, curElemVars, scv, /*previous time level?*/false);

        storage -= prevStorage;
        storage /= timeLoop_->timeStepSize();

        residual[scv.localDofIndex()] += storage;
    }

    /*!
     * \brief Compute the source local residual, i.e. the deviation of the
     *        source term from zero.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param curElemVars The variables for all local dofs of the element at the current time level
     * \param scv The sub control volume the source term is integrated over
     */
    template<class SubControlVolume>
    void evalSource(ElementResidualVector& residual,
                    const Problem& problem,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVariables& curElemVars,
                    const SubControlVolume& scv) const
    {
        //! Compute source with the model specific residual
        NumEqVector source = this->asImp().sourceIntegral(fvGeometry, curElemVars, scv);
        //! subtract source from local rate (sign convention in user interface)
        residual[scv.localDofIndex()] -= source;
    }

    /*!
     * \brief Compute the fluxes of the local residual
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element at the current  time level
     * \param elemBcTypes the boundary types for the boundary entities of an elements
     * \param scvf The sub control volume face the flux term is integrated over
     */
    template<class SubControlVolumeFace>
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVariables& elemVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const SubControlVolumeFace& scvf) const {}

    /*!
     * \brief Compute the fluxes of the local residual
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element at the current  time level
     * \param scvf The sub control volume face the flux term is integrated over
     */
    template<class SubControlVolumeFace>
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVariables& elemVars,
                         const SubControlVolumeFace& scvf) const
    {
        return asImp().evalFlux(problem, element, fvGeometry, elemVars, scvf);
    }

    //\}

    /*!
     * \name Interfaces accessed by local residual implementations
     */
    // \{

    //! the problem
    const Problem& problem() const
    { return *problem_; }

    //! the timeloop for instationary problems
    //! calling this for stationary leads to undefined behaviour
    const TimeLoop& timeLoop() const
    { return *timeLoop_; }

    //! returns true if the residual is stationary
    bool isStationary() const
    { return !timeLoop_; }

    // \}

protected:
    Implementation& asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_; //!< the timeloop for instationary problems
};

} // end namespace Dumux

#endif
