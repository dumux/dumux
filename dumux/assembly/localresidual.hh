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
#include <dumux/common/pointsources.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/concepts.hh>

namespace Dumux::Experimental {

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
    using ElementDiscretization = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Extrusion = Extrusion_t<GridGeometry>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! the container storing all element residuals
    using ElementResidualVector = ReservedBlockVector<NumEqVector, Dumux::Detail::LocalDofs::maxNumLocalDofs<ElementDiscretization>()>;

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
     * \param elemDisc The element discretization
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current  time level
     */
    ElementResidualVector evalStorage(const Element& element,
                                      const ElementDiscretization& elemDisc,
                                      const ElementVariables& prevElemVars,
                                      const ElementVariables& curElemVars) const
    {
        assert(!this->isStationary() && "no time loop set for storage term evaluation");

        // initialize the residual vector for all local dofs in this element
        ElementResidualVector residual(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        if constexpr (Concepts::FVElementDiscretization<ElementDiscretization>)
        {
            for (const auto& scv : scvs(elemDisc))
                this->asImp().evalStorage(residual, this->problem(), element, elemDisc, prevElemVars, curElemVars, scv);
        }

        // allow for additional contributions (e.g. hybrid CVFE / FE schemes)
        this->asImp().addToElementStorageResidual(residual, this->problem(), element, elemDisc, prevElemVars, curElemVars);

        return residual;
    }

    /*!
     * \brief Compute the flux and source
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element at the current time level
     */
    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const ElementDiscretization& elemDisc,
                                            const ElementVariables& elemVars) const
    {
        // initialize the residual vector for all local dofs in this element
        ElementResidualVector residual(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));

        if constexpr (Concepts::FVElementDiscretization<ElementDiscretization>)
        {
            // evaluate the volume terms (storage + source terms)
            // forward to the local residual specialized for the discretization methods
            for (const auto& scv : scvs(elemDisc))
                this->asImp().evalSource(residual, this->problem(), element, elemDisc, elemVars, scv);

            // forward to the local residual specialized for the discretization methods
            // TODO: Provide scvfs range only over interior scvfs
            for (auto&& scvf : scvfs(elemDisc))
                if(!scvf.boundary())
                    this->asImp().evalFlux(residual, this->problem(), element, elemDisc, elemVars, scvf);
        }

        // allow for additional contributions (e.g. hybrid CVFE / FE schemes)
        this->asImp().addToElementFluxAndSourceResidual(residual, this->problem(), element, elemDisc, elemVars);

        // add boundary flux contributions
        this->asImp().addBoundaryFluxIntegral(residual, this->problem(), elemDisc, elemVars);

        return residual;
    }

    //! add additional storage contributions (e.g. hybrid CVFE or FE schemes)
    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const ElementDiscretization& elemDisc,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {}

    //! add additional flux and source contributions (e.g. hybrid CVFE or FE schemes)
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& curElemVars) const
    {}

    //! add boundary flux contributions
    void addBoundaryFluxIntegral(ElementResidualVector& residual,
                                 const Problem& problem,
                                 const ElementDiscretization& elemDisc,
                                 const ElementVariables& elemVars) const
    {
        if(!elemDisc.hasBoundaryFaces())
            return;

        for (const auto& boundaryFace : boundaryFaces(elemDisc))
        {
            const auto& bcTypes = problem.boundaryTypes(elemDisc, boundaryFace);
            if(!bcTypes.hasFluxBoundary())
                continue;

            problem.addFEBoundaryFluxIntegral(residual, elemDisc, elemVars, boundaryFace, bcTypes);

            for(const auto& scvf : scvfs(elemDisc, boundaryFace))
                problem.addFVBoundaryFluxIntegral(residual, elemDisc, elemVars, scvf, bcTypes);
        }
    }

    // \}

    /*!
     * \name Model specific interface
     * \note The following methods are the model specific implementations of the local residual
     */
    // \{

    /*!
     * \brief Calculate the source term integral of the equation
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \param isPreviousTimeLevel If set to true, the storage integral term is evaluated on the previous time level.
     * \note has to be implemented by the model specific residual class
     *
     */
    template<class SubControlVolume>
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a storageIntegral method!");
    }

    /*!
     * \brief Calculate the source term integral of the equation
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    template<class SubControlVolume>
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
            source += qpData.weight() * problem.source(elemDisc, elemVars, qpData.ipData());

        source *= elemVars[scv].extrusionFactor();

        // add contribution from possible point sources
        const auto& pointSources = problem.pointSources();
        if (!pointSources.empty())
            for (const auto& context : pointSources.contexts(elemDisc, scv))
            {
                auto psValues = pointSources.eval(elemDisc, elemVars, context);
                source += psValues;
            }

        return source;
    }

    /*!
     * \brief Calculate the flux integral of the equation
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables associated with the element stencil
     * \param scvf The sub-control volume over which we integrate the flux
     *
     * \note has to be implemented by the model specific residual class
     *
     */
    template<class SubControlVolumeFace>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a fluxIntegral method!");
    }

    // \}

    /*!
     * \name Discretization specific interface
     * \note The following methods are the discretization specific wrapper methods
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
     * \param elemDisc The element discretization
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    template<class SubControlVolume>
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const ElementDiscretization& elemDisc,
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
        NumEqVector prevStorage = this->asImp().storageIntegral(elemDisc, prevElemVars, scv, /*previous time level?*/true);
        NumEqVector storage = this->asImp().storageIntegral(elemDisc, curElemVars, scv, /*previous time level?*/false);

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
     * \param elemDisc The element discretization
     * \param curElemVars The variables for all local dofs of the element at the current time level
     * \param scv The sub control volume the source term is integrated over
     */
    template<class SubControlVolume>
    void evalSource(ElementResidualVector& residual,
                    const Problem& problem,
                    const Element& element,
                    const ElementDiscretization& elemDisc,
                    const ElementVariables& curElemVars,
                    const SubControlVolume& scv) const
    {
        //! Compute source with the model specific residual
        NumEqVector source = this->asImp().sourceIntegral(elemDisc, curElemVars, scv);
        //! subtract source from local rate (sign convention in user interface)
        residual[scv.localDofIndex()] -= source;
    }

    /*!
     * \brief Compute the fluxes of the local residual
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element at the current time level
     * \param scvf The sub control volume face the flux term is integrated over
     */
    template<class SubControlVolumeFace>
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const ElementDiscretization& elemDisc,
                  const ElementVariables& elemVars,
                  const SubControlVolumeFace& scvf) const {}

    /*!
     * \brief Compute the fluxes of the local residual
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element at the current time level
     * \param scvf The sub control volume face the flux term is integrated over
     */
    template<class SubControlVolumeFace>
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const ElementDiscretization& elemDisc,
                         const ElementVariables& elemVars,
                         const SubControlVolumeFace& scvf) const
    {
        return asImp().evalFlux(problem, element, elemDisc, elemVars, scvf);
    }

    //\}

    /*!
     * \name Interfaces accessed by local residual implementations
     */
    // \{

    //! the problem
    const Problem& problem() const
    { return *problem_; }

    //! the time loop for instationary problems
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
    const TimeLoop* timeLoop_; //!< the time loop for instationary problems
};

} // end namespace Dumux::Experimental

#endif
