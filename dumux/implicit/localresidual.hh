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
 * \brief Calculates the element-wise residual of fully-implicit models.
 */
#ifndef DUMUX_IMPLICIT_LOCAL_RESIDUAL_HH
#define DUMUX_IMPLICIT_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/common/capabilities.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual matrix for models
 *        using a fully implicit discretization.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class ImplicitLocalResidual
{
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    // copying the local residual class is not a good idea
    ImplicitLocalResidual(const ImplicitLocalResidual &) = delete;

    // the default constructor
    ImplicitLocalResidual() = default;

    /*!
     * \brief Initialize the local residual.
     *
     * This assumes that all objects of the simulation have been fully
     * allocated but not necessarily initialized completely.
     *
     * \param problem The representation of the physical problem to be
     *             solved.
     */
    void init(Problem &problem)
    { problemPtr_ = &problem; }


    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting residual information.
     */
    // \{

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    void eval(const Element &element)
    {
        //! Set the context for the coupling manager
        problemPtr_->couplingManager().setCouplingContext(element);

        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(problem().model().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(problem().model().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, problem().model().curSol());

        auto prevElemVolVars = localView(problem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, problem().model().prevSol());

        auto elemFluxVarsCache = localView(problem().model().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem(), element, fvGeometry);

        asImp_().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, bcTypes, elemFluxVarsCache);
    }

    /*!
     * \brief Compute the storage term for the current solution.
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param element The DUNE Codim<0> entity for which the storage
     *                term ought to be calculated
     */
    void evalStorage(const Element &element)
    {
        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(problem().model().globalFvGeometry());
        fvGeometry.bindElement(element);

        auto curElemVolVars = localView(problem().model().curGlobalVolVars());
        curElemVolVars.bindElement(element, fvGeometry, problem().model().curSol());

        auto prevElemVolVars = localView(problem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, problem().model().prevSol());

        asImp_().evalStorage_(fvGeometry, prevElemVolVars, curElemVolVars);
    }

    // !
    //  * \brief Compute the flux term for the current solution.
    //  *
    //  * \param element The DUNE Codim<0> entity for which the residual
    //  *                ought to be calculated
    //  * \param curVolVars The volume averaged variables for all
    //  *                   sub-contol volumes of the element

    // void evalFluxes(const Element &element)
    // {
    //     elemPtr_ = &element;

    //     // make sure FVElementGeometry and volume variables are bound to the element
    //     problem().model().fvGeometries_().bind(element);
    //     problem().model().curVolVars_().bind(element);
    //     problem().model().prevVolVars_().bindElement(element);
    //     problem().model().fluxVariablesCache_().bindElement(element);

    //     ElementBoundaryTypes bcTypes;
    //     bcTypes.update(problem(), element, fvGeometry_());

    //     residual_.resize(fvGeometry_().numScv);
    //     residual_ = 0;

    //     bcTypesPtr_ = &bcTypes;
    //     asImp_().evalFluxes_();
    // }

    // \}


    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the previous
     *                   time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    void eval(const Element &element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& prevElemVolVars,
              const ElementVolumeVariables& curElemVolVars,
              const ElementBoundaryTypes &bcTypes,
              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // resize the vectors for all terms
        auto numScv = fvGeometry.numScv();
        residual_.resize(numScv);
        storageTerm_.resize(numScv);

        residual_ = 0.0;
        storageTerm_ = 0.0;

        asImp_().evalVolumeTerms_(element, fvGeometry, prevElemVolVars, curElemVolVars, bcTypes);
        asImp_().evalFluxes_(element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache);
        asImp_().evalBoundary_(element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     *
     */
    PrimaryVariables computeSource(const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolume &scv)
    {
        PrimaryVariables source(0);

        // add contributions from volume flux sources
        source += this->problem().source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += this->problem().scvPointSources(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Returns the local residual for all sub-control
     *        volumes of the element.
     */
    const ElementSolutionVector &residual() const
    { return residual_; }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element.
     *
     * \param localScvIdx The local index of the sub-control volume
     */
    const PrimaryVariables &residual(const int localScvIdx) const
    { return residual_[localScvIdx]; }

    /*!
     * \brief Returns the storage term for all sub-control volumes of the
     *        element.
     */
    const ElementSolutionVector &storageTerm() const
    { return storageTerm_; }

    /*!
     * \brief Returns the storage term for a given sub-control volumes
     *        of the element.
     */
    const PrimaryVariables &storageTerm(const int localScvIdx) const
    { return storageTerm_[localScvIdx]; }

    /*!
     * \brief Return the problem we are solving. Only call this after init()!
     */
    const Problem& problem() const
    { return *problemPtr_; }

    /*!
     * \brief Return the problem we are solving. Only call this after init()!
     */
    Problem& problem()
    { return *problemPtr_; }

    // \}

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    PrimaryVariables evalFlux_(const Element &element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        return asImp_().computeFlux_(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    /*!
     * \brief Set the local residual to the storage terms of all
     *        sub-control volumes of the current element.
     */
    void evalStorage_(const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& prevElemVolVars,
                      const ElementVolumeVariables& curElemVolVars)
    {
        storageTerm_.resize(fvGeometry.numScv());
        storageTerm_ = 0;

        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localScvIdx = scv.indexInElement();
            const auto& volVars = curElemVolVars[scv];
            storageTerm_[localScvIdx] = asImp_().computeStorage(scv, volVars);
            storageTerm_[localScvIdx] *= scv.volume() * volVars.extrusionFactor();
        }
    }

    // PrimaryVariables evalSource_(const Element& element,
    //                              const FVElementGeometry& fvGeometry)
    // {
    //     PrimaryVariables source(0);
    //     const auto& fvGeometry = fvGeometry.fvElementGeometry();
    //     for (auto&& scv : scvs(fvGeometry))
    //     {
    //         source += this->problem().source(element, scv);

    //         // add contribution from possible point sources
    //         source += this->problem().scvPointSources(element, scv);
    //     }

    //     return source;
    // }

    /*!
     * \brief Add the change the source term for stationary problems
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    template<class P = Problem>
    typename std::enable_if<Dumux::Capabilities::isStationary<P>::value, void>::type
    evalVolumeTerms_(const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const ElementBoundaryTypes &bcTypes)
    {
        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localScvIdx = scv.indexInElement();
            auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();

            // subtract the source term from the local rate
            PrimaryVariables source = asImp_().computeSource(element, fvGeometry, curElemVolVars, scv);
            source *= scv.volume()*curExtrusionFactor;

            residual_[localScvIdx] -= source;
        }
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    template<class P = Problem>
    typename std::enable_if<!Dumux::Capabilities::isStationary<P>::value, void>::type
    evalVolumeTerms_(const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const ElementBoundaryTypes &bcTypes)
    {
        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localScvIdx = scv.indexInElement();

            const auto& curVolVars = curElemVolVars[scv];
            const auto& prevVolVars = prevElemVolVars[scv];

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // We might need a more explicit way for
            // doing the time discretization...
            PrimaryVariables prevStorage = asImp_().computeStorage(scv, prevVolVars);
            PrimaryVariables curStorage = asImp_().computeStorage(scv, curVolVars);

            prevStorage *= prevVolVars.extrusionFactor();
            curStorage *= curVolVars.extrusionFactor();

            storageTerm_[localScvIdx] = std::move(curStorage);
            storageTerm_[localScvIdx] -= std::move(prevStorage);
            storageTerm_[localScvIdx] *= scv.volume();
            storageTerm_[localScvIdx] /= problem().timeManager().timeStepSize();

            // add the storage term to the residual
            residual_[localScvIdx] += storageTerm_[localScvIdx];

            // subtract the source term from the local rate
            PrimaryVariables source = asImp_().computeSource(element, fvGeometry, curElemVolVars, scv);
            source *= scv.volume()*curVolVars.extrusionFactor();

            residual_[localScvIdx] -= source;
        }
    }

protected:
    ElementSolutionVector storageTerm_;
    ElementSolutionVector residual_;

private:
    Problem* problemPtr_;
};

}

#endif
