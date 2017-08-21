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
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementResidualVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting / postprocessing residual information.
     */
    // \{

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    ElementResidualVector eval(const Problem& problem,
                               const Element &element,
                               const FVGridGeometry& fvGridGeometry,
                               const GridVariables& gridVariables,
                               const SolutionVector& curSol,
                               const SolutionVector& prevSol) const
    {
        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, prevSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, curElemVolVars);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem, element, fvGeometry);

        return asImp_().eval(problem, element, fvGeometry, prevElemVolVars, curElemVolVars, bcTypes, elemFluxVarsCache);
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
    // ElementResidualVector evalStorage(const Problem& problem, const Element &element)
    // {
    //     // make sure FVElementGeometry and volume variables are bound to the element
    //     auto fvGeometry = localView(problem.model().globalFvGeometry());
    //     fvGeometry.bindElement(element);

    //     auto curElemVolVars = localView(problem.model().curGlobalVolVars());
    //     curElemVolVars.bindElement(element, fvGeometry, problem.model().curSol());

    //     ElementResidualVector storage(fvGeometry.numScv());
    //     storage.resize(fvGeometry.numScv(), 0.0);

    //     // calculate the amount of conservation each quantity inside
    //     // all sub control volumes
    //     for (auto&& scv : scvs(fvGeometry))
    //     {
    //         auto localScvIdx = scv.indexInElement();
    //         const auto& volVars = elemVolVars[scv];
    //         storage[localScvIdx] = asImp_().computeStorage(scv, volVars);
    //         storage[localScvIdx] *= scv.volume() * volVars.extrusionFactor();
    //     }

    //     return storage;
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
     * \param problem The problem to solve

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
    ElementResidualVector eval(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& prevElemVolVars,
                               const ElementVolumeVariables& curElemVolVars,
                               const ElementBoundaryTypes &bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp_().evalStorage(residual, problem, element, fvGeometry, curElemVolVars, prevElemVolVars, scv);
            asImp_().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
        }

        for (auto&& scvf : scvfs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp_().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
            asImp_().evalBoundary(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
        }

        return residual;
    }

    // \}


    /*!
     * \name Model specific interface
     * \note The following method are the model specific implementations of the actual residual
     */
    // \{

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     * \note has to be implemented by the model specific residual class
     *
     */
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a storage method!");
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    ResidualVector computeSource(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolume &scv) const
    {
        ResidualVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     * \note has to be implemented by the model specific residual class
     *
     */
    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!");
    }

    // \}

    /*!
     * \name Discretization specific interface
     * \note The following method are the discretization specific wrapper methods
     */
    // \{

    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& curElemVolVars,
                     const ElementVolumeVariables& prevElemVolVars,
                     const SubControlVolume& scv) const
    {
        // const auto& curVolVars = curElemVolVars[scv];
        // const auto& prevVolVars = prevElemVolVars[scv];

        // // mass balance within the element. this is the
        // // \f$\frac{m}{\partial t}\f$ term if using implicit
        // // euler as time discretization.
        // //
        // // We might need a more explicit way for
        // // doing the time discretization...

        // //! Compute storage with the model specific storage residual
        // ResidualVector prevStorage = asImp_().computeStorage(problem, scv, prevVolVars);
        // ResidualVector storage = asImp_().computeStorage(problem, scv, curVolVars);

        // prevStorage *= prevVolVars.extrusionFactor();
        // storage *= curVolVars.extrusionFactor();

        // storage -= prevStorage;
        // storage *= scv.volume();
        // storage /= problem.timeManager().timeStepSize();

        // residual[scv.indexInElement()] += storage;
    }

    void evalSource(ElementResidualVector& residual,
                    const Problem& problem,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& curElemVolVars,
                    const SubControlVolume& scv) const
    {
        //! Compute source with the model specific storage residual
        const auto& curVolVars = curElemVolVars[scv];
        ResidualVector source = asImp_().computeSource(problem, element, fvGeometry, curElemVolVars, scv);
        source *= scv.volume()*curVolVars.extrusionFactor();

        //! subtract source from local rate (sign convention in user interface)
        residual[scv.indexInElement()] -= source;
    }

    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const {}

    void evalBoundary(ElementResidualVector& residual,
                      const Problem& problem,
                      const Element& element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const ElementBoundaryTypes& elemBcTypes,
                      const ElementFluxVariablesCache& elemFluxVarsCache,
                      const SubControlVolumeFace& scvf) const {}

    // \}

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace Dumux

#endif
