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
 * \brief Calculates the element-wise residual of finite-volume models.
 */
#ifndef DUMUX_FV_LOCAL_RESIDUAL_HH
#define DUMUX_FV_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \brief Element-wise calculation of the residual matrix for models
 *        using a fully implicit discretization.
 *
 * \note This class defines the interface used by the assembler using
 *       static polymorphism.
 */
template<class TypeTag>
class FVLocalResidual
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
    using ElementResidualVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, NumEqVector)>;
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! the constructor for stationary problems
    FVLocalResidual() : prevSol_(nullptr) {}

    //! the constructor for instationary problems
    FVLocalResidual(std::shared_ptr<TimeLoop> timeLoop)
    : timeLoop_(timeLoop)
    , prevSol_(nullptr)
    {}

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting / postprocessing residual information.
     */
    // \{

    /*!
     * \brief Compute the storage term for the current solution.
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param element The DUNE Codim<0> entity for which the storage
     *                term ought to be calculated
     */
    ElementResidualVector evalStorage(const Problem& problem,
                                      const Element &element,
                                      const FVGridGeometry& fvGridGeometry,
                                      const GridVariables& gridVariables,
                                      const SolutionVector& sol) const
    {
        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bind(element, fvGeometry, sol);

        ElementResidualVector storage(fvGeometry.numScv());

        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localScvIdx = scv.indexInElement();
            const auto& volVars = elemVolVars[scv];
            storage[localScvIdx] = asImp().computeStorage(scv, volVars);
            storage[localScvIdx] *= scv.volume() * volVars.extrusionFactor();
        }

        return storage;
    }

    // \}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero for instationary problems.
     * \param problem The problem to solve

     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                    sub-control volumes of the element at the previous
     *                    time level
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
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        assert(prevSol_ && "no solution set for storage term evaluation");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());
        residual = 0.0;

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp().evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scv);
            asImp().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
        }

        for (auto&& scvf : scvfs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
        }

        return residual;
    }

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     * \param problem The problem to solve

     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                    sub-control volumes of the element at the previous
     *                    time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    ElementResidualVector evalStorage(const Problem& problem,
                                      const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& prevElemVolVars,
                                      const ElementVolumeVariables& curElemVolVars,
                                      const ElementBoundaryTypes &bcTypes,
                                      const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        assert(prevSol_ && "no solution set for storage term evaluation");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());
        residual = 0.0;

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp().evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scv);
        }

        return residual;
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero for stationary problem.
     * \param problem The problem to solve

     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    ElementResidualVector eval(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& curElemVolVars,
                               const ElementBoundaryTypes &bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());
        residual = 0.0;

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
        }

        for (auto&& scvf : scvfs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            asImp().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
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
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // mass balance within the element. this is the
        // \f$\frac{m}{\partial t}\f$ term if using implicit or explicit
        // euler as time discretization.
        //
        // We might need a more explicit way for
        // doing the time discretization...

        //! Compute storage with the model specific storage residual
        ResidualVector prevStorage = asImp().computeStorage(problem, scv, prevVolVars);
        ResidualVector storage = asImp().computeStorage(problem, scv, curVolVars);

        prevStorage *= prevVolVars.extrusionFactor();
        storage *= curVolVars.extrusionFactor();

        storage -= prevStorage;
        storage *= scv.volume();
        storage /= timeLoop_->timeStepSize();

        residual[scv.indexInElement()] += storage;
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
        ResidualVector source = asImp().computeSource(problem, element, fvGeometry, curElemVolVars, scv);
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

    ResidualVector evalFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        return asImp().evalFlux(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic storage derivative");
    }

    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic source derivative");
    }

    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) != DiscretizationMethods::Box, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic flux derivative for cell-centered models");
    }

    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::Box, void>
    addFluxDerivatives(JacobianMatrix& A,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic flux derivative for box models");
    }

    template<class PartialDerivativeMatrices>
    void addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& curElemVolVars,
                                     const ElementFluxVariablesCache& elemFluxVarsCache,
                                     const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic Dirichlet flux derivative");
    }

    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "analytic Robin flux derivative");
    }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { prevSol_ = &u; }

    /*!
     * \brief Return the solution that has been set as the previous one.
     */
    const SolutionVector& prevSol() const
    {
        assert(prevSol_ && "no solution set for storage term evaluation");
        return *prevSol_;
    }

    /*!
     * \brief If no solution has been set, we treat the problem as stationary.
     */
    bool isStationary() const
    { return !prevSol_; }

    // \}
protected:
    TimeLoop& timeLoop()
    { return *timeLoop_; }

    const TimeLoop& timeLoop() const
    { return *timeLoop_; }

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    std::shared_ptr<TimeLoop> timeLoop_;
    const SolutionVector* prevSol_;
};

} // end namespace Dumux

#endif
