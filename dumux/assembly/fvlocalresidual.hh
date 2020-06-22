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
 * \brief The element-wise residual for finite volume schemes
 */
#ifndef DUMUX_FV_LOCAL_RESIDUAL_HH
#define DUMUX_FV_LOCAL_RESIDUAL_HH

#include <dune/common/exceptions.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise residual for finite volume schemes
 * \note This class defines the interface used by the assembler using
 *       static polymorphism. Implementations are specialized for a certain discretization scheme
 */
template<class TypeTag>
class FVLocalResidual
{
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! the container storing all element residuals
    using ElementResidualVector = ReservedBlockVector<NumEqVector, FVElementGeometry::maxNumElementScvs>;

    //! the constructor
    FVLocalResidual(const Problem* problem,
                    const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
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
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the storage
     *                term ought to be calculated
     * \param gridGeometry The finite-volume grid geometry
     * \param gridVariables The grid variables (volume and flux variables)
     * \param sol The solution vector
     */
    ElementResidualVector evalStorage(const Problem& problem,
                                      const Element &element,
                                      const GridGeometry& gridGeometry,
                                      const GridVariables& gridVariables,
                                      const SolutionVector& sol) const
    {
        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(gridVariables.curGridVolVars());
        elemVolVars.bind(element, fvGeometry, sol);

        ElementResidualVector storage(fvGeometry.numScv());

        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localScvIdx = scv.localDofIndex();
            const auto& volVars = elemVolVars[scv];
            storage[localScvIdx] = asImp().computeStorage(problem, scv, volVars);
            storage[localScvIdx] *= Extrusion::volume(scv) * volVars.extrusionFactor();
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
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     */
    ElementResidualVector evalStorage(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& prevElemVolVars,
                                      const ElementVolumeVariables& curElemVolVars) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (auto&& scv : scvs(fvGeometry))
            asImp().evalStorage(residual, this->problem(), element, fvGeometry, prevElemVolVars, curElemVolVars, scv);

        return residual;
    }

    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache,
                                            const ElementBoundaryTypes &bcTypes) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (auto&& scv : scvs(fvGeometry))
            asImp().evalSource(residual, this->problem(), element, fvGeometry, elemVolVars, scv);

        // forward to the local residual specialized for the discretization methods
        for (auto&& scvf : scvfs(fvGeometry))
            asImp().evalFlux(residual, this->problem(), element, fvGeometry, elemVolVars, bcTypes, elemFluxVarsCache, scvf);

        return residual;
    }

    // \}


    /*!
     * \name Model specific interface
     * \note The following method are the model specific implementations of the local residual
     */
    // \{

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
     * \note has to be implemented by the model specific residual class
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a storage method!");
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scvf The sub-control volume over which we integrate the flux
     * \param elemFluxVarsCache the flux variable caches for the element's flux stencils
     *
     * \note has to be implemented by the model specific residual class
     *
     */
    NumEqVector computeFlux(const Problem& problem,
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

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
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
        // TODO: We might need a more explicit way for
        // doing the time discretization...

        //! Compute storage with the model specific storage residual
        NumEqVector prevStorage = asImp().computeStorage(problem, scv, prevVolVars);
        NumEqVector storage = asImp().computeStorage(problem, scv, curVolVars);

        prevStorage *= prevVolVars.extrusionFactor();
        storage *= curVolVars.extrusionFactor();

        storage -= prevStorage;
        storage *= Extrusion::volume(scv);
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
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the source term is integrated over
     */
    void evalSource(ElementResidualVector& residual,
                    const Problem& problem,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& curElemVolVars,
                    const SubControlVolume& scv) const
    {
        //! Compute source with the model specific storage residual
        const auto& curVolVars = curElemVolVars[scv];
        NumEqVector source = asImp().computeSource(problem, element, fvGeometry, curElemVolVars, scv);
        source *= Extrusion::volume(scv)*curVolVars.extrusionFactor();

        //! subtract source from local rate (sign convention in user interface)
        residual[scv.localDofIndex()] -= source;
    }

    /*!
     * \brief Compute the source local residual, i.e. the deviation of the
     *        source term from zero and add to the residual
     * \note This is implemented for different discretization schemes
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume averaged variables for all
     *                    sub-control volumes of the element at the current  time level
     * \param elemBcTypes the boundary types for the boundary entities of an elements
     * \param elemFluxVarsCache The flux variable caches for the element stencil
     * \param scvf The sub control volume face the flux term is integrated over
     */
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const {}

    /*!
     * \brief Compute the flux local residual, i.e. the deviation of the
     *        flux term from zero.
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param elemFluxVarsCache The flux variable caches for the element stencil
     * \param scvf The sub control volume face the flux term is integrated over
     */
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        return asImp().evalFlux(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }

    //\}

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! Compute the derivative of the storage residual
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

    //! Compute the derivative of the source residual
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

    //! Compute the derivative of the flux residual
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod != DiscretizationMethod::box, void>
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

    //! Compute the derivative of the flux residual for the box method
    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::box, void>
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

    //! Compute the derivative of the Dirichlet flux residual for cell-centered schemes
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

    //! Compute the derivative of Robin type boundary conditions ("solution dependent Neumann")
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

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_; //!< the timeloop for instationary problems
};

} // end namespace Dumux

#endif
