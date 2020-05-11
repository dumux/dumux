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
 * \brief The element-wise residual for finite element schemes
 */
#ifndef DUMUX_FE_LOCAL_RESIDUAL_HH
#define DUMUX_FE_LOCAL_RESIDUAL_HH

#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/discretization/fem/ipdata.hh>
#include <dumux/discretization/fem/elementboundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise residual for finite element schemes
 * \note This class defines the interface used by the assembler using
 *       static polymorphism. Implementations should inherit from this class
 *       and implement the required functions.
 * \todo TODO: This currently assumes Standard Galerkin type fem schemes as
 *             well as non-composite function space bases.
 */
template<class GridVariables, class Impl>
class FELocalResidual
{
    using Problem = typename GridVariables::Problem;

    using GridVarsLocalView = typename GridVariables::LocalView;
    using IpVariables = typename GridVariables::IntegrationPointVariables;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename PrimaryVariables::value_type;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using TimeLoop = TimeLoopBase<Scalar>;
    using NumEqVector = typename ProblemTraits<Problem>::NumEqVector;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
    using ElemBoundaryTypes = FEElementBoundaryTypes<BoundaryTypes>;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();

public:
    //! the container storing the residual on all dofs of an element
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    // The flux term has size dim per equation
    using FluxTerm = Dune::FieldMatrix<Scalar, numEq, dim>;

    //! the constructor for instationary setups
    //! \note The grid geometry/grid variables local views are expected to be bound to the given element
    FELocalResidual(const Element& element,
                    const FEElementGeometry& feGeometry,
                    const GridVarsLocalView& gridVarsLocalView,
                    const TimeLoop* timeLoop)
    : element_(&element)
    , feGeometry_(&feGeometry)
    , gridVarsLocalView_(&gridVarsLocalView)
    , timeLoop_(timeLoop)
    {
        elemBcTypes_.update(problem_(), element, feGeometry);
        setDefaultIntegrationOrders_();
    }

    //! the constructor for stationary setups
    //! \note The grid geometry/grid variables local views are expected to be bound to the given element
    FELocalResidual(const Element& element,
                    const FEElementGeometry& feGeometry,
                    const GridVarsLocalView& gridVarsLocalView)
    : FELocalResidual(element, feGeometry, gridVarsLocalView, nullptr)
    {}

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting / postprocessing residual information.
     */
    // \{

    /*!
     * \brief Compute the storage term for a given current solution.
     * \todo TODO This interface is defined in FVLocalResidual,
     *            does this make sense to have here?
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param element The grid element
     * \param gridGeometry The grid geometry
     * \param gridVariables The grid variables (primary/secondary variables)
     * \param sol The solution vector
     */
    template<class SolutionVector>
    ElementResidualVector evalStorage(const Element& element,
                                      const GridGeometry& gridGeometry,
                                      const GridVariables& gridVariables,
                                      const SolutionVector& sol) const
    { DUNE_THROW(Dune::NotImplemented, "Storage term assemblyfor FEM"); }

    // \}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the element local residual for instationary problems.
     * \todo TODO Doc me.
     * \todo TODO Better way to distinguish implicit/explicit
     *            That depends on the time scheme developments on master
     */
    template<class ElementSolution>
    ElementResidualVector eval(const ElementSolution& prevElemSol,
                               const ElementSolution& curElemSol,
                               bool doImplicit = true) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");

        const auto& basisLocalView = feGeometry_->feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector residual(numLocalDofs);

        // choose elemSol to evaluate sources/fluxes depending on implicit/explicit
        const auto& elemSol = doImplicit ? curElemSol : prevElemSol;

        const auto& geometry = element().geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            FEIntegrationPointData ipData(geometry, quadPoint.position(), localBasis);

            // calculate variables for the previous and the current solution at the ip
            IpVariables curIpVars, prevIpVars;
            curIpVars.update(curElemSol, problem_(), element(), ipData);
            prevIpVars.update(prevElemSol, problem_(), element(), ipData);

            // evaluate storage term
            auto storage = asImp_().computeStorage(curElemSol, ipData, curIpVars);
            storage -= asImp_().computeStorage(prevElemSol, ipData, prevIpVars);
            storage /= timeLoop_->timeStepSize();

            // choose secondary variables depending on implicit/explicit
            const auto& ipVars = doImplicit ? curIpVars : prevIpVars;

            // evaluate source term contribution
            const auto source = asImp_().computeSource(elemSol, ipData, ipVars);

            // evaluate flux term contribution
            const auto flux = asImp_().computeFlux(elemSol, ipData, ipVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                    residual[i][eqIdx] -= qWeight*ipVars.extrusionFactor()
                                          *( (storage[eqIdx]+source[eqIdx])*ipData.shapeValue(i)
                                             + flux[eqIdx]*ipData.gradN(i) );
        }

        // add contribution from neumann segments
        residual += evalNeumannSegments(elemSol);

        return residual;
    }

    /*!
     * \brief Compute the element local residual for stationary problems
     * \todo elemSol The element solution vector
     */
    template<class ElementSolution>
    ElementResidualVector eval(const ElementSolution& elemSol) const
    {
        const auto& basisLocalView = feGeometry_->feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector residual(numLocalDofs);

        const auto& geometry = element().geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            FEIntegrationPointData ipData(geometry, quadPoint.position(), localBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            IpVariables ipVars;
            ipVars.update(elemSol, problem_(), element(), ipData);

            // evaluate source term contribution
            const auto source = asImp_().computeSource(elemSol, ipData, ipVars);

            // evaluate flux term contribution
            const auto flux = asImp_().computeFlux(elemSol, ipData, ipVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                    residual[i][eqIdx] -= qWeight*ipVars.extrusionFactor()
                                          *(source[eqIdx]*ipData.shapeValue(i)
                                            + flux[eqIdx]*ipData.gradN(i));
        }

        // add contribution from neumann segments
        residual += evalNeumannSegments(elemSol);

        return residual;
    }

    /*!
     * \brief Compute the contributions of Neumann boundary conditions on the local residual.
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    ElementResidualVector evalNeumannSegments(const ElementSolution& elemSol) const
    {
        const auto& basisLocalView = feGeometry_->feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector residual(numLocalDofs);
        residual = 0.0;

        // skip the rest if element is not connected to Neumann boundaries
        if (!elemBcTypes_.hasNeumann())
            return residual;

        // integrate Neumann boundary contribution
        const auto geometry = element().geometry();
        for (const auto& is : intersections(feGeometry_->gridGeometry().gridView(), element()))
        {
            // only handle faces on the boundary
            if (!is.boundary())
                continue;

            // TODO: If none of the dofs living on this intersection has neumann BC defined, skip rest

            // select quadrature rule for intersection faces (dim-1)
            auto insideGeom = is.geometryInInside();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(insideGeom.type(), intOrderBoundary_);
            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of inside element
                auto local = insideGeom.global(quadPoint.position());

                // evaluate basis functions of all element vertices for quadrature point
                FEIntegrationPointData ipData(geometry, local, localBasis);

                // evaluate primary/secondary variables at integration point
                IpVariables ipVars;
                ipVars.update(elemSol, problem_(), element(), ipData);

                // evaluate neumann boundary condition
                const auto neumannFlux = problem_().neumann(element(), is, elemSol, ipData, ipVars);

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight();
                qWeight *= is.geometry().integrationElement(quadPoint.position());
                qWeight *= ipVars.extrusionFactor();

                // add entries to residual vector
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                    for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        if (elemBcTypes_[i].isNeumann(eqIdx))
                            residual[i][eqIdx] += ipData.shapeValue(i)*qWeight*neumannFlux[eqIdx];
            }
        }

        return residual;
    }

    // \}

    /*!
     * \name Model specific interface
     * \note The following method are the model specific implementations of the local residual
     */
    // \{

    /*!
     * \brief Calculate the storage term of the equation
     * \param elemSol The element solution vector
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
     template<class ElementSolution, class IpData>
     NumEqVector computeStorage(const ElementSolution& elemSol,
                                const IpData& ipData,
                                const IpVariables& ipVars) const
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a storage method!"); }

    /*!
     * \brief Calculate the source term of the equation
     * \param elemSol The element solution vector
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
     template<class ElementSolution, class IpData>
     NumEqVector computeSource(const ElementSolution& elemSol,
                               const IpData& ipData,
                               const IpVariables& ipVars) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem_().source(element(), feGeometry(), elemSol, ipData, ipVars);

        // TODO: add contribution from possible point sources

        return source;
    }

    /*!
     * \brief Calculate the flux term of the equation
     * \param elemSol The element solution vector
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
    template<class ElementSolution, class IpData>
    FluxTerm computeFlux(const ElementSolution& elemSol,
                         const IpData& ipData,
                         const IpVariables& ipVars) const
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!"); }

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces

    //\}

    /*!
     * \name Return functions for underlying data structures
     */
    // \{

    /*!
     * \brief Return a reference to the underlying grid element
     */
    const Element& element() const
    { return *element_; }

    /*!
     * \brief Return the local view on the grid geometry
     */
    const FEElementGeometry& feGeometry() const
    { return *feGeometry_; }

    /*!
     * \brief Return reference to the grid variables
     */
    const GridVarsLocalView& gridVariablesLocalView() const
    { return *gridVarsLocalView_; }

    /*!
     * \brief Return reference to the element boundary types
     */
    const ElemBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    /*!
     * \brief The timeloop for instationary problems
     * \note Calling this for stationary residuals leads to undefined behaviour
     */
    const TimeLoop& timeLoop() const
    {
        assert(!isStationary() && "No time loop set!");
        return *timeLoop_;
    }

    /*!
     * \brief Returns true if the residual is stationary
     */
    bool isStationary() const
    { return !timeLoop_; }

    /*!
     * \brief Set the integration order to be used for volume integrals.
     * \param intOrder The integration order
     */
    void setIntegrationOrder(unsigned int intOrder)
    {
        intOrder_ = intOrder;
    }

    /*!
     * \brief Set the integration order to be used for boundary integrals.
     * \param intOrderBoundary The boundary integration order
     */
    void setBoundaryIntegrationOrder(unsigned int intOrderBoundary)
    {
        intOrderBoundary_ = intOrderBoundary;
    }

    // \}

protected:
    //! return reference to implementation
    Impl& asImp_()
    { return *static_cast<Impl*>(this); }

    //! return const reference to implementation
    const Impl& asImp_() const
    { return *static_cast<const Impl*>(this); }

    //! return reference to the underlying problem
    const Problem& problem_() const
    { return gridVariablesLocalView().gridVariables().problem(); }

private:
    //! obtains the integration orders from the input file or default values
    void setDefaultIntegrationOrders_()
    {
        static const auto basisOrder = feGeometry_->feBasisLocalView().tree().finiteElement().localBasis().order();
        static const auto intOrder = getParamFromGroup<unsigned int>(problem_().paramGroup(), "Assembly.FEIntegrationOrder", basisOrder+1);
        static const auto bOrder = getParamFromGroup<unsigned int>(problem_().paramGroup(), "Assembly.FEBoundaryIntegrationOrder", intOrder);

        intOrder_ = intOrder;
        intOrderBoundary_ = bOrder;
    }

    const Element* element_;                     //!< pointer to the element for which the residual is computed
    const FEElementGeometry* feGeometry_;        //!< the local view on the finite element grid geometry
    const GridVarsLocalView* gridVarsLocalView_; //!< the local view on the grid variables
    const TimeLoop* timeLoop_;                   //!< the time loop for instationary problems
    ElemBoundaryTypes elemBcTypes_;              //!< the boundary types defined for all dofs of the element

    unsigned int intOrder_;               //!< Integration order used for volume integrals
    unsigned int intOrderBoundary_;       //!< Integration order used for integration of boundary conditions
};

} // end namespace Dumux

#endif
