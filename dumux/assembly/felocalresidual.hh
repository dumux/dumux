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

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/fem/ipdata.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise residual for finite element schemes
 * \note This class defines the interface used by the assembler using
 *       static polymorphism. Implementations should inherit from this class
 *       and implement the required
 */
template<class TypeTag>
class FELocalResidual
{
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using AnsatzSpaceBasis = typename GridGeometry::FEBasis;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SecondaryVariables = typename GridVariables::SecondaryVariables;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using TimeLoop = TimeLoopBase< GetPropType<TypeTag, Properties::Scalar> >;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();

    // The flux term has size dim per equation
    using FluxTerm = Dune::FieldMatrix<Scalar, numEq, dim>;

public:
    //! the container storing the residual on all dofs of an element
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    //! the constructor
    FELocalResidual(const Problem* problem,
                    const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
    , intOrder_(getParamFromGroup<Scalar>(problem->paramGroup(), "Assembly.FEIntegrationOrder"))
    {}

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting / postprocessing residual information.
     */
    // \{

    /*!
     * \brief Compute the storage term for a given current solution.
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the storage
     *                term ought to be calculated
     * \param feGridGeometry The finite-element grid geometry
     * \param gridVariables The grid variables (secondary variables)
     * \param sol The solution vector
     */
    ElementResidualVector evalStorage(const Problem& problem,
                                      const Element &element,
                                      const GridGeometry& gridGeometry,
                                      const GridVariables& gridVariables,
                                      const SolutionVector& sol) const
    {
        DUNE_THROW(Dune::NotImplemented, "Storage term assembly");
    }

    // \}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the element local residual for instationary problems
     *        with ansatz and trial space being identical.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template<class ElementSolution>
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& prevElemSol,
                               const ElementSolution& curElemSol,
                               bool doImplicit = true) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        DUNE_THROW(Dune::NotImplemented, "Instationary FEM assembly");
    }

    /*!
     * \brief Compute the element local residual for instationary problems
     *        with a given local view on the trial space.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template< class ElementSolution, class TrialSpaceBasisLocalView >
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& prevElemSol,
                               const ElementSolution& curElemSol,
                               const TrialSpaceBasisLocalView& trialLocalView,
                               bool doImplicit = true) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        DUNE_THROW(Dune::NotImplemented, "Instationary FEM assembly");
    }

    /*!
     * \brief Compute the element local residual for stationary problems
     *        with a given local view on the trial space.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template< class ElementSolution, class TrialSpaceBasisLocalView >
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol,
                               const TrialSpaceBasisLocalView& trialLocalView) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using TrialLocalBasis = typename TrialSpaceBasisLocalView::Tree::FiniteElement::Traits::LocalBasisType;

        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;
        using IpDataTrial = FEIntegrationPointData<GlobalPosition, TrialLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();

        const auto& trialLocalBasis = trialLocalView.tree().finiteElement().localBasis();
        const auto trialNumLocalDofs = trialLocalBasis.size();

        ElementResidualVector residual(trialNumLocalDofs);

        const auto& geometry = element.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            IpDataAnsatz ipDataAnsatz(geometry, quadPoint.position(), ansatzLocalBasis);
            IpDataTrial ipDataTrial(geometry, quadPoint.position(), trialLocalBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables secVars;
            secVars.update(elemSol, problem(), element, ipDataAnsatz);

            // evaluate source term contribution
            const auto source = asImp_().computeSource(this->problem(), element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // evaluate flux term contribution
            const auto flux = asImp_().computeFlux(this->problem(), element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < trialNumLocalDofs; ++i)
                    residual[i][eqIdx] -= qWeight*secVars.extrusionFactor()
                                          *(source[eqIdx]*ipDataTrial.shapeValue(i)
                                            + flux[eqIdx]*ipDataTrial.gradN(i));
        }

        // add contribution from neumann segments
        residual += evalNeumannSegments_(element, feGeometry, elemSol, trialLocalView);

        return residual;
    }

    /*!
     * \brief Compute the element local residual for stationary problems
     *        with ansatz and trial space being identical.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template<class ElementSolution>
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();
        const auto ansatzNumLocalDofs = ansatzLocalBasis.size();

        ElementResidualVector residual(ansatzNumLocalDofs);

        const auto& geometry = element.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            IpDataAnsatz ipDataAnsatz(geometry, quadPoint.position(), ansatzLocalBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables secVars;
            secVars.update(elemSol, problem(), element, ipDataAnsatz);

            // evaluate source term contribution
            const auto source = asImp_().computeSource(this->problem(), element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // evaluate flux term contribution
            const auto flux = asImp_().computeFlux(this->problem(), element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < ansatzNumLocalDofs; ++i)
                    residual[i][eqIdx] -= qWeight*secVars.extrusionFactor()
                                          *(source[eqIdx]*ipDataAnsatz.shapeValue(i)
                                            + flux[eqIdx]*ipDataAnsatz.gradN(i));
        }

        // add contribution from neumann segments
        residual += evalNeumannSegments_(element, feGeometry, elemSol);

        return residual;
    }

    // \}

    /*!
     * \name Model specific interface
     * \note The following method are the model specific implementations of the local residual
     */
    // \{

    //! TODO: Storage interface

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param feGeometry The finite-element geometry
     * \param elemSol The element solution vector
     * \param ipData The trial and ansatz space shape function values/gradients
     *               evaluated at the integration point
     * \param secVars The secondary variables evaluated at the integration point
     */
     template<class IpData, class ElementSolution>
     NumEqVector computeSource(const Problem& problem,
                               const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol,
                               const IpData& ipData,
                               const SecondaryVariables& secVars) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, feGeometry, elemSol, ipData, secVars);

        // TODO: add contribution from possible point sources

        return source;
    }

    /*!
     * \brief Calculate the flux term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param feGeometry The finite-element geometry
     * \param elemSol The element solution vector
     * \param ipData The trial and ansatz space shape function values/gradients
     *               evaluated at the integration point
     * \param secVars The secondary variables evaluated at the integration point
     */
    template<class IpData, class ElementSolution>
    FluxTerm computeFlux(const Problem& problem,
                         const Element& element,
                         const FEElementGeometry& feGeometry,
                         const ElementSolution& elemSol,
                         const IpData& ipData,
                         const SecondaryVariables& secVars) const
    {
        DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!");
    }

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces

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

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:

    /*!
     * \brief Compute the contributions of Neumann boundary conditions on the
     *        local residual with a given local view on the trial space.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template< class ElementSolution, class TrialSpaceBasisLocalView >
    ElementResidualVector evalNeumannSegments_(const Element& element,
                                               const FEElementGeometry& feGeometry,
                                               const ElementSolution& elemSol,
                                               const TrialSpaceBasisLocalView& trialLocalView) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using TrialLocalBasis = typename TrialSpaceBasisLocalView::Tree::FiniteElement::Traits::LocalBasisType;

        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;
        using IpDataTrial = FEIntegrationPointData<GlobalPosition, TrialLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();

        const auto& trialLocalBasis = trialLocalView.tree().finiteElement().localBasis();
        const auto trialNumLocalDofs = trialLocalBasis.size();

        ElementResidualVector residual(trialNumLocalDofs);

        // integrate boundary contribution
        const auto geometry = element.geometry();
        for (const auto& is : intersections(feGeometry.gridGeometry().gridView(), element))
        {
            // only handle faces on the boundary
            if (!is.boundary())
                continue;

            // only treat faces with neumann boundary conditions
            auto bcTypes = problem().boundaryTypes(element, is);
            if (!bcTypes.hasNeumann())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            auto insideGeom = is.geometryInInside();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(insideGeom.type(), intOrder_-1);

            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of inside element
                auto local = insideGeom.global(quadPoint.position());

                // evaluate basis functions of all all element vertices for quadrature point
                IpDataAnsatz ipDataAnsatz(geometry, local, ansatzLocalBasis);
                IpDataTrial ipDataTrial(geometry, local, trialLocalBasis);

                // evaluate secondary variables
                SecondaryVariables secVars;
                secVars.update(elemSol, this->problem(), element, ipDataAnsatz);

                // evaluate neumann boundary condition
                auto neumannFlux = problem().neumann(element, is, elemSol, ipDataAnsatz, secVars);

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight();
                qWeight *= is.geometry().integrationElement(quadPoint.position());
                qWeight *= secVars.extrusionFactor();

                // add entries to residual vector
                for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    for (unsigned int i = 0; i < trialNumLocalDofs; ++i)
                        if (bcTypes.isNeumann(eqIdx))
                            residual[i][eqIdx] += ipDataTrial.shapeValue(i)*qWeight*neumannFlux[eqIdx];
            }
        }

        return residual;
    }

    /*!
     * \brief Compute the contributions of Neumann boundary conditions on the
     *        local residual with trial and ansatz space being identical.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template< class ElementSolution >
    ElementResidualVector evalNeumannSegments_(const Element& element,
                                               const FEElementGeometry& feGeometry,
                                               const ElementSolution& elemSol) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();
        const auto ansatzNumLocalDofs = ansatzLocalBasis.size();

        ElementResidualVector residual(ansatzNumLocalDofs);

        // integrate boundary contribution
        const auto geometry = element.geometry();
        for (const auto& is : intersections(feGeometry.gridGeometry().gridView(), element))
        {
            // only handle faces on the boundary
            if (!is.boundary())
                continue;

            // only treat faces with neumann boundary conditions
            auto bcTypes = problem().boundaryTypes(element, is);
            if (!bcTypes.hasNeumann())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            auto insideGeom = is.geometryInInside();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(insideGeom.type(), intOrder_-1);

            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of inside element
                auto local = insideGeom.global(quadPoint.position());

                // evaluate basis functions of all all element vertices for quadrature point
                IpDataAnsatz ipDataAnsatz(geometry, local, ansatzLocalBasis);

                // evaluate secondary variables
                SecondaryVariables secVars;
                secVars.update(elemSol, this->problem(), element, ipDataAnsatz);

                // evaluate neumann boundary condition
                auto neumannFlux = problem().neumann(element, is, elemSol, ipDataAnsatz, secVars);

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight();
                qWeight *= is.geometry().integrationElement(quadPoint.position());
                qWeight *= secVars.extrusionFactor();

                // add entries to residual vector
                for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    for (unsigned int i = 0; i < ansatzNumLocalDofs; ++i)
                        if (bcTypes.isNeumann(eqIdx))
                            residual[i][eqIdx] += ipDataAnsatz.shapeValue(i)*qWeight*neumannFlux[eqIdx];
            }
        }

        return residual;
    }

    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_; //!< the timeloop for instationary problems
    int intOrder_; //!< The order used for the quadrature rule to evaluate the residual
};

} // end namespace Dumux

#endif
