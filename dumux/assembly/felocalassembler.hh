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
 * \ingroup FEMDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (box method)
 */
#ifndef DUMUX_FEM_LOCAL_ASSEMBLER_HH
#define DUMUX_FEM_LOCAL_ASSEMBLER_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/assembly/entitycolor.hh>

#include <dumux/discretization/fem/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A base class for all fem assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The local assembler implementation
 * \tparam useImplicitAssembly Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag,
         class Assembler,
         class Implementation,
         bool useImplicitAssembly>
class FELocalAssemblerBase
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using AnsatzSpaceBasis = typename GridGeometry::FEBasis;
    using AnsatzSpaceBasisLocalView = typename AnsatzSpaceBasis::LocalView;

    using TrialSpaceBasis = typename Assembler::TrialSpaceBasis;
    using TrialSpaceBasisLocalView = typename TrialSpaceBasis::LocalView;

    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using SolutionVector = typename Assembler::ResidualType;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    /*!
     * \brief The constructor. Delegates to the general constructor.
     */
    explicit FELocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol)
    : FELocalAssemblerBase(assembler,
                           element,
                           curSol,
                           localView(assembler.gridGeometry()),
                           ElementSolution(),
                           ElementSolution(),
                           assembler.localResidual(),
                           element.partitionType() == Dune::GhostEntity)
    {}

    /*!
     * \brief The constructor. General version explicitly expecting each argument.
     */
    explicit FELocalAssemblerBase(const Assembler& assembler,
                                  const Element& element,
                                  const SolutionVector& curSol,
                                  const FEElementGeometry& feGeometry,
                                  const ElementSolution& curElemSol,
                                  const ElementSolution& prevElemSol,
                                  const LocalResidual& localResidual,
                                  const bool elementIsGhost)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
    , feGeometry_(feGeometry)
    , curElemSol_(curElemSol)
    , prevElemSol_(prevElemSol)
    , trialSpaceLocalView_(assembler.trialSpaceBasis().localView())
    , localResidual_(localResidual)
    , elementIsGhost_(elementIsGhost)
    {}

    /*!
     * \brief Returns true if the assembler considers implicit assembly.
     */
    static constexpr bool isImplicit()
    { return useImplicitAssembly; }

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element.
     */
    ElementResidualVector evalLocalResidual() const
    {
        if (!isImplicit())
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (elementIsGhost())
        {
            auto residual = ElementResidualVector(feGeometry_.feBasisLocalView().size());
            residual = 0.0;
            return residual;
        }

        // evaluate stationary local residua
        if (this->assembler().isStationaryProblem())
        {
            if (this->assembler().isStandardGalerkin())
                return localResidual_.eval(element_, feGeometry_, curElemSol_);
            else
                return localResidual_.eval(element_, feGeometry_, curElemSol_, trialSpaceLocalView_);
        }

        // instationary local residua
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Instationary FEM assembly");
        }
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
        const auto& gridGeometry = this->assembler().gridGeometry();

        // bind the local views
        feGeometry_.bind(element);
        if (!Assembler::isStandardGalerkin())
            trialSpaceLocalView_.bind(element);

        // element boundary types
        elemBcTypes_.update(this->problem(), element, feGeometry_);

        if (isImplicit())
        {
            curElemSol_ = elementSolution(element, curSol, gridGeometry);
            if (!this->assembler().isStationaryProblem())
                prevElemSol_ = elementSolution(element, prevSol, gridGeometry);
        }
        else
        {
            curElemSol_ = elementSolution(element, curSol, gridGeometry);
            prevElemSol_ = elementSolution(element, prevSol, gridGeometry);
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

    //! The finite volume geometry
    FEElementGeometry& feGeometry()
    { return feGeometry_; }

    //! The finite volume geometry
    const FEElementGeometry& feGeometry() const
    { return feGeometry_; }

    //! The current element solution
    ElementSolution& curElemSol()
    { return curElemSol_; }

    //! The current element solution
    const ElementSolution& curElemSol() const
    { return curElemSol_; }

    //! The element solution of the previous time step
    ElementSolution& prevElemSol()
    { return prevElemSol_; }

    //! The element solution of the provious time step
    const ElementSolution& prevElemSol() const
    { return prevElemSol_; }

    //! The local residual for the current element
    LocalResidual& localResidual()
    { return localResidual_; }

    //! The element's boundary types
    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    //! The element's boundary types
    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    //! The local residual for the current element
    const LocalResidual& localResidual() const
    { return localResidual_; }

    //! Return the local view on the trial space (standard galerkin case)
    template<bool sg = Assembler::isStandardGalerkin(), std::enable_if_t<sg, int> = 0>
    const AnsatzSpaceBasisLocalView& trialSpaceBasisLocalView() const
    { return feGeometry().feBasisLocalView(); }

    //! Return the local view on the trial space
    template<bool sg = Assembler::isStandardGalerkin(), std::enable_if_t<!sg, int> = 0>
    const TrialSpaceBasisLocalView& trialSpaceBasisLocalView() const
    { return trialSpaceLocalView_; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Assembler& assembler_;   //!< access pointer to assembler instance
    const Element& element_;       //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution

    FEElementGeometry feGeometry_;
    ElementSolution curElemSol_;
    ElementSolution prevElemSol_;
    ElementBoundaryTypes elemBcTypes_;
    TrialSpaceBasisLocalView trialSpaceLocalView_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_;         //!< whether the element's partitionType is ghost
};

/*!
 * \ingroup Assembly
 * \ingroup FEMDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (FEM methods)
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true>
class FELocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup FEMDiscretization
 * \brief Implementation of the FEM local assembler for implicit time discretization
 *        and numeric differentiation
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler>
class FELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>
: public FELocalAssemblerBase<TypeTag, Assembler, FELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>, true>
{
    using ThisType = FELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>;
    using ParentType = FELocalAssemblerBase<TypeTag, Assembler, ThisType, true>;

    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:

    //! Pull up the parent's constructor
    using ParentType::ParentType;

    //! Pull up parent's type definitions
    using typename ParentType::ElementResidualVector;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        this->asImp_().bindLocalViews();

        const auto eIdxGlobal = this->assembler().gridGeometry().elementMapper().index(this->element());
        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->asImp_().evalLocalResidual();

            const auto& tsLocalView = this->trialSpaceBasisLocalView();
            for (unsigned int i = 0; i < tsLocalView.size(); ++i)
                res[tsLocalView.index(i)] += residual[i];
        }
        else if (!this->elementIsGhost())
        {
            const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, partialReassembler);

            const auto& tsLocalView = this->trialSpaceBasisLocalView();
            for (unsigned int i = 0; i < tsLocalView.size(); ++i)
                res[tsLocalView.index(i)] += residual[i];
        }
        else
        {
            const auto& tsLocalView = this->trialSpaceBasisLocalView();
            const auto& tsFiniteElement = tsLocalView.tree().finiteElement();
            const auto numLocalDofs = tsFiniteElement.localBasis().size();

            for (unsigned int i = 0; i < numLocalDofs; ++i)
            {
                const auto& localKey = tsFiniteElement.localCoefficients().localKey(i);

                if (isGhostEntity_(this->element(), localKey.subEntity(), localKey.codim()))
                    continue;

                // set main diagonal entries for the entity
                const auto rowIdx = tsLocalView.index(i);

                typedef typename JacobianMatrix::block_type BlockType;
                BlockType &J = jac[rowIdx][rowIdx];
                for (int j = 0; j < BlockType::rows; ++j)
                    J[j][j] = 1.0;

                // set residual for the entity
                res[rowIdx] = 0;
            }
        }

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[dofIdx][eqIdx] = this->curElemSol()[localDofIdx][pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;

            // TODO: Periodic constraints
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            auto& row = jac[dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto residual = this->evalLocalResidual();

        const auto& tsLocalView = this->trialSpaceBasisLocalView();
        for (unsigned int i = 0; i < tsLocalView.size(); ++i)
            res[tsLocalView.index(i)] += residual[i];

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[dofIdx][eqIdx] = this->curElemSol()[localDofIdx][pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        this->asImp_().evalDirichletBoundaries(applyDirichlet);
        // TODO: internal constraints!
        // this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries(ApplyDirichletFunctionType applyDirichlet)
    {
        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        if (this->elemBcTypes().hasDirichlet())
        {
            const auto& gridGeometry = this->feGeometry().gridGeometry();
            for (const auto& is : intersections(gridGeometry.gridView(), this->element()))
            {
                const auto bcTypes = this->elemBcTypes()[is.indexInInside()];
                if (bcTypes.hasOutflow())
                    DUNE_THROW(Dune::NotImplemented, "Support for outflow BCs in FEM models");

                if (!bcTypes.hasDirichlet())
                    continue;

                // get reference elements and finite element
                const auto& tsLocalView = this->trialSpaceBasisLocalView();
                const auto& fe = tsLocalView.tree().finiteElement();

                const auto& element = this->element();
                const auto& eg = element.geometry();
                const auto refElement = ReferenceElements::general(eg.type());

                // handle Dirichlet boundaries
                for (unsigned int localDofIdx = 0; localDofIdx < tsLocalView.size(); localDofIdx++)
                {
                    const auto dofIdx = tsLocalView.index(localDofIdx);

                    const auto& localKey = fe.localCoefficients().localKey(localDofIdx);
                    auto subEntity = localKey.subEntity();
                    auto codim = localKey.codim();

                    // skip interior dofs
                    if (codim == 0)
                        continue;

                    bool found = false;
                    // try to find this local dof (on entity with known codim) on the current intersection
                    for (int j = 0; j < refElement.size(is.indexInInside(), 1, codim); j++)
                    {
                        // If j-th sub entity is the sub entity corresponding to local dof, continue and assign BC
                        if (subEntity == refElement.subEntity(is.indexInInside(), 1, j, codim))
                        {
                            // get global coordinate of this degree of freedom
                            const auto globalPos = eg.global(refElement.position(subEntity, codim));

                            // value of dirichlet BC
                            const auto dirichletValues = this->problem().dirichlet(element, is, globalPos);

                            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            {
                                if (bcTypes.isDirichlet(eqIdx))
                                {
                                    const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                                    assert(0 <= pvIdx && pvIdx < numEq);
                                    applyDirichlet(dirichletValues, localDofIdx, dofIdx, eqIdx, pvIdx);
                                }
                            }

                            // we found the dof
                            found = true; break;
                        }

                        // stop search after we found the dof
                        if (found) break;
                    }
                }
            }
        }
    }

protected:
    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        // get some aliases for convenience
        const auto& tsLocalView = this->trialSpaceBasisLocalView();
        const auto& curSol = this->curSol();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto& elemSol = this->curElemSol();

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(tsLocalView.size());

        // calculation of the derivatives
        for (unsigned int localI = 0; localI < tsLocalView.size(); ++localI)
        {
            // dof index and corresponding actual pri vars
            const auto globalI = tsLocalView.index(localI);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[localI][pvIdx] = priVar;
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals,
                                                          elemSol[localI][pvIdx],
                                                          partialDerivs,
                                                          origResiduals,
                                                          eps_(elemSol[localI][pvIdx], pvIdx),
                                                          numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (unsigned int localJ = 0; localJ < tsLocalView.size(); ++localJ)
                {
                    const auto globalJ = tsLocalView.index(localJ);

                    // don't add derivatives for green entities
                    if (!partialReassembler || partialReassembler->dofColor(globalJ) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof
                            // 'col'.
                            A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                        }
                    }
                }

                // restore the original element solution
                elemSol[localI][pvIdx] = curSol[globalI][pvIdx];

                // TODO additional dof dependencies
            }
        }

        return origResiduals;
    }

private:

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity (specialization for dim == 1)
     */
    template<class Element, std::enable_if_t<int(Element::Geometry::mydimension) == 1, int> = 0>
    bool isGhostEntity_(const Element& element, unsigned int subEntityIdx, unsigned int codim) const
    {
        if ( codim == 0 && isGhostEntity_(element.template subEntity<0>(subEntityIdx)) )
            return true;
        else if ( codim == 1 && isGhostEntity_(element.template subEntity<1>(subEntityIdx)) )
            return true;
        return false;
    }

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity (specialization for dim == 2)
     */
    template<class Element, std::enable_if_t<int(Element::Geometry::mydimension) == 2, int> = 0>
    bool isGhostEntity_(const Element& element, unsigned int subEntityIdx, unsigned int codim) const
    {
        if ( codim == 0 && isGhostEntity_(element.template subEntity<0>(subEntityIdx)) )
            return true;
        else if ( codim == 1 && isGhostEntity_(element.template subEntity<1>(subEntityIdx)) )
            return true;
        else if ( codim == 2 && isGhostEntity_(element.template subEntity<1>(subEntityIdx)) )
            return true;
        return false;
    }

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity (specialization for dim == 3)
     */
    template<class Element, std::enable_if_t<int(Element::Geometry::mydimension) == 3, int> = 0>
    bool isGhostEntity_(const Element& element, unsigned int subEntityIdx, unsigned int codim) const
    {
        if ( codim == 0 && isGhostEntity_(element.template subEntity<0>(subEntityIdx)) )
            return true;
        else if ( codim == 1 && isGhostEntity_(element.template subEntity<1>(subEntityIdx)) )
            return true;
        else if ( codim == 2 && isGhostEntity_(element.template subEntity<2>(subEntityIdx)) )
            return true;
        else if ( codim == 3 && isGhostEntity_(element.template subEntity<3>(subEntityIdx)) )
            return true;
        return false;
    }

    /*!
     * \brief Returns true if an entity is ghost entity
     */
    template<class Entity>
    bool isGhostEntity_(const Entity& e) const
    { return e.partitionType() == Dune::InteriorEntity || e.partitionType() == Dune::BorderEntity; }

}; // FELocalAssembler with numeric differentiation and implicit time discretization

} // end namespace Dumux

#endif
