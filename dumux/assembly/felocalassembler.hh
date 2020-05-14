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
 * \brief An assembler for Jacobian and residual contribution
 *        per element (finite element method).
 */
#ifndef DUMUX_FE_LOCAL_ASSEMBLER_HH
#define DUMUX_FE_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>
#include <dumux/discretization/fem/elementsolution.hh>
#include "diffmethod.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A base class for all fem assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The local assembler implementation
 * \tparam useImplicitAssembly Specifies whether the time discretization is implicit or not not (i.e. explicit)
 * \todo TODO: This assumes Standard-Galerkin discretizations and non-composite function spaces.
 *             Could this be generalized to other spaces and Petrov-Galerkin?
 */
template<class Assembler>
class FELocalAssembler
{
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename GridVariables::GridGeometry;

    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridVarsLocalView = typename GridVariables::LocalView;

    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numEq = PrimaryVariables::size();

public:
    using LocalResidual = typename Assembler::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualType = typename Assembler::ResidualType;

    /*!
     * \brief Construct the local assembler from an assembler.
     */
    explicit FELocalAssembler(const Assembler& assembler)
    : assembler_(assembler)
    , feGeometry_(localView(assembler.gridGeometry()))
    , gridVarsLocalView_(localView(assembler.gridVariables()))
    {}

    //! \todo TODO Remove this function and substitute by proper time integration framework
    static constexpr bool isImplicit() { return true; }

    /*!
     * \brief Bind this local assembler to an element of the grid.
     * \param element The grid element
     * \param sol A solution vector with the primary variables at the dofs of the grid
     */
    void bind(const Element& element)
    {
        const auto& problem = assembler().problem();
        const auto& gridGeometry = problem.gridGeometry();

        feGeometry_.bind(element);
        gridVarsLocalView_.bind(element, feGeometry_);
        curElemSol_ = elementSolution(element, gridVarsLocalView().gridVariables().solutionVector(), gridGeometry);

        if (assembler().isStationaryProblem())
            localResidual_ = std::make_unique<LocalResidual>(element, feGeometry_, gridVarsLocalView_);
        else
        {
            localResidual_ = std::make_unique<LocalResidual>(element, feGeometry_, gridVarsLocalView_, &assembler().timeLoop());
            prevElemSol_ = elementSolution(element, assembler().prevSol(), gridGeometry);
        }

        elementIsGhost_ = element.partitionType() == Dune::GhostEntity;
    }

    /*!
     * \brief Evaluate the complete local residual for the current element.
     */
    ElementResidualVector evalLocalResidual() const
    {
        // error handling for incompatible combination
        if (!isImplicit() && assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // residual of ghost elements is zero
        if (elementIsGhost_)
        {
            ElementResidualVector residual(feGeometry_.feBasisLocalView().size());
            residual = 0.0;
            return residual;
        }

        return assembler().isStationaryProblem() ? localResidual_->eval(curElemSol_)
                                                 : localResidual_->eval(prevElemSol_, curElemSol_, isImplicit());
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac,
                                     ResidualType& res,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        const auto& element = localResidual_->element();
        const auto eIdxGlobal = this->assembler().gridGeometry().elementMapper().index(element);

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = evalLocalResidual();
            const auto& localView = localResidual_->feGeometry().feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else if (!elementIsGhost_)
        {
            const auto residual = assembleJacobianAndResidual_(jac, partialReassembler);
            const auto& localView = localResidual_->feGeometry().feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else
        {
            const auto& localView = localResidual_->feGeometry().feBasisLocalView();
            const auto& finiteElement = localView.tree().finiteElement();
            const auto numLocalDofs = finiteElement.localBasis().size();

            for (unsigned int i = 0; i < numLocalDofs; ++i)
            {
                const auto& localKey = finiteElement.localCoefficients().localKey(i);
                if (!isGhostEntity_(localKey.subEntity(), localKey.codim()))
                    continue;

                // set main diagonal entries for the entity
                const auto rowIdx = localView.index(i);

                // TODO: use auto and range based for loop!
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
            res[dofIdx][eqIdx] = curElemSol_[localDofIdx][pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;

            // TODO: Periodic constraints
        };

        enforceDirichletConstraints_(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);

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

        enforceDirichletConstraints_(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(ResidualType& res)
    {
        const auto residual = evalLocalResidual();

        const auto& localView = feGeometry_.feBasisLocalView();
        for (unsigned int i = 0; i < localView.size(); ++i)
            res[localView.index(i)] += residual[i];

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[dofIdx][eqIdx] = curElemSol_[localDofIdx][pvIdx] - dirichletValues[pvIdx];
        };

        enforceDirichletConstraints_(applyDirichlet);
    }

    //! The assembler
    const Assembler& assembler() const
    { return assembler_; }

    //! The finite volume geometry
    FEElementGeometry& feGeometry()
    { return feGeometry_; }

    //! The finite volume geometry
    const FEElementGeometry& feGeometry() const
    { return feGeometry_; }

    //! The local view on the grid variables
    GridVarsLocalView& gridVarsLocalView()
    { return gridVarsLocalView_; }

    //! The local view on the grid variables
    const GridVarsLocalView& gridVarsLocalView() const
    { return gridVarsLocalView_; }

    //! Return the underlying local residual
    const LocalResidual& localResidual() const
    { return *localResidual_; }

protected:
    //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints_(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        evalDirichletBoundaries_(applyDirichlet);
        // TODO: internal constraints!
        // this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries_(ApplyDirichletFunctionType applyDirichlet)
    {
        const auto& elemBcTypes = localResidual_->elemBcTypes();
        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        if (elemBcTypes.hasDirichlet())
        {
            const auto& localView = feGeometry_.feBasisLocalView();
            const auto& finiteElement = localView.tree().finiteElement();

            for (unsigned int localDofIdx = 0; localDofIdx < localView.size(); localDofIdx++)
            {
                const auto& bcTypes = elemBcTypes[localDofIdx];
                if (!bcTypes.hasDirichlet())
                    continue;

                const auto dofIdx = localView.index(localDofIdx);
                const auto& localKey = finiteElement.localCoefficients().localKey(localDofIdx);
                const auto subEntity = localKey.subEntity();
                const auto codim = localKey.codim();

                // values of dirichlet BCs
                PrimaryVariables dirichletValues = getDirichletValues_(subEntity, codim);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (bcTypes.isDirichlet(eqIdx))
                    {
                        const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                        assert(0 <= pvIdx && pvIdx < numEq);
                        applyDirichlet(dirichletValues, localDofIdx, dofIdx, eqIdx, pvIdx);
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns the Dirichlet boundary conditions for a sub entity of the element
     */
    PrimaryVariables getDirichletValues_(unsigned int subEntityIdx, unsigned int codim)
    {
        static constexpr int dim = Element::Geometry::mydimension;

        const auto& element = localResidual_->element();
        if (codim == 0)
            return getDirichletValues_(element.template subEntity<0>(subEntityIdx));
        if (codim == 1)
            return getDirichletValues_(element.template subEntity<1>(subEntityIdx));
        if constexpr (dim > 1)
            if (codim == 2)
                return getDirichletValues_(element.template subEntity<2>(subEntityIdx));
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return getDirichletValues_(element.template subEntity<3>(subEntityIdx));
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid codimension provided");
    }

    /*!
     * \brief Returns the Dirichlet boundary conditions for a sub entity of the element
     */
    template<class SubEntity>
    PrimaryVariables getDirichletValues_(const SubEntity& subEntity)
    { return assembler().problem().dirichlet(localResidual_->element(), subEntity); }

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity
     */
    bool isGhostEntity_(unsigned int subEntityIdx, unsigned int codim) const
    {
        static constexpr int dim = Element::Geometry::mydimension;

        const auto& element = localResidual_->element();
        if (codim == 0)
            return isGhostEntity_(element.template subEntity<0>(subEntityIdx));
        if (codim == 1)
            return isGhostEntity_(element.template subEntity<1>(subEntityIdx));
        if constexpr (dim > 1)
            if (codim == 2)
                return isGhostEntity_(element.template subEntity<2>(subEntityIdx));
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return isGhostEntity_(element.template subEntity<3>(subEntityIdx));
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid codimension provided");
    }

    /*!
     * \brief Returns true if an entity is ghost entity
     */
    template<class Entity>
    bool isGhostEntity_(const Entity& e) const
    { return e.partitionType() == Dune::GhostEntity; }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \return The element residual at the current solution.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidual_(JacobianMatrix& A,
                                                       const PartialReassembler* partialReassembler = nullptr)
    {
        if (diffMethod_ != DiffMethod::numeric)
            DUNE_THROW(Dune::NotImplemented, "analytic differentiation for FEM");

        using BlockType = typename JacobianMatrix::block_type;
        using PrimaryVariables = typename ElementResidualVector::value_type;
        using Scalar = typename PrimaryVariables::value_type;

        static constexpr unsigned int numPriVars = PrimaryVariables::size();
        static constexpr unsigned int numEq = BlockType::rows;
        static_assert(numPriVars == int(BlockType::cols), "row size mismatch with provided matrix block type.");

        // compute original residuals
        const auto origResiduals = evalLocalResidual();

        // create copy of element solution to undo deflections later
        const auto origElemSol = curElemSol_;

        ///////////////////////////////////////tsLocalVi///////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        const auto& localView = feGeometry().feBasisLocalView();
        ElementResidualVector partialDerivs(localView.size());
        for (unsigned int localI = 0; localI < localView.size(); ++localI)
        {
            // dof index and corresponding actual pri vars
            const auto globalI = localView.index(localI);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numPriVars; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the element solution and compute element residual
                    curElemSol_[localI][pvIdx] = priVar;
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{assembler_.problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(assembler_.problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, curElemSol_[localI][pvIdx], partialDerivs,
                                                          origResiduals, eps_(curElemSol_[localI][pvIdx], pvIdx),
                                                          numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (unsigned int localJ = 0; localJ < localView.size(); ++localJ)
                {
                    const auto globalJ = localView.index(localJ);

                    // don't add derivatives for green entities
                    if (!partialReassembler || partialReassembler->dofColor(globalJ) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of the
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof 'col'
                            A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                        }
                    }
                }

                // restore the original element solution
                curElemSol_[localI][pvIdx] = origElemSol[localI][pvIdx];

                // TODO additional dof dependencies
            }
        }

        return origResiduals;
    }

private:
    const Assembler& assembler_;          //!< reference to assembler instance
    FEElementGeometry feGeometry_;        //!< element-local view on the grid geometry
    GridVarsLocalView gridVarsLocalView_; //!< element-local view on the grid variables
    bool elementIsGhost_;          //!< whether the element's partitionType is ghost

    ElementSolution curElemSol_;   //!< element solution based on the solution of the current time step
    ElementSolution prevElemSol_;  //!< element solution based on the solution of the last time step

    std::unique_ptr<LocalResidual> localResidual_; //!< the local residual evaluating the equations per element
    DiffMethod diffMethod_{DiffMethod::numeric}; //!< the differentiation method (numeric, analytic, ...)
    // TODO: how should diff method come in? Upon construction or via setter? Or as argument to assemble() ?
};

} // end namespace Dumux

#endif
