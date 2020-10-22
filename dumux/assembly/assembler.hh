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
 * \brief Assembler class for residuals and jacobian matrices
 *        for grid-based numerical schemes.
 */
#ifndef DUMUX_ASSEMBLER_HH
#define DUMUX_ASSEMBLER_HH

#include <cassert>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include "localassembler.hh"
#include "jacobianpattern.hh"

namespace Dumux {

//! Default types used for the linear system
template<class SolutionVector> struct DefaultLinearSystemTraits;

//! Default linear system types for Dune::BlockVector
template<class PrimaryVariables>
struct DefaultLinearSystemTraits<Dune::BlockVector<PrimaryVariables>>
{
private:
    static constexpr int numEq = PrimaryVariables::size();
    using Scalar = typename PrimaryVariables::value_type;
    using MatrixBlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;

public:
    using ResidualVector = Dune::BlockVector<PrimaryVariables>;
    using JacobianMatrix = Dune::BCRSMatrix<MatrixBlockType>;
};

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for grid-based numerical schemes
 * \tparam LO The local operator (evaluation of the terms of the equation)
 * \tparam diffMethod The differentiation method to compute derivatives
 * \tparam LST The linear system traits (types used for jacobians and residuals)
 */
template< class LO, DiffMethod diffMethod,
          class LST = DefaultLinearSystemTraits<typename LO::GridVariables::SolutionVector> >
class Assembler
{
    using ThisType = Assembler<LO, diffMethod, LST>;
    using GG = typename LO::GridVariables::GridGeometry;
    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the types used for the linear system
    using Scalar = typename LO::GridVariables::Scalar;
    using JacobianMatrix = typename LST::JacobianMatrix;
    using ResidualVector = typename LST::ResidualVector;
    using ResidualType = ResidualVector;

    //! export the local operator type
    using LocalOperator = LO;

    //! the local operator states the type of variables needed for evaluation
    using GridVariables = typename LO::GridVariables;

    //! export the underlying grid geometry type
    using GridGeometry = GG;

    //! export a grid-independent alias of the variables
    using Variables = GridVariables;

    //! export type used for solution vectors
    using SolutionVector = typename GridVariables::SolutionVector;

    /*!
     * \brief The Constructor from a grid geometry.
     * \param gridGeometry A grid geometry instance
     * \note This assembler class is, after construction, defined for a specific equation
     *       (given by the template argument of the LocalOperator) and a specific grid
     *       geometry - which defines the connectivity of the degrees of freedom of the
     *       underlying discretization scheme on a particular grid. The evaluation point,
     *       consisting of a particular solution/variables/parameters may vary, and therefore,
     *       an instance of the grid variables is passed to the assembly functions.
     */
    Assembler(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    /*!
     * \brief The Constructor from a grid geometry.
     * \param gridGeometry A grid geometry instance
     * \param isImplicit boolean to indicate whether an implicit or explicit pattern should be used
     * \todo TODO replace bool with time stepping scheme once available
     */
    Assembler(std::shared_ptr<const GridGeometry> gridGeometry,
              bool isImplicit)
    : gridGeometry_(gridGeometry)
    , isImplicit_(isImplicit)
    {}

    /*!
     * \brief Assembles the Jacobian matrix and the residual around the given evaluation point
     *        which is determined by the grid variables, containing all quantities required
     *        to evaluate the equations to be assembled.
     * \param gridVariables The variables corresponding to the given solution state
     * \note We assume the grid geometry on which the grid variables are defined
     *       to be the same as the one used to instantiate this class
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        resetJacobian_(partialReassembler);
        resetResidual_();

        // TODO: Remove this assert?
        assert(gridVariables.gridGeometry().numDofs() == gridGeometry().numDofs());

        using LocalAssembler = Dumux::LocalAssembler<ThisType, diffMethod>;
        assemble_([&](const Element& element)
        {
            auto fvGeometry = localView(gridGeometry());
            auto elemVars = localView(gridVariables);

            fvGeometry.bind(element);
            elemVars.bind(element, fvGeometry);

            LocalAssembler localAssembler(element, fvGeometry, elemVars);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, partialReassembler);
        });

        // TODO: Put these into discretization-specific helpers?
        enforceDirichletConstraints_(gridVariables, *jacobian_, *residual_);
        enforceInternalConstraints_(gridVariables, *jacobian_, *residual_);
        enforcePeriodicConstraints_(gridVariables, *jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the Jacobian matrix of the discrete system of equations
     *        around a given state represented by the grid variables object.
     */
    void assembleJacobian(const GridVariables& gridVariables)
    {
        resetJacobian_();

        // TODO: Remove this assert?
        assert(gridVariables.gridGeometry().numDofs() == gridGeometry().numDofs());

        using LocalAssembler = Dumux::LocalAssembler<ThisType, diffMethod>;
        assemble_([&](const Element& element)
        {
            auto fvGeometry = localView(gridGeometry());
            auto elemVars = localView(gridVariables);

            fvGeometry.bind(element);
            elemVars.bind(element, fvGeometry);

            LocalAssembler localAssembler(element, fvGeometry, elemVars);
            localAssembler.assembleJacobianAndResidual(*jacobian_);
        });

        // TODO: Put these into discretization-specific helpers?
        enforceDirichletConstraints_(gridVariables, *jacobian_, *residual_);
        enforceInternalConstraints_(gridVariables, *jacobian_, *residual_);
        enforcePeriodicConstraints_(gridVariables, *jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the internal residual vector to store the result.
     */
    void assembleResidual(const GridVariables& gridVariables)
    {
        resetResidual_();
        assembleResidual(*residual_, gridVariables);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the provided residual vector to store the result.
     */
    void assembleResidual(ResidualVector& r, const GridVariables& gridVariables) const
    {
        using LocalAssembler = Dumux::LocalAssembler<ThisType, diffMethod>;
        assemble_([&](const Element& element)
        {
            auto fvGeometry = localView(gridGeometry());
            auto elemVars = localView(gridVariables);

            fvGeometry.bind(element);
            elemVars.bind(element, fvGeometry);

            LocalAssembler localAssembler(element, fvGeometry, elemVars);
            localAssembler.assembleResidual(r);
        });
    }

    //! Will become obsolete once the new linear solvers are available
    Scalar residualNorm(const GridVariables& gridVars) const
    {
        ResidualType residual(numDofs());
        assembleResidual(residual, gridVars);

        // for box communicate the residual with the neighboring processes
        if (GridGeometry::discMethod == DiscretizationMethod::box && gridView().comm().size() > 1)
        {
            using VertexMapper = typename GridGeometry::VertexMapper;
            VectorCommDataHandleSum<VertexMapper, SolutionVector, GridGeometry::GridView::dimension>
                sumResidualHandle(gridGeometry_->vertexMapper(), residual);
            gridView().communicate(sumResidualHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }


    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (jacobian_->buildMode() == JacobianMatrix::BuildMode::unknown)
            jacobian_->setBuildMode(JacobianMatrix::random);
        else if (jacobian_->buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        // TODO: Does the bool need to be at compile time in getJacobianPattern?
        const auto occupationPattern = isImplicit_ ? getJacobianPattern<true>(gridGeometry())
                                                   : getJacobianPattern<false>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry().numDofs(); }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    ResidualVector& residual()
    { return *residual_; }

protected:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if (!residual_)
        {
            residual_ = std::make_shared<ResidualVector>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if (!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    void enforceDirichletConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

    void enforceInternalConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

    void enforcePeriodicConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

private:
    //! the grid geometry on which it is assembled
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! states if an implicit of explicit scheme is used (affects jacobian pattern)
    bool isImplicit_ = true;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualVector> residual_;
};

} // namespace Dumux

#endif
