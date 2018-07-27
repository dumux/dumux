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
 * \brief Base class for fully-implicit finite-element models
 */
#ifndef DUMUX_IMPLICIT_FEM_MODEL_HH
#define DUMUX_IMPLICIT_FEM_MODEL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/bvector.hh>

#include "properties.hh"
#include "localresidual.hh"

#include <dumux/parallel/vertexhandles.hh>


namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class FemImplicitModel
{
    // The local jacobian needs to be able to modify objects during the assembly
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using Implementation = typename GET_PROP_TYPE(TypeTag, Model);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FEBasis = typename GET_PROP_TYPE(TypeTag, FeBasis);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianAssembler = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using LocalJacobian = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using NewtonMethod = typename GET_PROP_TYPE(TypeTag, NewtonMethod);

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dofCodim = dim; // TODO: for higher order FEM?

    using CoordScalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;
    using ReferenceElement = typename Dune::ReferenceElement<CoordScalar, dim>;

public:

    // copying a model is not a good idea
    FemImplicitModel(const FemImplicitModel&) = delete;

    /*!
     * \brief The constructor.
     */
    FemImplicitModel() : problemPtr_(nullptr) {}

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        // initialize the jacobian assembler first (holds the FEBasis)
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem);

        // initialize the local jacobian
        localJacobian_.init(problem);

        // apply initial solution
        uCur_.resize(asImp_().numDofs());
        asImp_().applyInitialSolution_();

        // set the solution of the "previous" time step to the initial solution.
        uPrev_ = uCur_;
    }

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param residual Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(SolutionVector &residual,
                          const SolutionVector &u)
    {
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(residual);
        curSol() = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param residual Stores the result
     */
    Scalar globalResidual(SolutionVector &residual)
    {
        residual = 0;
std::cout << "femModel globalResidual: " << std::endl;
        const auto& curSol = asImp_().curSol();
        const auto& prevSol = asImp_().prevSol();

        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();

        for (const auto& element : elements(gridView_()))
        {
            // prepare the element solutions etc...
            localView.bind(element);
            localIndexSet.bind(localView);

            auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
            ElementSolution curElemSol(numLocalDofs);
            ElementSolution prevElemSol(numLocalDofs);
            for (unsigned int i = 0; i < numLocalDofs; i++)
            {
                auto dofIdxGlobal = localIndexSet.index(i);
                curElemSol[i] = curSol[dofIdxGlobal];
                prevElemSol[i] = prevSol[dofIdxGlobal];
            }

            localResidual().eval(element, localView, localIndexSet, curElemSol, prevElemSol);

            for (unsigned int i = 0; i < localView.tree().size(); ++i)
                residual[localIndexSet.index(i)] += localResidual().residual(i);
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView_().comm().size() > 1)
            result2 = gridView_().comm().sum(result2);

        // add up the residuals on the process borders
        if (gridView_().comm().size() > 1)
        {
            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper> sumHandle(residual, vertexMapper());
            gridView_().communicate(sumHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param storage Stores the result
     */
    void globalStorage(PrimaryVariables &storage)
    {
        storage = 0;

        for (const auto& element : elements(gridView_(), Dune::Partitions::interior))
            storage += localResidual().evalStorage(element);

        if (gridView_().comm().size() > 1)
            storage = gridView_().comm().sum(storage);
    }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    const SolutionVector &curSol() const
    { return uCur_; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    const SolutionVector &prevSol() const
    { return uPrev_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianAssembler()
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    LocalJacobian &localJacobian()
    { return localJacobian_; }
    /*!
     * \copydoc localJacobian()
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }

    /*!
     * \brief Returns the local residual function.
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual()
     */
    const LocalResidual &localResidual() const
    { return localJacobian().localResidual(); }

    /*!
     * \brief Returns the maximum relative shift between two vectors of
     *        primary variables.
     *
     * \param priVars1 The first vector of primary variables
     * \param priVars2 The second vector of primary variables
     */
    Scalar relativeShiftAtDof(const PrimaryVariables &priVars1,
                              const PrimaryVariables &priVars2)
    {
        Scalar result = 0.0;
        for (int j = 0; j < numEq; ++j) {
            Scalar eqErr = std::abs(priVars1[j] - priVars2[j]);
            eqErr /= std::max<Scalar>(1.0, std::abs(priVars1[j] + priVars2[j])/2);

            result = std::max(result, eqErr);
        }
        return result;
    }

    /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver The non-linear solver
     * \param controller The controller which specifies the behaviour
     *                   of the non-linear solver
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i)
            Valgrind::CheckDefined(curSol()[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        int converged = solver.execute(controller);

        if (this->gridView_().comm().size() > 1)
        {
            converged = this->gridView_().comm().min(converged);
        }
        if (converged) {
            asImp_().updateSuccessful();
        }
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i) {
            Valgrind::CheckDefined(curSol()[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }

    /*!
     * \brief Check the plausibility of the current solution
     *
     *        This has to be done by the actual model, it knows
     *        best, what (ranges of) variables to check.
     *        This is primarily a hook
     *        which the actual model can overload.
     */
    void checkPlausibility() const {}

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primarily a hook
     *        which the actual model can overload.
     */
    void updateBegin() {}

    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateSuccessful() {}

    void newtonEndStep() {}

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        // make the current solution the previous one.
        uPrev_ = uCur_;
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.template serializeEntities<dofCodim>(asImp_(), this->gridView_());
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities<dofCodim>(asImp_(), this->gridView_());
        prevSol() = curSol();
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void serializeEntity(std::ostream &outstream,
                         const Entity &entity)
    {
        int dofIdxGlobal = vertexMapper().index(entity);

        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << dofIdxGlobal);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << curSol()[dofIdxGlobal][eqIdx] << " ";
        }
    }

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void deserializeEntity(std::istream &instream,
                           const Entity &entity)
    {
        int dofIdxGlobal = vertexMapper().index(entity);

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << dofIdxGlobal);
            instream >> curSol()[dofIdxGlobal][eqIdx];
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    std::size_t numDofs() const
    { return jacobianAssembler().numDofs(); }

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); }

    /*!
     * \brief Returns the global finite element basis
     *
     */
    const FEBasis& feBasis() const
    { return jacobianAssembler().feBasis(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        jacAsm_.template reset<JacobianAssembler>(0);
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    // TODO: WHAT IS THIS?
    // void updatePVWeights(const FVElementGeometry& fvGeometry) const {}

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param writer  The VTK multi writer object on which the fields should be added.
     * \param u       The solution function
     * \param deltaU  The delta of the solution function before and after the Newton update
     */
    template <class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
                                 const SolutionVector &u,
                                 const SolutionVector &deltaU)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        SolutionVector residual(u);
        asImp_().globalResidual(residual, u);

        // create the required scalar fields
        unsigned numDofs = asImp_().numDofs();

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            x[eqIdx] = writer.allocateManagedBuffer(numDofs);
            delta[eqIdx] = writer.allocateManagedBuffer(numDofs);
            def[eqIdx] = writer.allocateManagedBuffer(numDofs);
        }

        for (unsigned int dofIdxGlobal = 0; dofIdxGlobal < u.size(); dofIdxGlobal++)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                (*x[eqIdx])[dofIdxGlobal] = u[dofIdxGlobal][eqIdx];
                (*delta[eqIdx])[dofIdxGlobal] = - deltaU[dofIdxGlobal][eqIdx];
                (*def[eqIdx])[dofIdxGlobal] = residual[dofIdxGlobal][eqIdx];
            }
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            std::ostringstream oss;
            oss.str(""); oss << "x_" << eqIdx;
            writer.attachVertexData(*x[eqIdx], oss.str());
            oss.str(""); oss << "delta_" << eqIdx;
            writer.attachVertexData(*delta[eqIdx], oss.str());
            oss.str(""); oss << "defect_" << eqIdx;
            writer.attachVertexData(*def[eqIdx], oss.str());
        }

        asImp_().addOutputVtkFields(u, writer);
    }

    /*!
     * \brief Add the quantities of a time step which ought to be written to disk.
     *
     * This should be overwritten by the actual model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param sol The global vector of primary variable values.
     * \param writer The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {}

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * For the finite element method of order one, this is the vertex mapper
     * TODO AMG needs this method to be present, what do we do if order is > 1
     */
    template <class T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, FemBasisOrder) == 1, VertexMapper>::type &dofMapper() const
    {
        return problem_().vertexMapper();
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_().gridView(); }

protected:

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Reference to the local residal object
     */
    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        uCur_ = Scalar(0.0);

        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();

        // iterate through leaf grid and evaluate initial condition at each dof of the element
        for (const auto& element : elements(gridView_()))
        {
            // the local restrictions of the fe basis
            localView.bind(element);
            localIndexSet.bind(localView);

            // element geometry & reference element
            auto eg = element.geometry();
            const auto& refElement = ReferenceElements::general(eg.type());

            auto numLocalDofs = localView.tree().size();
            for (unsigned int i = 0; i < numLocalDofs; ++i)
            {
                const auto& localKey = localView.tree().finiteElement().localCoefficients().localKey(i);
                auto subEntity = localKey.subEntity();
                auto codim = localKey.codim();

                // get the global position of the dof
                auto globalDofPos = eg.global(refElement.position(subEntity, codim));

                // let the problem do the dirty work of nailing down the initial solution.
                uCur_[localIndexSet.index(i)] = problem_().initialAtPos(globalDofPos);
            }
        }
    }

    // the problem we want to solve. defines the constitutive relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;

    // Linearizes the problem at the current time step using the local jacobian
    std::shared_ptr<JacobianAssembler> jacAsm_;

    // cur is the current iterative solution, prev the converged solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

private:
    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
