/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for models using box discretization
 */
#ifndef DUMUX_BOX_MODEL_HH
#define DUMUX_BOX_MODEL_HH

#include "boxproperties.hh"
#include "boxpropertydefaults.hh"

#include "boxelementvolumevariables.hh"
#include "boxlocaljacobian.hh"
#include "boxlocalresidual.hh"

#include <dumux/parallel/vertexhandles.hh>


namespace Dumux
{

/*!
 * \defgroup BoxModel The core infrastructure of the box scheme
 */
/*!
 * \ingroup BoxModel
 * \defgroup BoxModels Physical models and problems which use the box scheme
 */


/*!
 * \ingroup BoxModel
 *
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class BoxModel
{
    typedef BoxModel<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DofMapper)) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianAssembler)) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVariables)) ElementVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;


    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        historySize = GET_PROP_VALUE(TypeTag, PTAG(TimeDiscHistorySize)),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) LocalResidual;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    // copying a model is not a good idea
    BoxModel(const BoxModel &);

public:
    /*!
     * \brief The constructor.
     */
    BoxModel()
    { }

    ~BoxModel()
    { delete jacAsm_;  }

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param prob The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &prob)
    {
        problemPtr_ = &prob;

        updateBoundaryTypes_();

        int nDofs = asImp_().numDofs();
        for (int historyIdx = 0; historyIdx < historySize; ++historyIdx)
            solution_[historyIdx].resize(nDofs);
        boxVolume_.resize(nDofs);

        localJacobian_.init(problem_());
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem_());

        asImp_().applyInitialSolution_();

        // resize the hint vectors
        if (enableHints_()) {
            int nVerts = gridView_().size(dim);
            for (int historyIdx = 0; historyIdx < historySize; ++historyIdx)
                hints_[historyIdx].resize(nVerts);
            hintsUsable_.resize(nVerts);
            std::fill(hintsUsable_.begin(),
                      hintsUsable_.end(),
                      false);
        }

        // also set the solution of the "previous" time step to the
        // initial solution.
        solution_[/*historyIdx=*/1] = solution_[/*historyIdx=*/0];
    }

    void setAllHints(ElementVariables &elemVars) const
    {
        if (!enableHints_())
            return;

        for (int historyIdx = 0; historyIdx < historySize; ++historyIdx)
            setHints(elemVars, historyIdx);
    }

    void setHints(ElementVariables &elemVars, int historyIdx) const
    {
        if (!enableHints_())
            return;

        int numScv = elemVars.numScv();
        for (int scvIdx = 0; scvIdx < numScv; ++scvIdx) {
            setHintsScv(elemVars, historyIdx, scvIdx);
        }
    }

    void setHintsScv(ElementVariables &elemVars, int historyIdx, int scvIdx) const
    {
        int globalIdx = vertexMapper().map(elemVars.element(), scvIdx, dim);
        
        if (!hintsUsable_[globalIdx])
            elemVars.volVars(scvIdx, historyIdx).setHint(NULL);
        else
            elemVars.volVars(scvIdx, historyIdx)
                .setHint(&hints_[historyIdx][globalIdx]);
    };

    void updateAllHints(const ElementVariables &elemVars) const
    {
        if (!enableHints_())
            return;

        for (int historyIdx = 0; historyIdx < historySize; ++historyIdx) {
            updateHints(elemVars, historyIdx);
        }
    };

    void updateHints(const ElementVariables &elemVars, int historyIdx) const
    {
        if (!enableHints_())
            return;
        
        for (int scvIdx = 0; scvIdx < elemVars.numScv(); ++scvIdx) {
            int globalIdx = vertexMapper().map(elemVars.element(), scvIdx, dim);
            hints_[historyIdx][globalIdx] = elemVars.volVars(scvIdx, historyIdx);
            hintsUsable_[globalIdx] = true;
        }
    };

    void shiftHints()
    {
        for (int historyIdx = 0; historyIdx < historySize - 1; ++ historyIdx)
            hints_[historyIdx + 1] = hints_[historyIdx];
    };

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param dest Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(SolutionVector &dest,
                          const SolutionVector &u)
    {
        SolutionVector tmp(solution(/*historyIdx=*/0));
        solution(/*historyIdx=*/0) = u;
        Scalar res = globalResidual(dest);
        solution(/*historyIdx=*/0) = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param dest Stores the result
     */
    Scalar globalResidual(SolutionVector &dest)
    {
        dest = 0;

        ElementVariables elemVars(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemVars.updateAll(*elemIt);
            localResidual().eval(elemVars);

            for (int scvIdx = 0; scvIdx < elemVars.numScv(); ++scvIdx) {
                int globalI = vertexMapper().map(*elemIt, scvIdx, dim);
                dest[globalI] += localResidual().residual(scvIdx);
            }
        };

        // calculate the square norm of the residual
        Scalar result2 = dest.two_norm2();
        result2 = gridView().comm().sum(result2);

        // add up the residuals on the process borders
        if (gridView().comm().size() > 1) {
            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumHandle(dest, vertexMapper());
            gridView().communicate(sumHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param dest Stores the result
     */
    void globalStorage(PrimaryVariables &dest)
    {
        dest = 0;
        
        ElementVariables elemVars(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemVars.updateFVElemGeom(*elemIt);
            elemVars.updateScvVars(/*historyIdx=*/0);
            localResidual().evalStorage(elemVars, /*historyIdx=*/0);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                dest += localResidual().storageTerm()[i];
        };

        dest = gridView_().comm().sum(dest);
    }

    /*!
     * \brief Returns the volume \f$\mathrm{[m^3]}\f$ of a given control volume.
     *
     * \param globalIdx The global index of the control volume's
     *                  associated vertex
     */
    Scalar boxVolume(int globalIdx) const
    { return boxVolume_[globalIdx][0]; }

    /*!
     * \brief Reference to the solution at a given history index as a block vector.
     */
    const SolutionVector &solution(int historyIdx) const
    { return solution_[historyIdx]; }

    /*!
     * \brief Reference to the solution at a given history index as a block vector.
     */
    SolutionVector &solution(int historyIdx)
    { return solution_[historyIdx]; }

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
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param vertIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int vertIdx, int pvIdx) const
    {
        Scalar absPv = std::abs(this->solution(/*historyIdx=*/1)[vertIdx][pvIdx]);
        return 1.0/std::max(absPv, 1.0);
    }

    /*!
     * \brief Returns the relative error between two vectors of
     *        primary variables.
     *
     * \param vertexIdx The global index of the control volume's
     *                  associated vertex
     * \param pv1 The first vector of primary variables
     * \param pv2 The second vector of primary variables
     *
     * \todo The vertexIdx argument is pretty hacky. it is required by
     *       models with pseudo primary variables (i.e. the primary
     *       variable switching models). the clean solution would be
     *       to access the pseudo primary variables from the primary
     *       variables.
     */
    Scalar relativeErrorVertex(int vertexIdx,
                               const PrimaryVariables &pv1,
                               const PrimaryVariables &pv2)
    {
        Scalar result = 0.0;
        for (int j = 0; j < numEq; ++j) {
            //Scalar weight = asImp_().primaryVarWeight(vertexIdx, j);
            //Scalar eqErr = std::abs(pv1[j] - pv2[j])*weight;
            Scalar eqErr = std::abs(pv1[j] - pv2[j]);
            eqErr /= std::max<Scalar>(1.0, std::abs(pv1[j] + pv2[j])/2);

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
        for (size_t i = 0; i < solution(/*historyIdx=*/0).size(); ++i)
            Valgrind::CheckDefined(solution(/*historyIdx=*/0)[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.execute(controller);
        if (converged) {
            asImp_().updateSuccessful();
        }
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < solution(/*historyIdx=*/0).size(); ++i) {
            Valgrind::CheckDefined(solution(/*historyIdx=*/0)[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }


    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primary a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    {
        updateBoundaryTypes_();
    }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    { };

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        hints_[/*historyIdx=*/0] = hints_[/*historyIdx=*/1];

        jacAsm_->reassembleAll();
    };

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
        solution_[/*historyIdx=*/1] = solution_[/*historyIdx=*/0];

        // shift the hints by one position in the history
        shiftHints();
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
    { res.template serializeEntities<dim>(asImp_(), this->gridView_()); }

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
        res.template deserializeEntities<dim>(asImp_(), this->gridView_());
        solution_[/*historyIdx=*/1] = solution_[/*historyIdx=*/0];
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param vert The DUNE Codim<dim> entity which's data should be
     *             serialized
     */
    void serializeEntity(std::ostream &outstream,
                         const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);

        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << solution_[/*historyIdx=*/0][vertIdx][eqIdx] << " ";
        }
    };

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param vert The DUNE Codim<dim> entity which's data should be
     *             deserialized
     */
    void deserializeEntity(std::istream &instream,
                           const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << vertIdx);
            instream >> solution_[/*historyIdx=*/0][vertIdx][eqIdx];
        }
    };

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    { return gridView_().size(dim); }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * This usually means a mapper for vertices.
     */
    const DofMapper &dofMapper() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); };

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler;
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Return a pointer to the BoundaryTypes for a given global
     *        vertex index or 0 if the vertex is not on the boundary.
     */
    const BoundaryTypes *boundaryTypes(int globalIdx) const
    {
        static BoundaryTypes dummy;
        int bvertIdx = boundaryVertexIndex_[globalIdx];
        if (bvertIdx < 0)
            return &dummy;
        return &boundaryTypes_[bvertIdx];
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    void updatePVWeights(const ElementVariables &elemVars) const
    { };

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

        SolutionVector globalResid(u);
        asImp_().globalResidual(globalResid, u);

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        //unsigned numElements = this->gridView_().size(0);

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        ScalarField* relError = writer.allocateManagedBuffer(numVertices);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            x[pvIdx] = writer.allocateManagedBuffer(numVertices);
            delta[pvIdx] = writer.allocateManagedBuffer(numVertices);
            def[pvIdx] = writer.allocateManagedBuffer(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                (*x[pvIdx])[globalIdx] = u[globalIdx][pvIdx];
                (*delta[pvIdx])[globalIdx] = 
                    - deltaU[globalIdx][pvIdx];
                (*def[pvIdx])[globalIdx] = globalResid[globalIdx][pvIdx];
            }

            PrimaryVariables uOld(u[globalIdx]);
            PrimaryVariables uNew(uOld - deltaU[globalIdx]);
            (*relError)[globalIdx] = asImp_().relativeErrorVertex(globalIdx, 
                                                                  uOld,
                                                                  uNew);
        }

        writer.attachVertexData(*relError, "relative Error");
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            writer.attachVertexData(*x[pvIdx], (boost::format("x_%i")%pvIdx).str().c_str());
            writer.attachVertexData(*delta[pvIdx], (boost::format("delta_%i")%pvIdx).str().c_str());
            writer.attachVertexData(*def[pvIdx], (boost::format("defect_%i")%pvIdx).str().c_str());
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
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);

        // global defect of the two auxiliary equations
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i) {
            x[i] = writer.allocateManagedBuffer(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int i = 0; i < numEq; ++i) {
                (*x[i])[globalIdx] = sol[globalIdx][i];
            }
        }

        for (int i = 0; i < numEq; ++i)
            writer.attachVertexData(*x[i], (boost::format("primaryVar%i")%i).str().c_str());
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_().gridView(); }

protected:
    static bool enableHints_()
    { return GET_PARAM(TypeTag, bool, EnableHints); }

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
     * \brief Updates the stuff which determines a vertex' or
     *        element's boundary type
     */
    void updateBoundaryTypes_()
    {
        // resize the vectors
        boundaryVertexIndex_.resize(numDofs());
        std::fill(boundaryVertexIndex_.begin(),
                  boundaryVertexIndex_.end(),
                  -1);

        int numBoundaryVertices = 0;
        ElementVariables elemVars(problem_());

        // loop over all elements of the grid
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            // do nothing if the element does not have boundary intersections
            if (!elemIt->hasBoundaryIntersections())
                continue;

            // retrieve the reference element for the current element
            const Element &elem = *elemIt;
            Dune::GeometryType geoType = elem.geometry().type();
            const ReferenceElement &refElem = ReferenceElements::general(geoType);
            
            elemVars.updateFVElemGeom(*elemIt);

            // loop over all intersections of the element
            IntersectionIterator isIt = gridView_().ibegin(elem);
            const IntersectionIterator &endIt = gridView_().iend(elem);
            for (; isIt != endIt; ++isIt)
            {
                // do nothing if the face is _not_ on the boundary
                if (!isIt->boundary())
                    continue;
                
                // loop over all vertices of the intersection
                int faceIdx = isIt->indexInInside();
                int numFaceVerts = refElem.size(faceIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    // find the local element index of the face's
                    // vertex
                    int scvIdx = refElem.subEntity(/*entityIdx=*/faceIdx,
                                                   /*entityCodim=*/1,
                                                   /*subEntityIdx=*/faceVertIdx,
                                                   /*subEntityCodim=*/dim);
                    int globalIdx = vertexMapper().map(*elemIt, scvIdx, /*codim=*/dim);
                    if (boundaryVertexIndex_[globalIdx] >= 0)
                        continue; // vertex has already been visited
                    
                    // add a BoundaryTypes object
                    if (boundaryTypes_.size() <= numBoundaryVertices)
                        boundaryTypes_.resize(numBoundaryVertices + 1);
                    BoundaryTypes &bTypes = boundaryTypes_[numBoundaryVertices];
                    boundaryVertexIndex_[globalIdx] = numBoundaryVertices;
                    ++numBoundaryVertices;

                    problem_().boundaryTypes(bTypes, elemVars, scvIdx);
                } // loop over intersection's vertices
            } // loop over intersections
        } // loop over elements
    };

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        SolutionVector &uCur = solution(/*historyIdx=*/0);
        uCur = Scalar(0.0);
        boxVolume_ = Scalar(0.0);

        FVElementGeometry fvElemGeom;

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator it = gridView_().template begin<0>();
        const ElementIterator &eendit = gridView_().template end<0>();
        for (; it != eendit; ++it) {
            // deal with the current element
            fvElemGeom.update(gridView_(), *it);

            // loop over all element vertices, i.e. sub control volumes
            for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; scvIdx++)
            {
                // map the local vertex index to the global one
                int globalIdx = vertexMapper().map(*it,
                                                   scvIdx,
                                                   dim);

                // let the problem do the dirty work of nailing down
                // the initial solution.
                PrimaryVariables initVal;
                Valgrind::SetUndefined(initVal);
                problem_().initial(initVal,
                                   *it,
                                   fvElemGeom,
                                   scvIdx);
                Valgrind::CheckDefined(initVal);

                // add up the initial values of all sub-control
                // volumes. If the initial values disagree for
                // different sub control volumes, the initial value
                // will be the arithmetic mean.
                initVal *= fvElemGeom.subContVol[scvIdx].volume;
                boxVolume_[globalIdx] += fvElemGeom.subContVol[scvIdx].volume;
                uCur[globalIdx] += initVal;
                Valgrind::CheckDefined(uCur[globalIdx]);
            }
        }

        // add up the primary variables and the volumes of the boxes
        // which cross process borders
        if (gridView().comm().size() > 1) {
            VertexHandleSum<Dune::FieldVector<Scalar, 1>,
                Dune::BlockVector<Dune::FieldVector<Scalar, 1> >,
                VertexMapper> sumVolumeHandle(boxVolume_, vertexMapper());
            gridView().communicate(sumVolumeHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);

            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumPVHandle(uCur, vertexMapper());
            gridView().communicate(sumPVHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

        // divide all primary variables by the volume of their boxes
        int n = gridView_().size(dim);
        for (int i = 0; i < n; ++i) {
            uCur[i] /= boxVolume(i);
        }
    }

    // the hint cache for the previous and the current volume
    // variables
    mutable std::vector<bool> hintsUsable_;
    mutable std::vector<VolumeVariables> hints_[historySize];

    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; };

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector solution_[historySize];

    // all the index of the BoundaryTypes object for a vertex
    std::vector<int> boundaryVertexIndex_;
    std::vector<BoundaryTypes> boundaryTypes_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;
};
}

#endif
