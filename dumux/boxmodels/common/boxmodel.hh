// $Id$
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

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        enableJacobianRecycling  = GET_PROP_VALUE(TypeTag, PTAG(EnableJacobianRecycling)),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) LocalResidual;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    enum { enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)) };

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

        updateBoundaryIndices_();

        int nDofs = asImp_().numDofs();
        uCur_.resize(nDofs);
        uPrev_.resize(nDofs);
        boxVolume_.resize(nDofs);

        localJacobian_.init(problem_());
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem_());

        asImp_().applyInitialSolution_();

        // also set the solution of the "previous" time step to the
        // initial solution.
        uPrev_ = uCur_;
    }

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
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(dest);
        curSol() = tmp;
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

        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            localResidual().eval(*elemIt);

            for (int i = 0; i < elemIt->template count<dim>(); ++i) {
                int globalI = vertexMapper().map(*elemIt, i, dim);
                dest[globalI] += localResidual().residual(i);
            }
        };

        Scalar result = dest.two_norm();
        /*
        Scalar result = 0;
        for (int i = 0; i < (*tmp).size(); ++i) {
            for (int j = 0; j < numEq; ++j)
                result += std::abs((*tmp)[i][j]);
        }
        */
        return result;
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

        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            localResidual().evalStorage(*elemIt);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                dest += localResidual().residual(i);
        };

        if (gridView_().comm().size() > 1)
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
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param vertIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int vertIdx, int pvIdx) const
    {
        return 1.0/std::max(std::abs(this->prevSol()[vertIdx][pvIdx]), 1.0);
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
            Scalar weight = asImp_().primaryVarWeight(vertexIdx, j);
            Scalar eqErr = std::abs(pv1[j] - pv2[j])*weight;

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

        bool converged = solver.execute(controller);
        if (converged)
            asImp_().updateSuccessful();
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
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primary a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    { }


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
        uCur_ = uPrev_;
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
        prevSol() = curSol();
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
            outstream << curSol()[vertIdx][eqIdx] << " ";
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
            instream >> curSol()[vertIdx][eqIdx];
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
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        SolutionVector globalResid(u);
        asImp_().globalResidual(globalResid, u);

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        //unsigned numElements = this->gridView_().size(0);

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i) {
            x[i] = writer.template createField<Scalar, 1>(numVertices);
            delta[i] = writer.template createField<Scalar, 1>(numVertices);
            def[i] = writer.template createField<Scalar, 1>(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int i = 0; i < numEq; ++i) {
                (*x[i])[globalIdx] = u[globalIdx][i];
                (*delta[i])[globalIdx] = - deltaU[globalIdx][i];
                (*def[i])[globalIdx] = globalResid[globalIdx][i];
            }
        }

        for (int i = 0; i < numEq; ++i) {
            writer.addVertexData(x[i], (boost::format("x_%i")%i).str().c_str());
            writer.addVertexData(delta[i], (boost::format("delta_%i")%i).str().c_str());
            writer.addVertexData(def[i], (boost::format("defect_%i")%i).str().c_str());
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
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);

        // global defect of the two auxiliary equations
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i) {
            x[i] = writer.template createField<Scalar, 1>(numVertices);
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
            writer.addVertexData(x[i], (boost::format("primaryVar%i")%i).str().c_str());
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns true if the vertex with 'globalVertIdx' is
     *        located on the grid's boundary.
     *
     * \param globalVertIdx The global index of the control volume's
     *                      associated vertex
     */
    bool onBoundary(int globalVertIdx) const
    { return boundaryIndices_.count(globalVertIdx) > 0; }

    /*!
     * \brief Returns true if a vertex is located on the grid's
     *        boundary.
     *
     * \param vertex The DUNE Codim<dim> entity associated with the
     *               control volume
     */
    bool onBoundary(const Vertex &vertex) const
    { return onBoundary(vertexMapper().map(vertex)); }

    /*!
     * \brief Returns true if a vertex is located on the grid's
     *        boundary.
     *
     * \param elem A DUNE Codim<0> entity which contains the control
     *             volume's associated vertex.
     * \param vIdx The local vertex index inside elem
     */
    bool onBoundary(const Element &elem, int vIdx) const
    { return onBoundary(vertexMapper().map(elem, vIdx, dim)); }

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
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        uCur_ = Scalar(0.0);
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
                uCur_[globalIdx] += initVal;
                Valgrind::CheckDefined(uCur_[globalIdx]);
            }
        }
        // divide all primary variables by the volume of their boxes
        int n = gridView_().size(dim);
        for (int i = 0; i < n; ++i) {
            uCur_[i] /= boxVolume(i);
        }
    }

    /*!
     * \brief Find all indices of boundary vertices.
     *
     * For this we need to loop over all intersections (which is slow
     * in general). If the DUNE grid interface would provide a
     * onBoundary() method for entities this could be done in a much
     * nicer way (actually this would not be necessary)
     */
    void updateBoundaryIndices_()
    {
        boundaryIndices_.clear();
        ElementIterator eIt = gridView_().template begin<0>();
        ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            Dune::GeometryType geoType = eIt->geometry().type();
            const ReferenceElement &refElem = ReferenceElements::general(geoType);

            IntersectionIterator isIt = gridView_().ibegin(*eIt);
            IntersectionIterator isEndIt = gridView_().iend(*eIt);
            for (; isIt != isEndIt; ++isIt) {
                if (!isIt->boundary())
                    continue;
                // add all vertices on the intersection to the set of
                // boundary vertices
                int faceIdx = isIt->indexInInside();
                int numFaceVerts = refElem.size(faceIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    int elemVertIdx = refElem.subEntity(faceIdx,
                                                        1,
                                                        faceVertIdx,
                                                        dim);
                    int globalVertIdx = vertexMapper().map(*eIt, elemVertIdx, dim);
                    boundaryIndices_.insert(globalVertIdx);
                }
            }
        }
    }

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

    // the set of all indices of vertices on the boundary
    std::set<int> boundaryIndices_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;
};
}

#endif
