// $Id: boxjacobian.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the jacobian of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_JACOBIAN_HH
#define DUMUX_BOX_JACOBIAN_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/boxmodels/boxscheme/fvelementgeometry.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <boost/format.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/matrix.hh>

#include "boxproperties.hh"


namespace Dumux
{
/*!
 * \ingroup BoxScheme
 * \brief Element-wise caculation of the jacobian matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxJacobian
{
private:
    typedef BoxJacobian<TypeTag>                                 ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))       Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))         Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))      GridView;

    enum {
        numEq     = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim       = GridView::dimension,
        dimWorld  = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype                CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity               Element;
    typedef typename GridView::template Codim<0>::Iterator             ElementIterator;
    typedef typename Element::EntityPointer                            ElementPointer;

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;
    typedef typename RefElemProp::ReferenceElement              ReferenceElement;

    typedef typename GridView::IntersectionIterator                    IntersectionIterator;
    typedef typename Element::Geometry                                 Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry))   FVElementGeometry;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::DofEntityMapper         DofEntityMapper;
    typedef typename SolutionTypes::VertexMapper            VertexMapper;
    typedef typename SolutionTypes::ElementMapper           ElementMapper;
    typedef typename SolutionTypes::SolutionVector          SolutionVector;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;
    typedef typename SolutionTypes::JacobianAssembler       JacobianAssembler;


    typedef typename GET_PROP(TypeTag, PTAG(VertexData))::type VertexData;
    typedef typename std::vector<VertexData>                   VertexDataArray;

    typedef std::vector<Dumux::BoundaryTypes<numEq> > BoundaryTypeArray;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq>  MatrixBlock;
    typedef Dune::Matrix<MatrixBlock>                LocalBlockMatrix;

public:
    BoxJacobian(Problem &problem)
        : problem_(problem),
          gridView_(problem.gridView()),
          curElementPtr_(* ++gridView_.template begin<0>()),
          curElementGeom_(gridView_)
    {
    }

    /*!
     * \brief Compute the global residual right hand side
     *        of an equation we would like to have zero.
     */
    void evalGlobalResidual(SolutionVector &residual)
    {
        residual = 0;
        SolutionOnElement tmpSol, tmpSolOld;
        SolutionOnElement localResid;
        localResid.resize(12);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            this->setCurrentElement(*elemIt);
            this->restrictToElement(tmpSol, this->problem().model().curSol());
            this->restrictToElement(tmpSolOld, this->problem().model().prevSol());
            this->asImp_().setCurrentSolution(tmpSol);
            this->asImp_().setPreviousSolution(tmpSolOld);

            asImp_().evalLocalResidual(localResid);

            for (int i = 0; i < elemIt->template count<dim>(); ++i) {
                int globalI = this->problem().model().vertexMapper().map(*elemIt, i, dim);
                residual[globalI] += localResid[i];
            }
        };
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element, given a local solution 'localU'.
     */
    void assemble(const Element &element)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        int numVertices = curElementGeom_.numVertices;
        SolutionOnElement localU(numVertices);
        restrictToElement(localU, problem_.model().curSol());
        asImp_().assemble_(element, localU);
    }

    /*!
     * \brief Compute the local residual, i.e. the right hand side
     *        of an equation we would like to have zero.
     */
    void evalLocalResidual(SolutionOnElement &residual)
    {
        // reset residual
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int j = 0; j < numEq; ++j)
                residual[i][j] = Scalar(0);
        }

        asImp_().evalFluxes_(residual);
        asImp_().evalVolumeTerms_(residual);

        // add the neumann fluxes
        asImp_().addNeumannFluxes_(residual);
        // set the defect of the equations used for dirichlet
        // conditions to 0
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int j = 0; j < numEq; j++) {
                if (this->bctype[i].isDirichlet(j))
                    residual[i][j] = 0;
            }
        }

#if HAVE_VALGRIND
        for (int i=0; i < curElementGeom_.numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Restrict the global function 'globalFn' to the vertices
     *        of the current element, save the result to 'dest'.
     */
    void restrictToElement(SolutionOnElement &dest,
                           const SolutionVector &globalSol) const
    {
        const DofEntityMapper &dofMapper = this->problem_.model().dofEntityMapper();
        // we assert that the i-th shape function is
        // associated to the i-th vert of the element.
        int n = curElement_().template count<dim>();
        dest.resize(n);
        for (int i = 0; i < n; i++) {
            dest[i] = globalSol[dofMapper.map(curElement_(), i, dim)];
        }
    }

    /*!
     * \brief Set the current grid element.
     */
    void setCurrentElement(const Element &element)
    {
        if (curElementPtr_ != ElementPointer(element)) {
            // update the FV element geometry if the current element
            // is to be changed
            curElementPtr_ = element;
            curElementGeom_.update(curElement_());

            // resize the array for the local jacobian matrix to the
            // current number of degrees of freedom
            int n = curElementGeom_.numVertices;
            A.setSize(n, n);
            bctype.resize(n);
            updateBoundaryTypes_();
        }
    };

    /*!
     * \brief Initialize the static data of all elements. The current
     *        solution is not yet valid when this method is
     *        called. (updateStaticData() is called as soon the
     *        current solution is valid.)
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void initStaticData()
    { }

    /*!
     * \brief Update the static data of all elements with the current solution
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void updateStaticData(SolutionVector &curSol, SolutionVector &oldSol)
    { };

    /*!
     * \brief Update the model specific vertex data of a whole
     *        element.
     */
    void updateElementData_(VertexDataArray &dest, const SolutionOnElement &sol, bool isOldSol)
    {
        int n = sol.size();
        dest.resize(n);

#ifdef ENABLE_VALGRIND
        for (int i = 0; i < dest.size(); ++i)
            Valgrind::SetUndefined(dest[i]);
#endif // ENABLE_VALGRIND

        int numVertices = this->curElement_().template count<dim>();
        for (int vertIdx = 0; vertIdx < numVertices; vertIdx++) {
            dest[vertIdx].update(sol[vertIdx],
                                 curElement_(),
                                 curFvElementGeometry(),
                                 vertIdx,
                                 problem(),
                                 isOldSol);
        }
    }


    /*!
     * \brief This returns the finite volume geometric information of the current
     *        element. <b>Use this method with great care!</b>
     */
    const FVElementGeometry &curFvElementGeometry() const
    { return curElementGeom_; }

    /*!
     * \brief Set current local solution
     */
    void setCurrentSolution(const SolutionOnElement &sol)
    {
        curElemDat_.resize(sol.size());
        asImp_().updateElementData_(curElemDat_, sol, false);
    }

    /*!
     * \brief Set local solution of the last time step
     */
    void setPreviousSolution(const SolutionOnElement &sol)
    {
        prevElemDat_.resize(sol.size());
        asImp_().updateElementData_(prevElemDat_, sol, true);
    }

    /*!
     * \brief Vary a single component of a single vert of the
     *        local solution for the current element.
     *
     * This method is a optimization, since if varying a single
     * component at a degree of freedom not the whole element cache
     * needs to be recalculated. (Updating the element cache is very
     * expensive since material laws need to be evaluated.)
     */
    void deflectCurrentSolution(SolutionOnElement &curSol,
                                int vertIdx,
                                int eqIdx,
                                Scalar value)
    {
        // stash away the orignial vertex data so that we do not need
        // to re calculate them if restoreCurSolution() is called.
        curVertexDataStash_ = curElemDat_[vertIdx];

        // recalculate the vertex data for the box which should be
        // changed
        curSol[vertIdx][eqIdx] = value;

        Valgrind::SetUndefined(curElemDat_[vertIdx]);
        curElemDat_[vertIdx].update(curSol[vertIdx],
                                    curElement_(),
                                    curFvElementGeometry(),
                                    vertIdx,
                                    problem(),
                                    false);
        Valgrind::CheckDefined(curElemDat_[vertIdx]);
    }

    /*!
     * \brief Restore the local jacobian to the state before
     *        deflectCurSolution() was called.
     *
     * This only works if deflectSolution was only called with
     * (vert, component) as arguments.
     */
    void restoreCurrentSolution(SolutionOnElement &curSol,
                                int vertIdx,
                                int eqIdx,
                                Scalar origValue)
    {
        curSol[vertIdx][eqIdx] = origValue;
        curElemDat_[vertIdx] = curVertexDataStash_;
    };

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem()
    { return problem_; };

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem() const
    { return problem_; };

    /*!
     * \brief Returns a reference to the model.
     */
    Model &model()
    { return problem_.model(); };

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model() const
    { return problem_.model(); };

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper() const
    { return model().vertexMapper(); };

    /*!
     * \brief Returns a reference to the element mapper.
     */
    const ElementMapper &elementMapper() const
    { return model().elementMapper(); };


    const MatrixBlock &mat(int i, int j) const
    { return A[i][j]; }

    const VertexDataArray& currentElementData() const
    {
        return curElemDat_;
    }

protected:
    const Element &curElement_() const
    { return *curElementPtr_; }

    // The problem we would like to solve
    Problem &problem_;

    // The grid view we are dealing with
    const GridView gridView_;

    ElementPointer    curElementPtr_;
    FVElementGeometry curElementGeom_;

    // current and previous element data. (this is model specific.)
    VertexDataArray  curElemDat_;
    VertexDataArray  prevElemDat_;

    // temporary variable to store the variable vertex data
    VertexData   curVertexDataStash_;

    void updateBoundaryTypes_()
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        int numVerts = curElement_().template count<dim>();
        for (int i = 0; i < numVerts; ++i)
            this->bctype[i].reset();

        // evaluate boundary conditions
        IntersectionIterator isIt = gridView_.template ibegin(curElement_());
        const IntersectionIterator &endIt = gridView_.template iend(curElement_());
        for (; isIt != endIt; ++isIt)
        {
            // Ignore non- boundary faces.
            if (!isIt->boundary())
                continue;

            // Set the bctype for all vertices of the face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 faceVertIdx++)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);
                // set the boundary types
                problem_.boundaryTypes(this->bctype[elemVertIdx],
                                       curElement_(),
                                       curFvElementGeometry(),
                                       *isIt,
                                       elemVertIdx,
                                       boundaryFaceIdx);
                this->bctype[elemVertIdx].checkWellPosed();
                Valgrind::CheckDefined(this->bctype[elemVertIdx]);
            }
        }
    };

    void addNeumannFluxes_(SolutionOnElement &result)
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        // evaluate boundary conditions for all intersections of
        // the current element
        IntersectionIterator isIt = gridView_.template ibegin(curElement_());
        const IntersectionIterator &endIt = gridView_.template iend(curElement_());
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
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

                if (!this->bctype[elemVertIdx].hasNeumann())
                    // the current boundary segment does not have any
                    // equation where a neumann condition should be
                    // applied
                    continue;

                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);

                // add the neuman fluxes of a single boundary segment
                addSingleNeumannSegment_(result[elemVertIdx],
                                         isIt,
                                         elemVertIdx,
                                         boundaryFaceIdx);
            }
        }
    }

    // handle boundary conditions for a single
    // sub-control volume face
    void addSingleNeumannSegment_(PrimaryVarVector &result,
                                  const IntersectionIterator &isIt,
                                  int scvIdx,
                                  int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVarVector values(0.0);

        problem_.neumann(values,
                         curElement_(),
                         curElementGeom_,
                         *isIt,
                         scvIdx,
                         boundaryFaceIdx);
        values *= curElementGeom_.boundaryFace[boundaryFaceIdx].area;
        Valgrind::CheckDefined(values);

        result += values;
    }

    void assemble_(const Element &element, SolutionOnElement& localU)
    {
        int numVertices = curElementGeom_.numVertices;

        // restrict the previous global solution to the current element
        SolutionOnElement localOldU(numVertices);
        restrictToElement(localOldU, problem_.model().prevSol());

        this->asImp_().setCurrentSolution(localU);
        this->asImp_().setPreviousSolution(localOldU);

        // approximate the local stiffness matrix numerically
        // TODO: do this analytically if possible
        SolutionOnElement partialStiffness(numVertices);
        for (int j = 0; j < numVertices; j++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().assemblePartialStiffness_(partialStiffness,
                                                   localU,
                                                   j,
                                                   pvIdx);

                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalStiffness_(j,
                                      pvIdx,
                                      partialStiffness);
            }
        }
    };

    void evalFluxes_(SolutionOnElement &residual)
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int k = 0; k < curElementGeom_.numEdges; k++)
        {
            int i = curElementGeom_.subContVolFace[k].i;
            int j = curElementGeom_.subContVolFace[k].j;

            PrimaryVarVector flux;
            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, k);
            Valgrind::CheckDefined(flux);

            // subtract fluxes from the local mass rates of the
            // respective sub control volume adjacent to the face.
            for (int eq = 0; eq < numEq; ++ eq) {
                residual[i][eq] -= flux[eq];
                residual[j][eq] += flux[eq];
            }
        }
    }

    void evalVolumeTerms_(SolutionOnElement &residual)
    {
        // evaluate the volume terms (storage + source terms)
        for (int i=0; i < curElementGeom_.numVertices; i++)
        {
            PrimaryVarVector massContrib(0), tmp(0);

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            this->asImp_().computeStorage(massContrib, i, false);
            this->asImp_().computeStorage(tmp, i, true);

            massContrib -= tmp;
            massContrib *=
                curElementGeom_.subContVol[i].volume
                /
                problem_.timeManager().timeStepSize();

            for (int j = 0; j < numEq; ++j)
                residual[i][j] += massContrib[j];

            // subtract the source term from the local rate
            PrimaryVarVector source;
            this->asImp_().computeSource(source, i);
            source *= curElementGeom_.subContVol[i].volume;

            for (int j = 0; j < numEq; ++j) {
                residual[i][j] -= source[j];

                // make sure that only defined quantities where used
                // to calculate the residual.
                Valgrind::CheckDefined(residual[i][j]);
            }
        }
    }

    /*!
     * \brief Update the stiffness matrix for all equations on all
     *        vertices of the current element with the partial
     *        derivative to the pvIdx's primary variable at the
     *        vertexIdx-th vertex.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than central differences ought to be used.
     */
    void assemblePartialStiffness_(SolutionOnElement &dest,
                                   SolutionOnElement &elemSol,
                                   int vertexIdx,
                                   int pvIdx)
    {
        Scalar eps = asImp_().numericEpsilon_(elemSol, vertexIdx, pvIdx);
        Scalar uJ = elemSol[vertexIdx][pvIdx];

        // vary the pvIdx-th primary variable at the element's j-th
        // vertex and calculate the residual, don't include the
        // boundary conditions
        asImp_().deflectCurrentSolution(elemSol, vertexIdx, pvIdx, uJ + eps);
        asImp_().evalLocalResidual(dest);
        asImp_().restoreCurrentSolution(elemSol, vertexIdx, pvIdx, uJ);

        asImp_().deflectCurrentSolution(elemSol, vertexIdx, pvIdx, uJ - eps);
        SolutionOnElement tmp(curElementGeom_.numVertices);
        asImp_().evalLocalResidual(tmp);
        asImp_().restoreCurrentSolution(elemSol, vertexIdx, pvIdx, uJ);

        // central differences
        dest -= tmp;
        dest /= 2*eps;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < dest.size(); ++i)
            Valgrind::CheckDefined(dest[i]);
#endif
    }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemSol    The current solution on the element
     * \param vertexIdx  The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon_(const SolutionOnElement &elemSol,
                           int vertIdx,
                           int pvIdx) const
    { return 1e-11*(std::abs(elemSol[vertIdx][pvIdx]) + 1); }

    /*!
     * \brief Updates the current local stiffness matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex j .
     */
    void updateLocalStiffness_(int j,
                               int pvIdx,
                               const SolutionOnElement &stiffness)
    {
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[i][j][eqIdx][pvIdx] is the approximate rate of
                // change of the residual of equation 'eqIdx' at
                // vertex 'i' depending on the primary variable
                // 'pvIdx' at vertex 'j'.
                this->A[i][j][eqIdx][pvIdx] = stiffness[i][eqIdx];
                Valgrind::CheckDefined(this->A[i][j][eqIdx][pvIdx]);
            }
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    LocalBlockMatrix  A;
    BoundaryTypeArray bctype;
};
}

#endif
