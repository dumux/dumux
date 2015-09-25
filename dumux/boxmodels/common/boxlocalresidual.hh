// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
 * \brief Caculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>
#include <dune/grid/common/geometry.hh>

#include <dumux/common/valgrind.hh>

#include "boxproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef BoxLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>  LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq>  MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &);

public:
    BoxLocalResidual()
    { }

    ~BoxLocalResidual()
    { }

    /*!
     * \brief Initialize the local residual.
     *
     * This assumes that all objects of the simulation have been fully
     * allocated but not necessarrily initialized completely.
     *
     * \param prob The representation of the physical problem to be
     *             solved.
     */
    void init(Problem &prob)
    { problemPtr_ = &prob; }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    void eval(const Element &element)
    {
        FVElementGeometry fvElemGeom;
        fvElemGeom.update(gridView_(), element);
        fvElemGeomPtr_ = &fvElemGeom;

        ElementVolumeVariables volVarsPrev, volVarsCur;
        volVarsPrev.update(problem_(),
                           element,
                           fvElemGeom_(),
                           true /* oldSol? */);
        volVarsCur.update(problem_(),
                          element,
                          fvElemGeom_(),
                          false /* oldSol? */);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvElemGeom_());

        // this is pretty much a HACK because the internal state of
        // the problem is not supposed to be changed during the
        // evaluation of the residual. (Reasons: It is a violation of
        // abstraction, makes everything more prone to errors and is
        // not thread save.) The real solution are context objects!
        problem_().updateCouplingParams(element);

        asImp_().eval(element, fvElemGeom_(), volVarsPrev, volVarsCur, bcTypes);
    }
    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero. Gets a solution vector computed by PDELab
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param elementSolVector The local solution for the element using PDELab ordering
     */
    template<typename ElemSolVectorType>
    void evalPDELab(const Element &element, const ElemSolVectorType& elementSolVector)
    {
        FVElementGeometry fvGeom;
        fvGeom.update(gridView_(), element);
        ElementVolumeVariables volVarsPrev, volVarsCur;
        volVarsPrev.update(problem_(),
                           element,
                           fvGeom,
                           true /* oldSol? */);
        volVarsCur.updatePDELab(problem_(),
                          element,
                          fvGeom,
                          elementSolVector);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeom);

        // this is pretty much a HACK because the internal state of
        // the problem is not supposed to be changed during the
        // evaluation of the residual. (Reasons: It is a violation of
        // abstraction, makes everything more prone to errors and is
        // not thread save.) The real solution are context objects!
        problem_().updateCouplingParams(element);

        asImp_().eval(element, fvGeom, volVarsPrev, volVarsCur, bcTypes);
    }

    /*!
     * \brief Compute the storage term for the current solution.
     *
     * This can be used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param element The DUNE Codim<0> entity for which the storage
     *                term ought to be calculated
     */
    void evalStorage(const Element &element)
    {
        elemPtr_ = &element;

        FVElementGeometry fvElemGeom;
        fvElemGeom.update(gridView_(), element);
        fvElemGeomPtr_ = &fvElemGeom;

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvElemGeom_());
        bcTypesPtr_ = &bcTypes;

        // no previous volume variables!
        prevVolVarsPtr_ = 0;

        ElementVolumeVariables volVars;
        volVars.update(problem_(), element, fvElemGeom_(), false);
        curVolVarsPtr_ = &volVars;

        residual_.resize(fvElemGeom_().numVertices);
        residual_ = 0;

        asImp_().evalStorage_();
    }

    /*!
     * \brief Compute the flux term for the current solution.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param curVolVars The volume averaged variables for all
     *                   sub-contol volumes of the element
     */
    void evalFluxes(const Element &element,
                    const ElementVolumeVariables &curVolVars)
    {
        elemPtr_ = &element;

        FVElementGeometry fvElemGeom;
        fvElemGeom.update(gridView_(), element);
        fvElemGeomPtr_ = &fvElemGeom;

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvElemGeom_());

        residual_.resize(fvElemGeom_().numVertices);
        residual_ = 0;

        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = 0;
        curVolVarsPtr_ = &curVolVars;
        asImp_().evalFluxes_();
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeom The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                   sub-contol volumes of the element at the previous
     *                   time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-contol volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    void eval(const Element &element,
              const FVElementGeometry &fvGeom,
              const ElementVolumeVariables &prevVolVars,
              const ElementVolumeVariables &curVolVars,
              const ElementBoundaryTypes &bcTypes)
    {
        //Valgrind::CheckDefined(fvGeom);
        Valgrind::CheckDefined(prevVolVars);
        Valgrind::CheckDefined(curVolVars);

#if HAVE_VALGRIND
        for (int i=0; i < fvGeom.numVertices; i++) {
            Valgrind::CheckDefined(prevVolVars[i]);
            Valgrind::CheckDefined(curVolVars[i]);
        }
#endif // HAVE_VALGRIND

        elemPtr_ = &element;
        fvElemGeomPtr_ = &fvGeom;
        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = &prevVolVars;
        curVolVarsPtr_ = &curVolVars;

        // reset residual
        int numVerts = fvElemGeom_().numVertices;
        residual_.resize(numVerts);
        residual_ = 0;

        asImp_().evalFluxes_();
        asImp_().evalVolumeTerms_();

        // evaluate the boundary
        asImp_().evalBoundary_();

#if HAVE_VALGRIND
        for (int i=0; i < fvElemGeom_().numVertices; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element.
     */
    const ElementSolutionVector &residual() const
    { return residual_; }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element.
     *
     * \param scvIdx The local index of the sub-control volume
     *               (i.e. the element's local vertex index)
     */
    const PrimaryVariables &residual(int scvIdx) const
    { return residual_[scvIdx]; }

    /*!
     * \brief Returns the local residual for a given sub-control
     *        volume of the element - as a reference!
     *
     * \param scvIdx The local index of the sub-control volume
     *               (i.e. the element's local vertex index)
     */
    PrimaryVariables& residualReference(int scvIdx)
    { return residual_[scvIdx]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    /*!
     * \brief Evaluate all boundary conditions on the current element.
     */
    void evalBoundary_()
    {
        if (bcTypes_().hasNeumann())
            asImp_().evalNeumann_();
        if (bcTypes_().hasDirichlet())
            asImp_().evalDirichlet_();
    }

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_()
    {
        PrimaryVariables tmp;
        for (int i = 0; i < fvElemGeom_().numVertices; ++i) {
            const BoundaryTypes &bcTypes = bcTypes_(i);
            if (! bcTypes.hasDirichlet())
                continue;

            // ask the problem for the dirichlet values
            const VertexPointer vPtr = elem_().template subEntity<dim>(i);
            asImp_().problem_().dirichlet(tmp, *vPtr);

            // set the dirichlet conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isDirichlet(eqIdx))
                    continue;
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                this->residual_[i][eqIdx] =
                    curPrimaryVars_(i)[pvIdx] - tmp[pvIdx];
            };
        };
    }

    /*!
     * \brief Add all Neumann boundary conditions to the local
     *        residual.
     */
    void evalNeumann_()
    {
        Dune::GeometryType geoType = elem_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        IntersectionIterator isIt = gridView_().ibegin(elem_());
        const IntersectionIterator &endIt = gridView_().iend(elem_());
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

                int boundaryFaceIdx =
                    fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                evalNeumannSegment_(isIt,
                                    elemVertIdx,
                                    boundaryFaceIdx);
            }
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             int scvIdx,
                             int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables values(0.0);
        const BoundaryTypes &bcTypes = bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            problem_().neumann(values,
                               elem_(),
                               fvElemGeom_(),
                               *isIt,
                               scvIdx,
                               boundaryFaceIdx);
            values *= fvElemGeom_().boundaryFace[boundaryFaceIdx].area;
            Valgrind::CheckDefined(values);
            residual_[scvIdx] += values;
        }
    }

    /*!
     * \brief Add the flux terms to the local residual of all
     *        sub-control volumes of the current element.
     */
    void evalFluxes_()
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int k = 0; k < fvElemGeom_().numEdges; k++)
        {
            int i = fvElemGeom_().subContVolFace[k].i;
            int j = fvElemGeom_().subContVolFace[k].j;

            PrimaryVariables flux;
            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, k);
            Valgrind::CheckDefined(flux);

            // The balance equation for a finite volume is:
            //
            // dStorage/dt = Flux + Source
            //
            // where the 'Flux' and the 'Source' terms represent the
            // mass per second which _ENTER_ the finite
            // volume. Re-arranging this, we get
            //
            // dStorage/dt - Source - Flux = 0
            //
            // Since the flux calculated by computeFlux() goes _OUT_
            // of sub-control volume i and _INTO_ sub-control volume
            // j, we need to add the flux to finite volume i and
            // subtract it from finite volume j
            residual_[i] += flux;
            residual_[j] -= flux;
        }
    }

    /*!
     * \brief Set the local residual to the storage terms of all
     *        sub-control volumes of the current element.
     */
    void evalStorage_()
    {
        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int i=0; i < fvElemGeom_().numVertices; i++) {
            residual_[i] = 0.0;
            this->asImp_().computeStorage(residual_[i], i, false /*isOldSol*/);
            residual_[i] *= fvElemGeom_().subContVol[i].volume;
        }
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    void evalVolumeTerms_()
    {
        // evaluate the volume terms (storage + source terms)
        for (int i=0; i < fvElemGeom_().numVertices; i++)
        {
            PrimaryVariables dStorage_dt(0), tmp(0);

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            this->asImp_().computeStorage(dStorage_dt, i, false);
            this->asImp_().computeStorage(tmp, i, true);

            dStorage_dt -= tmp;
            dStorage_dt *=
                fvElemGeom_().subContVol[i].volume
                /
                problem_().timeManager().timeStepSize();
            residual_[i] += dStorage_dt;

            // subtract the source term from the local rate
            PrimaryVariables source(0.0);
            this->asImp_().computeSource(source, i);
            source *= fvElemGeom_().subContVol[i].volume;
            residual_[i] -= source;

            // make sure that only defined quantities where used
            // to calculate the residual.
            Valgrind::CheckDefined(residual_[i]);
        }
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return *problemPtr_; };

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); };

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the current element.
     */
    const Element &elem_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }

    /*!
     * \brief Returns a reference to the current element's finite
     *        volume geometry.
     */
    const FVElementGeometry &fvElemGeom_() const
    {
        Valgrind::CheckDefined(fvElemGeomPtr_);
        return *fvElemGeomPtr_;
    }

    /*!
     * \brief Returns a reference to the primary variables of the i'th
     *        sub-control volume of the current element.
     */
    const PrimaryVariables &curPrimaryVars_(int i) const
    {
        return curVolVars_(i).primaryVars();
    }

    /*!
     * \brief Returns a reference to the current volume variables of
     *        all sub-control volumes of the current element.
     */
    const ElementVolumeVariables &curVolVars_() const
    {
        Valgrind::CheckDefined(curVolVarsPtr_);
        return *curVolVarsPtr_;
    }

    /*!
     * \brief Returns a reference to the volume variables of the i-th
     *        sub-control volume of the current element.
     */
    const VolumeVariables &curVolVars_(int i) const
    {
        return curVolVars_()[i];
    }

    /*!
     * \brief Returns a reference to the previous time step's volume
     *        variables of all sub-control volumes of the current
     *        element.
     */
    const ElementVolumeVariables &prevVolVars_() const
    {
        Valgrind::CheckDefined(prevVolVarsPtr_);
        return *prevVolVarsPtr_;
    }

    /*!
     * \brief Returns a reference to the previous time step's volume
     *        variables of the i-th sub-control volume of the current
     *        element.
     */
    const VolumeVariables &prevVolVars_(int i) const
    {
        return prevVolVars_()[i];
    }

    /*!
     * \brief Returns a reference to the boundary types of all
     *        sub-control volumes of the current element.
     */
    const ElementBoundaryTypes &bcTypes_() const
    {
        Valgrind::CheckDefined(bcTypesPtr_);
        return *bcTypesPtr_;
    }

    /*!
     * \brief Returns a reference to the boundary types of the i-th
     *        sub-control volume of the current element.
     */
    const BoundaryTypes &bcTypes_(int i) const
    {
        return bcTypes_()[i];
    }

protected:
    ElementSolutionVector residual_;

    // The problem we would like to solve
    Problem *problemPtr_;

    const Element *elemPtr_;
    const FVElementGeometry *fvElemGeomPtr_;

    // current and previous secondary variables for the element
    const ElementVolumeVariables *prevVolVarsPtr_;
    const ElementVolumeVariables *curVolVarsPtr_;

    const ElementBoundaryTypes *bcTypesPtr_;
};

}

#endif
