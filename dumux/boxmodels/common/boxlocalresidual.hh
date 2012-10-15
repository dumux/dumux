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
 * \brief Calculates the residual of models based on the box scheme element-wise.
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
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the residual matrix for models
 *        based on the box scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

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
     * allocated but not necessarily initialized completely.
     *
     * \param problem The representation of the physical problem to be
     *             solved.
     */
    void init(Problem &problem)
    { problemPtr_ = &problem; }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    void eval(const Element &element)
    {
        FVElementGeometry fvGeometry;

        fvGeometry.update(gridView_(), element);
        fvElemGeomPtr_ = &fvGeometry;

        ElementVolumeVariables volVarsPrev, volVarsCur;
        // update the hints
        model_().setHints(element, volVarsPrev, volVarsCur);

        volVarsPrev.update(problem_(),
                           element,
                           fvGeometry_(),
                           true /* oldSol? */);
        volVarsCur.update(problem_(),
                          element,
                          fvGeometry_(),
                          false /* oldSol? */);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());

        // this is pretty much a HACK because the internal state of
        // the problem is not supposed to be changed during the
        // evaluation of the residual. (Reasons: It is a violation of
        // abstraction, makes everything more prone to errors and is
        // not thread save.) The real solution are context objects!
        problem_().updateCouplingParams(element);

        asImp_().eval(element, fvGeometry_(), volVarsPrev, volVarsCur, bcTypes);
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

        FVElementGeometry fvGeometry;
        fvGeometry.update(gridView_(), element);
        fvElemGeomPtr_ = &fvGeometry;

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());
        bcTypesPtr_ = &bcTypes;

        // no previous volume variables!
        prevVolVarsPtr_ = 0;

        ElementVolumeVariables volVars;

        // update the hints
        model_().setHints(element, volVars);

        // calculate volume current variables
        volVars.update(problem_(), element, fvGeometry_(), false);
        curVolVarsPtr_ = &volVars;

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

        FVElementGeometry fvGeometry;
        fvGeometry.update(gridView_(), element);
        fvElemGeomPtr_ = &fvGeometry;

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());

        residual_.resize(fvGeometry_().numVertices);
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
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the previous
     *                   time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    void eval(const Element &element,
              const FVElementGeometry &fvGeometry,
              const ElementVolumeVariables &prevVolVars,
              const ElementVolumeVariables &curVolVars,
              const ElementBoundaryTypes &bcTypes)
    {
        Valgrind::CheckDefined(prevVolVars);
        Valgrind::CheckDefined(curVolVars);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry.numVertices; i++) {
            prevVolVars[i].checkDefined();
            curVolVars[i].checkDefined();
        }
#endif // HAVE_VALGRIND

        elemPtr_ = &element;
        fvElemGeomPtr_ = &fvGeometry;
        bcTypesPtr_ = &bcTypes;
        prevVolVarsPtr_ = &prevVolVars;
        curVolVarsPtr_ = &curVolVars;

        // resize the vectors for all terms
        int numVerts = fvGeometry_().numVertices;
        residual_.resize(numVerts);
        storageTerm_.resize(numVerts);

        residual_ = 0.0;
        storageTerm_ = 0.0;

        asImp_().evalFluxes_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numVertices; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND

        asImp_().evalVolumeTerms_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numVertices; i++) {
            Valgrind::CheckDefined(residual_[i]);
        }
#endif // HAVE_VALGRIND

        // evaluate the boundary conditions
        asImp_().evalBoundary_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numVertices; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Returns the local residual for all sub-control
     *        volumes of the element.
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
    const PrimaryVariables &residual(const int scvIdx) const
    { return residual_[scvIdx]; }

    /*!
     * \brief Returns the storage term for all sub-control volumes of the
     *        element.
     */
    const ElementSolutionVector &storageTerm() const
    { return storageTerm_; }

    /*!
     * \brief Returns the storage term for a given sub-control volumes
     *        of the element.
     */
    const PrimaryVariables &storageTerm(const int scvIdx) const
    { return storageTerm_[scvIdx]; }

protected:
    Implementation &asImp_()
    {
        assert(static_cast<Implementation*>(this) != 0);
        return *static_cast<Implementation*>(this);
    }

    const Implementation &asImp_() const
    {
        assert(static_cast<const Implementation*>(this) != 0);
        return *static_cast<const Implementation*>(this);
    }

    /*!
     * \brief Evaluate the boundary conditions
     *        of the current element.
     */
    void evalBoundary_()
    {
        if (bcTypes_().hasNeumann() || bcTypes_().hasOutflow())
            asImp_().evalBoundaryFluxes_();
#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numVertices; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND

        if (bcTypes_().hasDirichlet())
            asImp_().evalDirichlet_();
    }

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_()
    {
        PrimaryVariables dirichletValues(0);
        for (int scvIdx = 0; scvIdx < fvGeometry_().numVertices; ++scvIdx) {
            const BoundaryTypes &bcTypes = bcTypes_(scvIdx);
            if (! bcTypes.hasDirichlet())
                continue;

            // ask the problem for the dirichlet values
            const VertexPointer vPtr = element_().template subEntity<dim>(scvIdx);
            Valgrind::SetUndefined(dirichletValues);
            asImp_().problem_().dirichlet(dirichletValues, *vPtr);

            // set the dirichlet conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isDirichlet(eqIdx))
                    continue;
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(dirichletValues[pvIdx]);

                residual_[scvIdx][eqIdx] =
                        curPriVar_(scvIdx, pvIdx) - dirichletValues[pvIdx];

                storageTerm_[scvIdx][eqIdx] = 0.0;
            }
        }
    }

    /*!
     * \brief Add all Neumann and outflow boundary conditions to the local
     *        residual.
     */
    void evalBoundaryFluxes_()
    {
        Dune::GeometryType geoType = element_().geometry().type();
        const ReferenceElement &refElement = ReferenceElements::general(geoType);

        IntersectionIterator isIt = gridView_().ibegin(element_());
        const IntersectionIterator &endIt = gridView_().iend(element_());
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElement.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int scvIdx = refElement.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);

                int boundaryFaceIdx =
                    fvGeometry_().boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                asImp_().evalNeumannSegment_(isIt,
                                             scvIdx,
                                             boundaryFaceIdx);
                // evaluate the outflow conditions at the boundary face
                // ATTENTION: This is so far a beta version that is only for the 2p2c and 2p2cni model
                //              available and not thoroughly tested.
                asImp_().evalOutflowSegment_(isIt,
                                             scvIdx,
                                             boundaryFaceIdx);
            }
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables neumannFlux(0.0);
        const BoundaryTypes &bcTypes = bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            Valgrind::SetUndefined(neumannFlux);
            problem_().boxSDNeumann(neumannFlux,
                                    element_(),
                                    fvGeometry_(),
                                    *isIt,
                                    scvIdx,
                                    boundaryFaceIdx,
                                    curVolVars_());
            neumannFlux *=
                fvGeometry_().boundaryFace[boundaryFaceIdx].area
                * curVolVars_(scvIdx).extrusionFactor();
            Valgrind::CheckDefined(neumannFlux);

            // set the neumann conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isNeumann(eqIdx))
                    continue;
                residual_[scvIdx][eqIdx] += neumannFlux[eqIdx];
            }
        }
    }

    /*!
    * \brief Add outflow boundary conditions for a single sub-control
    *        volume face to the local residual.
    *
    * \param isIt   The intersection iterator of current element
    * \param scvIdx The index of the considered face of the sub-control volume
    * \param boundaryFaceIdx The index of the considered boundary face of the sub control volume
    */
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                            const int scvIdx,
                            const int boundaryFaceIdx)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
        // deal with outflow boundaries
        if (bcTypes.hasOutflow())
        {
            //calculate outflow fluxes
            PrimaryVariables values(0.0);
            asImp_().computeFlux(values, boundaryFaceIdx, true);
            Valgrind::CheckDefined(values);
            
            for (int equationIdx = 0; equationIdx < numEq; ++equationIdx)
            {
                if (!bcTypes.isOutflow(equationIdx) )
                    continue;
                // deduce outflow
                this->residual_[scvIdx][equationIdx] += values[equationIdx];
            }
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
        for (int scvfIdx = 0; scvfIdx < fvGeometry_().numEdges; scvfIdx++)
        {
            int i = fvGeometry_().subContVolFace[scvfIdx].i;
            int j = fvGeometry_().subContVolFace[scvfIdx].j;

            PrimaryVariables flux;

            Valgrind::SetUndefined(flux);
            asImp_().computeFlux(flux, scvfIdx);
            Valgrind::CheckDefined(flux);

            Scalar extrusionFactor =
                (curVolVars_(i).extrusionFactor()
                 + curVolVars_(j).extrusionFactor())
                / 2;
            flux *= extrusionFactor;

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
        storageTerm_.resize(fvGeometry_().numVertices);
        storageTerm_ = 0;

        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int scvIdx = 0; scvIdx < fvGeometry_().numVertices; scvIdx++) {
            Valgrind::SetUndefined(storageTerm_[scvIdx]);
            asImp_().computeStorage(storageTerm_[scvIdx], scvIdx, /*isOldSol=*/false);
            storageTerm_[scvIdx] *=
                fvGeometry_().subContVol[scvIdx].volume
                * curVolVars_(scvIdx).extrusionFactor();
            Valgrind::CheckDefined(storageTerm_[scvIdx]);
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
        for (int scvIdx = 0; scvIdx < fvGeometry_().numVertices; scvIdx++)
        {
            Scalar extrusionFactor =
                curVolVars_(scvIdx).extrusionFactor();

            PrimaryVariables values(0.0);

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            Valgrind::SetUndefined(storageTerm_[scvIdx]);
            Valgrind::SetUndefined(values);
            asImp_().computeStorage(storageTerm_[scvIdx], scvIdx, false);
            asImp_().computeStorage(values, scvIdx, true);
            Valgrind::CheckDefined(storageTerm_[scvIdx]);
            Valgrind::CheckDefined(values);

            storageTerm_[scvIdx] -= values;
            storageTerm_[scvIdx] *=
                fvGeometry_().subContVol[scvIdx].volume
                / problem_().timeManager().timeStepSize()
                * extrusionFactor;
            residual_[scvIdx] += storageTerm_[scvIdx];

            // subtract the source term from the local rate
            Valgrind::SetUndefined(values);
            asImp_().computeSource(values, scvIdx);
            Valgrind::CheckDefined(values);
            values *= fvGeometry_().subContVol[scvIdx].volume * extrusionFactor;
            residual_[scvIdx] -= values;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(residual_[scvIdx]);
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
    const Element &element_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }
    
    /*!
     * \brief Returns a reference to the current element's finite
     *        volume geometry.
     */
    const FVElementGeometry &fvGeometry_() const
    {
        Valgrind::CheckDefined(fvElemGeomPtr_);
        return *fvElemGeomPtr_;
    }

    /*!
     * \brief Returns a reference to the primary variables of
     *           the last time step of the i'th
     *        sub-control volume of the current element.
     */
    const PrimaryVariables &prevPriVars_(const int i) const
    {
        return prevVolVars_(i).priVars();
    }

    /*!
     * \brief Returns a reference to the primary variables of the i'th
     *        sub-control volume of the current element.
     */
    const PrimaryVariables &curPriVars_(const int i) const
    {
        return curVolVars_(i).priVars();
    }

    /*!
     * \brief Returns the j'th primary of the i'th sub-control volume
     *        of the current element.
     */
    Scalar curPriVar_(const int i, const int j) const
    {
        return curVolVars_(i).priVar(j);
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
    const VolumeVariables &curVolVars_(const int i) const
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
    const VolumeVariables &prevVolVars_(const int i) const
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
    const BoundaryTypes &bcTypes_(const int i) const
    {
        return bcTypes_()[i];
    }

protected:
    ElementSolutionVector storageTerm_;
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
