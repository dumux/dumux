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
 * \brief Calculates the element-wise residual of fully-implicit models.
 */
#ifndef DUMUX_IMPLICIT_LOCAL_RESIDUAL_HH
#define DUMUX_IMPLICIT_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual matrix for models
 *        using a fully implicit discretization.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class ImplicitLocalResidual
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    // copying the local residual class is not a good idea
    ImplicitLocalResidual(const ImplicitLocalResidual &);

public:
    ImplicitLocalResidual()
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

        residual_.resize(fvGeometry_().numScv);
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
        for (unsigned int i = 0; i < prevVolVars.size(); i++) {
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
        int numScv = fvGeometry_().numScv;
        residual_.resize(numScv);
        storageTerm_.resize(numScv);

        residual_ = 0.0;
        storageTerm_ = 0.0;

        asImp_().evalFluxes_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numScv; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND

        asImp_().evalVolumeTerms_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numScv; i++) {
            Valgrind::CheckDefined(residual_[i]);
        }
#endif // HAVE_VALGRIND

        // evaluate the boundary conditions
        asImp_().evalBoundary_();

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numScv; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink in the sub-control volume for each phase
     * \param scvIdx The index of the sub-control volume
     *
     */
    void computeSource(PrimaryVariables &source, const int scvIdx)
    {
        this->problem_().solDependentSource(source,
                                            element_(),
                                            fvGeometry_(),
                                            scvIdx,
                                            curVolVars_());

        // add contribution from possible point sources
        this->problem_().scvPointSources(source,
                                         element_(),
                                         fvGeometry_(),
                                         scvIdx,
                                         curVolVars_());
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
        // Dirichlet boundary conditions are treated differently
        // depending on the spatial discretization: for box,
        // they are incorporated in a strong sense, whereas for
        // cell-centered, they are treated by means of fluxes.
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            if (bcTypes_().hasNeumann() || bcTypes_().hasOutflow())
                asImp_().evalBoundaryFluxes_();

            if (bcTypes_().hasDirichlet())
                asImp_().evalDirichlet_();
        }
        else
        {
            asImp_().evalBoundaryFluxes_();

            // additionally treat mixed D/N conditions in a strong sense
            if (bcTypes_().hasDirichlet())
                asImp_().evalDirichlet_();
        }

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < fvGeometry_().numScv; i++)
            Valgrind::CheckDefined(residual_[i]);
#endif // HAVE_VALGRIND
    }

    void evalDirichlet_()
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "evalDirichlet_ has been called but is not implemented");
    }

    /*!
     * \brief Set the local residual to the storage terms of all
     *        sub-control volumes of the current element.
     */
    void evalStorage_()
    {
        storageTerm_.resize(fvGeometry_().numScv);
        storageTerm_ = 0;

        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int scvIdx = 0; scvIdx < fvGeometry_().numScv; scvIdx++) {
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
        for (int scvIdx = 0; scvIdx < fvGeometry_().numScv; scvIdx++)
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
    { return *problemPtr_; }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); }

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
