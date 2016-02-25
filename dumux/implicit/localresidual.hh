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
    friend class typename GET_PROP_TYPE(TypeTag, LocalJacobian);
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
        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());

        asImp_().eval(element, fvGeometry_(), bcTypes);
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

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());
        bcTypesPtr_ = &bcTypes;

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
    void evalFluxes(const Element &element)
    {
        elemPtr_ = &element;

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem_(), element, fvGeometry_());

        residual_.resize(fvGeometry_().numScv);
        residual_ = 0;

        bcTypesPtr_ = &bcTypes;
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
              const ElementBoundaryTypes &bcTypes)
    {
        elemPtr_ = &element;
        bcTypesPtr_ = &bcTypes;

        // resize the vectors for all terms
        int numScv = fvGeometry_().numScv();
        residual_.resize(numScv);
        storageTerm_.resize(numScv);

        residual_ = 0.0;
        storageTerm_ = 0.0;

        evalFluxes_();

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
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink in the sub-control volume for each phase
     * \param scvIdx The index of the sub-control volume
     *
     */
    PrimaryVariables computeSource(const SubControlVolume &scv)
    {
        PrimaryVariables source(0);

        source += this->problem_().source(element_(), scv);

        // add contribution from possible point sources
        source += this->problem_().scvPointSources(element_(), scv);
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

    void evalFluxes_()
    {
        // calculate the mass flux over the scv faces and subtract
        for (auto&& scvFace : fvGeometry_().scvfs())
        {
            PrimaryVariables flux = asImp_().computeFlux_(scvFace);

            if (!isBox)
                residual_[0] += flux;
            else
            {
                const auto& insideScv = problem_().fvGeometries().subControlVolume(scvFace.insideScvIdx());
                const auto& outsideScv = problem_().fvGeometries().subControlVolume(scvFace.outsideScvIdx());

                residual[insideScv.indexInElement()] += flux;
                residual[outsideScv.indexInElement()] -= flux;
            }
        }
    }

    PrimaryVariables evalFlux_(const int scvFaceIdx)
    {
        auto&& scvFace = problem_().model().fvGeometries().subControlVolumeFace(scvFaceIdx);

        return evalFlux_(scvFace);
    }

    PrimaryVariables evalFlux_(const subControlVolumeFace &scvFace)
    {
        return asImp_().computeFlux_(scvFace);
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
        for (auto&& scv : fvGeometry_().scvs())
        {
            int scvIdx = scv.indexInElement();

            storageTerm_[scvIdx] = asImp_().computeStorage(scv, /*isOldSol=*/false);
            storageTerm_[scvIdx] *= scv.volume();
            storageTerm_[scvIdx] *= problem_().model().curVolVars(scv).extrusionFactor();
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
        for (auto&& scv : fvGeometry_().scvs())
        {
            int scvIdx = scv.indexInElement();
            Scalar prevExtrusionFactor = problem_().model().prevVolVars(scv).extrusionFactor();
            Scalar curExtrusionFactor = problem_().model().curVolVars(scv).extrusionFactor();

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // We might need a more explicit way for
            // doing the time discretization...
            PrimaryVariables prevStorage = asImp_().computeStorage(scv, true);
            PrimaryVariables curStorage = asImp_().computeStorage(scv, false);

            storageTerm_[scvIdx] = curStorage*curExtrusionFactor;
            storageTerm_[scvIdx] -= prevStorage*prevExtrusionFactor;
            storageTerm_[scvIdx] *= scv.volume();
            storageTerm_[scvIdx] /= problem_().timeManager().timeStepSize();

            // add the storage term to the residual
            residual_[scvIdx] += storageTerm_[scvIdx];

            // subtract the source term from the local rate
            PrimaryVariables source = asImp_().computeSource(scv);
            source *= scv.volume()*curExtrusionFactor;

            residual_[scvIdx] -= source;
        }
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
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
        return model_().fvGeometries(element_());
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

    const ElementBoundaryTypes *bcTypesPtr_;
};

}

#endif
