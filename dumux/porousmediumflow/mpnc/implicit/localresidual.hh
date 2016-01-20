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
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_HH

#include "fluxvariables.hh"
#include "diffusion/diffusion.hh"
#include "energy/localresidual.hh"
#include "mass/localresidual.hh"
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitLocalResidual
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the
 * MpNc flow.
 */
template<class TypeTag>
class MPNCLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    enum {enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum {numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {phase0NcpIdx = Indices::phase0NcpIdx};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, numEnergyEquations> EnergyResid;
    typedef MPNCLocalResidualMass<TypeTag, enableKinetic> MassResid;

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     *  \param scvIdx The SCV (sub-control-volume) index
     *
     */
    void computeStorage(PrimaryVariables &storage,
                        const unsigned int scvIdx,
                        const bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        storage =0;

        // compute mass and energy storage terms
        MassResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
        EnergyResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within all sub-control volumes of an
     *        element.
     *
     *        \param phaseStorage The conserved quantity within the phase in the whole domain
     *        \param element The finite element
     *        \param phaseIdx The index of the fluid phase
     */
    void addPhaseStorage(PrimaryVariables &phaseStorage,
                         const Element &element,
                         const unsigned int phaseIdx) const
    {
        // create a finite volume element geometry
        FVElementGeometry fvGeometry;
        fvGeometry.update(this->gridView_(), element);

        // calculate volume variables
        ElementVolumeVariables elemVolVars;
        this->model_().setHints(element, elemVolVars);
        elemVolVars.update(this->problem_(),
                           element,
                           fvGeometry,
                           /*useOldSolution=*/false);

        // calculate the phase storage for all sub-control volumes
        for (int scvIdx=0;
             scvIdx < fvGeometry.numScv;
             scvIdx++)
        {
            PrimaryVariables tmpPriVars(0.0);

            // compute mass and energy storage terms in terms of
            // averaged quantities
            MassResid::addPhaseStorage(tmpPriVars,
                                       elemVolVars[scvIdx],
                                       phaseIdx);
            EnergyResid::addPhaseStorage(tmpPriVars,
                                         elemVolVars[scvIdx],
                                         phaseIdx);

            // multiply with volume of sub-control volume
            tmpPriVars *= fvGeometry.subContVol[scvIdx].volume;

            // Add the storage of the current SCV to the total storage
            phaseStorage += tmpPriVars;
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     *   \param scvIdx The SCV (sub-control-volume) index
     *   \param source The source/sink in the sub-control volume for each component
     */
    void computeSource(PrimaryVariables &source,
                       const unsigned int scvIdx)
     {
        Valgrind::SetUndefined(source);
        ParentType::computeSource(source, scvIdx);

        const VolumeVariables &volVars = this->curVolVars_(scvIdx);

        PrimaryVariables tmp(0);
        MassResid::computeSource(tmp, volVars);
        source += tmp;
        Valgrind::CheckDefined(source);

        /*
         *      EnergyResid also called in the MassResid
         *      1) Makes some sense because energy is also carried by mass
         *      2) The mass transfer between the phases is needed.
         */
//        tmp = 0.;
//        EnergyResid::computeSource(tmp, volVars);
//        source += tmp;
//        Valgrind::CheckDefined(source);
     };


    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux,
                     const unsigned int fIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars(this->problem_(),
                               this->element_(),
                               this->fvGeometry_(),
                               fIdx,
                               this->curVolVars_(),
                               onBoundary);

        flux = 0.0;
        MassResid::computeFlux(flux, fluxVars, this->curVolVars_() );
        Valgrind::CheckDefined(flux);
/*
 *      EnergyResid also called in the MassResid
 *      1) Makes some sense because energy is also carried by mass
 *      2) The component-wise mass flux in each phase is needed.
 */
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     *        \param element The finite element
     */
    void eval(const Element &element)
    { ParentType::eval(element); }

    /*!
     * \brief Evaluate the local residual.
     *
     *      \param element The finite element
     *      \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *      \param prevElemVolVars The element volume variables of the previous timestep
     *      \param curElemVolVars The element volume variables of the current timestep
     *      \param bcType The types of the boundary conditions for all vertices of the element
     */
    void eval(const Element &element,
              const FVElementGeometry &fvGeometry,
              const ElementVolumeVariables &prevElemVolVars,
              const ElementVolumeVariables &curElemVolVars,
              const ElementBoundaryTypes &bcType)
    {
        ParentType::eval(element,
                         fvGeometry,
                         prevElemVolVars,
                         curElemVolVars,
                         bcType);

        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox)
            || !bcType.hasDirichlet())
        {
            for (int i = 0; i < this->fvGeometry_().numScv; ++i) {
                // add the two auxiliary equations, make sure that the
                // dirichlet boundary condition is conserved
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    if (!bcType[i].isDirichlet(phase0NcpIdx + phaseIdx))
                    {
                        this->residual_[i][phase0NcpIdx + phaseIdx] =
                            this->curVolVars_(i).phaseNcp(phaseIdx);
                    }
                }
            }
        }
    }


    /*!
     * \brief Add Dirichlet boundary conditions for a single intersection
     *
     * Sets the Dirichlet conditions in a strong sense, in contrast to
     * the general handling in CCLocalResidual.
     *
     *      \param isIt
     *      \param bcTypes
     */
    template <class IntersectionIterator>
    void evalDirichletSegment_(const IntersectionIterator &isIt,
                               const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the Dirichlet boundary fluxes
        PrimaryVariables values;
        Valgrind::SetUndefined(values);
        this->problem_().dirichlet(values, *isIt);
        Valgrind::CheckDefined(values);

        // set Dirichlet conditions in a strong sense
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                this->residual_[0][eqIdx]
                  = this->curPriVar_(0, pvIdx) - values[pvIdx];
            }
            else if (eqIdx >= phase0NcpIdx)
            {
                int phaseIdx = eqIdx - phase0NcpIdx;
                this->residual_[0][eqIdx] =
                    this->curVolVars_(0).phaseNcp(phaseIdx);
            }
        }
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace

#endif
