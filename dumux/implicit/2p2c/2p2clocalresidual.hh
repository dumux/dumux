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
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component fully implicit model.
 */

#ifndef DUMUX_2P2C_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_2P2C_LOCAL_RESIDUAL_BASE_HH

#include "2p2cproperties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component fully implicit model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the
 * two-phase two-component flow.
 */
template<class TypeTag>
class TwoPTwoCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
 protected:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    static constexpr unsigned int replaceCompEqIdx =
        GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);

    //! Property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

 public:
    /*!
     * \brief Constructor
     *
     * Sets the mass upwind weight.
     */
    TwoPTwoCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    /*!
     * \brief Evaluate the storage term of the current solution in a
     *        single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void evalPhaseStorage(const Element &element, const int phaseIdx)
    {
        FVElementGeometry fvGeometry;
        fvGeometry.update(this->gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(this->problem_(), element, fvGeometry);
        ElementVolumeVariables elemVolVars;
        elemVolVars.update(this->problem_(), element, fvGeometry, false);

        this->storageTerm_.resize(fvGeometry.numScv);
        this->storageTerm_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeometry;
        this->bcTypesPtr_ = &bcTypes;
        this->prevVolVarsPtr_ = 0;
        this->curVolVarsPtr_ = &elemVolVars;
        evalPhaseStorage_(phaseIdx);
    }

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The sub-control-volume index
     *  \param usePrevSol Based on usePrevSol solution of current or previous time step is used
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        storage = 0;
        if(useMoles) // mole-fraction formulation
        {
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (unsigned int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
                {
                    unsigned int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    storage[eqIdx] += volVars.molarDensity(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.moleFraction(phaseIdx, compIdx);
                 }
                 // this is only processed if one component mass balance equation
                 // is replaced by the total mass balance equation
                 if (replaceCompEqIdx < numComponents)
                     storage[replaceCompEqIdx] +=
                         volVars.molarDensity(phaseIdx)
                         * volVars.saturation(phaseIdx);
            }
            storage *= volVars.porosity();
        }
        else // mass-fraction formulation
        {
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (unsigned int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
                {
                    unsigned int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    storage[eqIdx] += volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.massFraction(phaseIdx, compIdx);
                }
                // this is only processed if one component mass balance equation
                // is replaced by the total mass balance equation
                if (replaceCompEqIdx < numComponents)
                    storage[replaceCompEqIdx] +=
                        volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx);
            }
            storage *= volVars.porosity();
        }
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fIdx The index of the sub-control-volume face
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, bool onBoundary=false) const
    {
        FluxVariables fluxVars(this->problem_(),
                               this->element_(),
                               this->fvGeometry_(),
                               fIdx,
                               this->curVolVars_(),
                               onBoundary);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////

        if(useMoles) // mole-fraction formulation
        {
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                // data attached to upstream and the downstream vertices
                // of the current phase
                const VolumeVariables &up =
                    this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
                const VolumeVariables &dn =
                    this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

                for (unsigned int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
                {
                    unsigned int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    // add advective flux of current component in current
                    // phase
                    if (massUpwindWeight_ > 0.0)
                        // upstream vertex
                        flux[eqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * massUpwindWeight_
                            * up.molarDensity(phaseIdx)
                            * up.moleFraction(phaseIdx, compIdx);
                    if (massUpwindWeight_ < 1.0)
                        // downstream vertex
                        flux[eqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * (1 - massUpwindWeight_)
                            * dn.molarDensity(phaseIdx)
                            * dn.moleFraction(phaseIdx, compIdx);

                    Valgrind::CheckDefined(fluxVars.volumeFlux(phaseIdx));
                    Valgrind::CheckDefined(up.molarDensity(phaseIdx));
                    Valgrind::CheckDefined(up.moleFraction(phaseIdx, compIdx));
                    Valgrind::CheckDefined(dn.molarDensity(phaseIdx));
                    Valgrind::CheckDefined(dn.moleFraction(phaseIdx, compIdx));
                }
                // flux of the total mass balance;
                // this is only processed if one component mass balance equation
                // is replaced by a total mass balance equation
                if (replaceCompEqIdx < numComponents)
                {
                    // upstream vertex
                    if (massUpwindWeight_ > 0.0)
                        flux[replaceCompEqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * massUpwindWeight_
                            * up.molarDensity(phaseIdx);
                    // downstream vertex
                    if (massUpwindWeight_ < 1.0)
                        flux[replaceCompEqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * (1 - massUpwindWeight_)
                            * dn.molarDensity(phaseIdx);
                    Valgrind::CheckDefined(fluxVars.volumeFlux(phaseIdx));
                    Valgrind::CheckDefined(up.molarDensity(phaseIdx));
                    Valgrind::CheckDefined(dn.molarDensity(phaseIdx));

                }

            }
        }
        else // mass-fraction formulation
        {
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                // data attached to upstream and downstream vertices
                // of the current phase
                const VolumeVariables &up =
                    this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
                const VolumeVariables &dn =
                    this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

                for (unsigned int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
                {
                    unsigned int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    // add advective flux of current component in current
                    // phase
                    if (massUpwindWeight_ > 0.0)
                        // upstream vertex
                        flux[eqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * massUpwindWeight_
                            * up.density(phaseIdx)
                            * up.massFraction(phaseIdx, compIdx);
                    if (massUpwindWeight_ < 1.0)
                        // downstream vertex
                        flux[eqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * (1 - massUpwindWeight_)
                            * dn.density(phaseIdx)
                            * dn.massFraction(phaseIdx, compIdx);

                    Valgrind::CheckDefined(fluxVars.volumeFlux(phaseIdx));
                    Valgrind::CheckDefined(up.density(phaseIdx));
                    Valgrind::CheckDefined(up.massFraction(phaseIdx, compIdx));
                    Valgrind::CheckDefined(dn.density(phaseIdx));
                    Valgrind::CheckDefined(dn.massFraction(phaseIdx, compIdx));
                }
                // flux of the total mass balance;
                // this is only processed if one component mass balance equation
                // is replaced by a total mass balance equation
                if (replaceCompEqIdx < numComponents)
                {
                    // upstream vertex
                    if (massUpwindWeight_ > 0.0)
                        flux[replaceCompEqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * massUpwindWeight_
                            * up.density(phaseIdx);
                    // downstream vertex
                    if (massUpwindWeight_ < 1.0)
                        flux[replaceCompEqIdx] +=
                            fluxVars.volumeFlux(phaseIdx)
                            * (1 - massUpwindWeight_)
                            * dn.density(phaseIdx);
                    Valgrind::CheckDefined(fluxVars.volumeFlux(phaseIdx));
                    Valgrind::CheckDefined(up.density(phaseIdx));
                    Valgrind::CheckDefined(dn.density(phaseIdx));

                }

            }
        }
    }

    /*!
     * \brief Evaluates the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const

    {
        if(useMoles) // mole-fraction formulation
        {
            // add diffusive flux of gas component in liquid phase
            Scalar tmp = - (fluxVars.moleFractionGrad(wPhaseIdx)*fluxVars.face().normal);
            tmp *=
                fluxVars.porousDiffCoeff(wPhaseIdx) *
                fluxVars.molarDensity(wPhaseIdx);
            // add the diffusive fluxes only to the component mass balance
            if (replaceCompEqIdx != contiNEqIdx)
                flux[contiNEqIdx] += tmp;
            if (replaceCompEqIdx != contiWEqIdx)
                flux[contiWEqIdx] -= tmp;

            // add diffusive flux of liquid component in non-wetting phase
            tmp = -(fluxVars.moleFractionGrad(nPhaseIdx)*fluxVars.face().normal);
            tmp *=
                fluxVars.porousDiffCoeff(nPhaseIdx) *
                fluxVars.molarDensity(nPhaseIdx);
            // add the diffusive fluxes only to the component mass balance
            if (replaceCompEqIdx != contiWEqIdx)
                flux[contiWEqIdx] += tmp;
            if (replaceCompEqIdx != contiNEqIdx)
                flux[contiNEqIdx] -= tmp;
        }
        else // mass-fraction formulation
        {
            // add diffusive flux of gas component in liquid phase
            Scalar tmp = - (fluxVars.moleFractionGrad(wPhaseIdx)*fluxVars.face().normal);
            tmp *=
                fluxVars.porousDiffCoeff(wPhaseIdx) *
                fluxVars.molarDensity(wPhaseIdx);
            // add the diffusive fluxes only to the component mass balance
            if (replaceCompEqIdx != contiNEqIdx)
                flux[contiNEqIdx] += tmp * FluidSystem::molarMass(nCompIdx);
            if (replaceCompEqIdx != contiWEqIdx)
                flux[contiWEqIdx] -= tmp * FluidSystem::molarMass(wCompIdx);

            // add diffusive flux of liquid component in non-wetting phase
            tmp = -(fluxVars.moleFractionGrad(nPhaseIdx)*fluxVars.face().normal);
            tmp *=
                fluxVars.porousDiffCoeff(nPhaseIdx) *
                fluxVars.molarDensity(nPhaseIdx);
            // add the diffusive fluxes only to the component mass balance
            if (replaceCompEqIdx != contiWEqIdx)
                flux[contiWEqIdx] += tmp * FluidSystem::molarMass(wCompIdx);
            if (replaceCompEqIdx != contiNEqIdx)
                flux[contiNEqIdx] -= tmp * FluidSystem::molarMass(nCompIdx);
        }
    }

    /*!
     * \brief Evaluates the source term
     *
     * \param source The source/sink in the sub-control volume
     * \param scvIdx The index of the sub-control volume
     *
     * Be careful what you use, mole or mass fraction! Think of the units!
     * If a total mass balance is used, the sum of both components has to be specified as source.
     */
    void computeSource(PrimaryVariables& source, const int scvIdx)
    {
        this->problem_().solDependentSource(source,
                                     this->element_(),
                                     this->fvGeometry_(),
                                     scvIdx,
                                     this->curVolVars_());
    }

 protected:
    void evalPhaseStorage_(const int phaseIdx)
    {
        if(useMoles) // mole-fraction formulation
        {
            // evaluate the storage terms of a single phase
            for (int i=0; i < this->fvGeometry_().numScv; i++) {
                PrimaryVariables &storage = this->storageTerm_[i];
                const ElementVolumeVariables &elemVolVars = this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[i];

                // compute storage term of all components within all phases
                storage = 0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    storage[eqIdx] += volVars.molarDensity(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.moleFraction(phaseIdx, compIdx);
                }

                storage *= volVars.porosity();
                storage *= this->fvGeometry_().subContVol[i].volume;
            }
        }
        else // mass-fraction formulation
        {
            // evaluate the storage terms of a single phase
            for (int i=0; i < this->fvGeometry_().numScv; i++) {
                PrimaryVariables &storage = this->storageTerm_[i];
                const ElementVolumeVariables &elemVolVars = this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[i];

                // compute storage term of all components within all phases
                storage = 0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    int eqIdx = (compIdx == wCompIdx) ? contiWEqIdx : contiNEqIdx;
                    storage[eqIdx] += volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.massFraction(phaseIdx, compIdx);
                }

                storage *= volVars.porosity();
                storage *= this->fvGeometry_().subContVol[i].volume;
            }
        }
    }

    /*!
     * \brief Returns the equation index of the first mass-balance equation
     *        of the component (used for loops)
     *
     * Returns the equation index of the first mass-balance equation
     * of the component (used for loops) if one component mass balance
     * is replaced by the total mass balance, this is the index
     * of the remaining component mass-balance equation.
     */
    unsigned int contiCompIdx1_() const {
        switch (replaceCompEqIdx)
        {
        case contiWEqIdx: return contiNEqIdx;
        case contiNEqIdx: return contiWEqIdx;
        default:          return 0;
        }
    }

    /*!
     * \brief Returns the equation index of the second mass balance
     *        of the component (used for loops)
     *
     * Returns the equation index of the second mass balance
     * of the component (used for loops)
     * if one component mass balance is replaced by the total mass balance
     * (replaceCompEqIdx < 2), this index is the same as contiCompIdx1().
     */
    unsigned int contiCompIdx2_() const {
        switch (replaceCompEqIdx)
        {
        case contiWEqIdx: return contiNEqIdx;
        case contiNEqIdx: return contiWEqIdx;
        default:          return numComponents-1;
        }
    }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

 private:
    Scalar massUpwindWeight_;
};

} // end namespace

#endif
