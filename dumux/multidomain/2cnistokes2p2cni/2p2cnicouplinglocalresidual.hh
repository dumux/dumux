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
 * \brief Extending the TwoPTwoCNILocalResidual by the required functions for
 *        a coupled application.
 */
#ifndef DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH
#define DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH

#include <dune/common/deprecated.hh>

#include <dumux/porousmediumflow/nonisothermal/implicit/localresidual.hh>
#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitLocalResidual
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief Extending the TwoPTwoCNILocalResidual by the required functions for
 *        a coupled application.
 */
template<class TypeTag>
class TwoPTwoCNICouplingLocalResidual : public NILocalResidual<TypeTag>
{
    typedef NILocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { useMoles = GET_PROP_VALUE(TypeTag, UseMoles) };
    enum {
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum {
        massBalanceIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),
        contiWEqIdx = Indices::contiWEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum {
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };
    enum { phaseIdx = nPhaseIdx }; // index of the phase for the phase flux calculation
    enum { compIdx = wCompIdx}; // index of the component for the phase flux calculation

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1> > ElementFluxVector;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \brief Implementation of the boundary evaluation for the Darcy model
     *
     * Evaluate one part of the Dirichlet-like coupling conditions for a single
     * sub-control volume face; rest is done in the local coupling operator
     */
    void evalBoundary_()
    {
        ParentType::evalBoundary_();

        typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
        typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());

        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvGeometry_().subContVol[scvIdx].inner)
                continue;

            // evaluate boundary conditions for the intersections of the current element
            for (const auto& intersection : Dune::intersections(this->gridView_(), this->element_()))
            {
                // handle only intersections on the boundary
                if (!intersection.boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int fIdx = intersection.indexInInside();
                const int numFaceVertices = refElement.size(fIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertexIdx = 0; faceVertexIdx < numFaceVertices; ++faceVertexIdx)
                {
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (refElement.subEntity(fIdx, 1, faceVertexIdx, dim) != scvIdx)
                        continue;

                    const VolumeVariables &volVars = this->curVolVars_()[scvIdx];

                    // set pressure as part of the momentum coupling
                    static_assert(GET_PROP_VALUE(TypeTag, Formulation) == TwoPTwoCFormulation::pnsw,
                                  "This coupling condition is only implemented for a pnsw formulation.");
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(massBalanceIdx))
                        this->residual_[scvIdx][massBalanceIdx] = volVars.pressure(nPhaseIdx);

                    // set mass/mole fraction for transported component
                    static_assert(GET_PROP_VALUE(TypeTag, NumComponents) == 2,
                                  "This coupling condition is only implemented for two components.");
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(contiWEqIdx))
                    {
                        if (useMoles)
                            this->residual_[scvIdx][contiWEqIdx] = volVars.moleFraction(nPhaseIdx, wCompIdx);
                        else
                            this->residual_[scvIdx][contiWEqIdx] = volVars.massFraction(nPhaseIdx, wCompIdx);
                    }

                    // set temperature
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(energyEqIdx))
                        this->residual_[scvIdx][energyEqIdx] = volVars.temperature();
                }
            }
        }
    }

    /*!
     * \brief Evaluates the time derivative of the storage term
     *
     * \param storage The vector in which the result is written
     * \param scvIdx The sub-control-volume index
     */
    void evalStorageDerivative(PrimaryVariables &storage, const int scvIdx) const
    {
        PrimaryVariables result;
        this->computeStorage(result, scvIdx, false);
        Valgrind::CheckDefined(result);
        storage = result;
        this->computeStorage(result, scvIdx, true);
        Valgrind::CheckDefined(result);
        storage -= result;

        storage *= this->fvGeometry_().subContVol[scvIdx].volume
                   / this->problem_().timeManager().timeStepSize()
                   * this->curVolVars_(scvIdx).extrusionFactor();
    }

    /*!
     * \brief Evaluates the time derivative of the phase storage
     */
    DUNE_DEPRECATED_MSG("evalPhaseStorageDerivative() is deprecated. Use evalStorageDerivative() instead.")
    Scalar evalPhaseStorageDerivative(const int scvIdx) const
    {
        Scalar result = computePhaseStorage(scvIdx, false);
        Scalar oldPhaseStorage = computePhaseStorage(scvIdx, true);
        Valgrind::CheckDefined(result);
        Valgrind::CheckDefined(oldPhaseStorage);

        result -= oldPhaseStorage;
        result *= this->fvGeometry_().subContVol[scvIdx].volume
                  / this->problem_().timeManager().timeStepSize()
                  * this->curVolVars_(scvIdx).extrusionFactor();
        Valgrind::CheckDefined(result);

        return result;
    }

    /*!
     * \brief Compute storage term of all components within all phases
     */
    DUNE_DEPRECATED_MSG("computePhaseStorage() is deprecated.")
    Scalar computePhaseStorage(const int scvIdx, bool usePrevSol) const
    {
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        Scalar phaseStorage = volVars.density(phaseIdx)
                * volVars.saturation(phaseIdx)
                * volVars.massFraction(phaseIdx, compIdx);
        phaseStorage *= volVars.porosity();

        return phaseStorage;
    }

    /*!
     * \brief Compute the fluxes merely for the flux output.
     *
     * \param scvIdx The sub-control-volume index
     */
    DUNE_DEPRECATED_MSG("evalComponentFluxes() has been deprecated.")
    PrimaryVariables evalComponentFluxes(const int scvIdx)
    {
        const unsigned int numScv = this->fvGeometry_().numScv;
        PrimaryVariables fluxes[numScv];
        for (int fIdx = 0; fIdx < numScv; fIdx++)
            fluxes[fIdx] = 0.0;

        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int fIdx = 0; fIdx < numScv; fIdx++)
        {
            FluxVariables fluxVars(this->problem_(),
                                   this->element_(),
                                   this->fvGeometry_(),
                                   fIdx,
                                   this->curVolVars_());

            int i = this->fvGeometry_().subContVolFace[fIdx].i;
            int j = this->fvGeometry_().subContVolFace[fIdx].j;

            PrimaryVariables flux(0.);
            PrimaryVariables fluxTemp(0.);
            this->computeAdvectiveFlux(fluxTemp, fluxVars);
            flux += fluxTemp;
            this->computeAdvectiveFlux(fluxTemp, fluxVars);
            flux += fluxTemp;

            const Scalar extrusionFactor = 0.5 * this->curVolVars_(i).extrusionFactor()
                                           + 0.5 * this->curVolVars_(j).extrusionFactor();
            flux *= extrusionFactor;

            fluxes[i] += flux;
            fluxes[j] -= flux;
        }
        return fluxes[scvIdx];
    }

    /*!
     * \brief Compute the fluxes within the different fluid phases. This is
     *        merely for the computation of flux output.
     */
    DUNE_DEPRECATED_MSG("evalPhaseFluxes() has been deprecated.")
    void evalPhaseFluxes()
    {
        elementFluxes_.resize(this->fvGeometry_().numScv);
        elementFluxes_ = 0.;

        Scalar flux(0.);

        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int fIdx = 0; fIdx < this->fvGeometry_().numScvf; fIdx++)
        {
            FluxVariables vars(this->problem_(),
                               this->element_(),
                               this->fvGeometry_(),
                               fIdx,
                               this->curVolVars_());

            int i = this->fvGeometry_().subContVolFace[fIdx].i;
            int j = this->fvGeometry_().subContVolFace[fIdx].j;

            const Scalar extrusionFactor =
                (this->curVolVars_(i).extrusionFactor()
                 + this->curVolVars_(j).extrusionFactor())
                / 2;
            flux = computeAdvectivePhaseFlux(vars) +
                    computeDiffusivePhaseFlux(vars);
            flux *= extrusionFactor;

            elementFluxes_[i] += flux;
            elementFluxes_[j] -= flux;
        }
    }

    /*!
     * \brief Returns the advective fluxes within the different phases.
     */
    DUNE_DEPRECATED_MSG("computeAdvectivePhaseFlux() is deprecated.")
    Scalar computeAdvectivePhaseFlux(const FluxVariables &fluxVars) const
    {
        Scalar advectivePhaseFlux = 0.;
        const Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up =
            this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
        const VolumeVariables &dn =
            this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

        if (massUpwindWeight > 0.0)
            // upstream vertex
            advectivePhaseFlux +=
                fluxVars.volumeFlux(phaseIdx)
                * massUpwindWeight
                * up.density(phaseIdx)
                * up.massFraction(phaseIdx, compIdx);
        if (massUpwindWeight < 1.0)
            // downstream vertex
            advectivePhaseFlux +=
                fluxVars.volumeFlux(phaseIdx)
                * (1 - massUpwindWeight)
                * dn.density(phaseIdx)
                * dn.massFraction(phaseIdx, compIdx);

        return advectivePhaseFlux;
    }

    /*!
     * \brief Returns the diffusive fluxes within the different phases.
     */
    DUNE_DEPRECATED_MSG("computeAdvectivePhaseFlux() is deprecated.")
    Scalar computeDiffusivePhaseFlux(const FluxVariables &fluxVars) const
    {
        // add diffusive flux of gas component in liquid phase
        Scalar diffusivePhaseFlux = fluxVars.moleFractionGrad(phaseIdx)*fluxVars.face().normal;
        diffusivePhaseFlux *= -fluxVars.porousDiffCoeff(phaseIdx) *
               fluxVars.molarDensity(phaseIdx);
        diffusivePhaseFlux *= FluidSystem::molarMass(compIdx);

        return diffusivePhaseFlux;
    }

    /*!
     * \brief Set the Dirichlet-like conditions for the coupling
     *        and replace the existing residual
     *
     * \param scvIdx Sub control vertex index for the coupling condition
     */
    DUNE_DEPRECATED_MSG("evalCouplingVertex_ is deprecated. Its functionality is now included in evalBoundary_.")
    void evalCouplingVertex_(const int scvIdx)
    {
        const VolumeVariables &volVars = this->curVolVars_()[scvIdx];

        // set pressure as part of the momentum coupling
        if (this->bcTypes_(scvIdx).isCouplingDirichlet(massBalanceIdx))
            this->residual_[scvIdx][massBalanceIdx] = volVars.pressure(nPhaseIdx);

        // set mass fraction;
        if (this->bcTypes_(scvIdx).isCouplingDirichlet(contiWEqIdx))
            this->residual_[scvIdx][contiWEqIdx] = volVars.massFraction(nPhaseIdx, wCompIdx);

        // set temperature as part of the coupling
        if (this->bcTypes_(scvIdx).isCouplingDirichlet(energyEqIdx))
            this->residual_[scvIdx][energyEqIdx] = volVars.temperature();
    }

    /*!
     * \brief Check if one of the boundary conditions is coupling.
     */
    DUNE_DEPRECATED_MSG("boundaryHasCoupling_ is deprecated")
    bool boundaryHasCoupling_(const BoundaryTypes& bcTypes) const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isCouplingInflow(eqIdx) || bcTypes.isCouplingOutflow(eqIdx))
                return true;
        return false;
    }

    /*!
     * \brief Check if one of the boundary conditions is Neumann.
     */
    DUNE_DEPRECATED_MSG("boundaryHasNeumann_ is unused in dumux and therefore deprecated")
    bool boundaryHasNeumann_(const BoundaryTypes& bcTypes) const
    {
        bool hasNeumann = false;
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isNeumann(eqIdx) || bcTypes.isNeumann(eqIdx))
                hasNeumann = true;
        return hasNeumann;
    }

    /*!
     * \brief Returns the fluxes of the specified sub control volume
     */
    DUNE_DEPRECATED_MSG("elementFluxes() is deprecated")
    Scalar elementFluxes(const int scvIdx)
    {
        return elementFluxes_[scvIdx];
    }

private:
    ElementFluxVector elementFluxes_;
};

} // namespace Dumux

#endif // DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH
