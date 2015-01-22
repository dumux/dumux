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

#include <dumux/implicit/nonisothermal/nilocalresidual.hh>
#include <dumux/implicit/2p2c/2p2cproperties.hh>

#define VELOCITY_OUTPUT 1 // uncomment this line if an output of the velocity is needed

namespace Dumux
{

/*!
 * \ingroup ImplicitLocalResidual
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
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
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum {
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
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
    //! \brief Constructor
    TwoPTwoCNICouplingLocalResidual()
    { };

    /*!
     * \brief Implementation of the boundary evaluation
     */
    void evalBoundary_()
    {
        ParentType::evalBoundary_();

        typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
        typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());

        // evaluate Dirichlet-like coupling conditions
        for (int idx = 0; idx < this->fvGeometry_().numScv; idx++)
        {
            // evaluate boundary conditions for the intersections of
            // the current element
            IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
            const IntersectionIterator &endIt = this->gridView_().iend(this->element_());
            for (; isIt != endIt; ++isIt)
            {
                // handle only intersections on the boundary
                if (!isIt->boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int fIdx = isIt->indexInInside();
                const int numFaceVertices = refElement.size(fIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertexIdx = 0; faceVertexIdx < numFaceVertices; ++faceVertexIdx)
                {
                    const int vIdx = refElement.subEntity(fIdx, 1, faceVertexIdx, dim);
                    const int boundaryFaceIdx = this->fvGeometry_().boundaryFaceIndex(fIdx, faceVertexIdx);
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (vIdx != idx)
                        continue;

                    //for the corner points, the boundary flux across the vertical non-coupling boundary faces
                    //has to be calculated to fulfill the mass balance
                    //convert suddomain intersection into multidomain intersection and check whether it is an outer boundary
                    if(!GridView::Grid::multiDomainIntersection(*isIt).neighbor()
                             && this->boundaryHasNeumann_(this->bcTypes_(idx)))
                     {
                         const DimVector& globalPos = this->fvGeometry_().subContVol[idx].global;
                         //problem specific function, in problem orientation of interface is known
                         if(this->problem_().isInterfaceCornerPoint(globalPos))
                         {
                             //check whether the second boundary node on the vertical edge is Neumann-Null
                             const int numVertices = refElement.size(dim);
                             bool evalBoundaryFlux = false;
                             for (int equationIdx = 0; equationIdx < numEq; ++equationIdx)
                             {
                                 for(int i= 0; i < numVertices; i++)
                                 {
                                     //if vertex is on boundary and not the coupling vertex: check whether a Neumann condition is set
                                     //in case of Neumann Null the boundary flux must not be calculated
                                     if(this->model_().onBoundary(this->element_(), i) && i!=vIdx)
                                         if (!this->bcTypes_(i).isNeumann(equationIdx) && !this->bcTypes_(i).isOutflow(equationIdx))
                                             evalBoundaryFlux = true;
                                 }

                                 PrimaryVariables values(0.0);
                                 //calculate the actual boundary fluxes and add to residual
                                 this->asImp_()->computeFlux(values, boundaryFaceIdx, true /*on boundary*/);
                                 if(evalBoundaryFlux)
                                     this->residual_[idx][equationIdx] += values[equationIdx];
                             }

                         }
                     }

                    if (boundaryHasCoupling_(this->bcTypes_(idx)))
                        evalCouplingVertex_(idx);
                }
            }
        }
    }

    /*!
     * \brief Evaluates the phase storage
     */
    Scalar evalPhaseStorage(const int scvIdx) const
    {
        Scalar phaseStorage = computePhaseStorage(scvIdx, false);
        Scalar oldPhaseStorage = computePhaseStorage(scvIdx, true);
        Valgrind::CheckDefined(phaseStorage);
        Valgrind::CheckDefined(oldPhaseStorage);

        phaseStorage -= oldPhaseStorage;
        phaseStorage *= this->fvGeometry_().subContVol[scvIdx].volume
                / this->problem_().timeManager().timeStepSize()
                * this->curVolVars_(scvIdx).extrusionFactor();
        Valgrind::CheckDefined(phaseStorage);

        return phaseStorage;
    }

    /*!
     * \brief Compute storage term of all components within all phases
     */
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
     * \brief Compute the fluxes within the different fluid phases. This is
     *        merely for the computation of flux output.
     */
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
    void evalCouplingVertex_(const int scvIdx)
    {
        const VolumeVariables &volVars = this->curVolVars_()[scvIdx];

        // set pressure as part of the momentum coupling
        if (this->bcTypes_(scvIdx).isCouplingOutflow(massBalanceIdx))
            this->residual_[scvIdx][massBalanceIdx] = volVars.pressure(nPhaseIdx);

        // set mass fraction; TODO: this is fixed to contiWEqIdx so far!
        if (this->bcTypes_(scvIdx).isCouplingOutflow(contiWEqIdx))
            this->residual_[scvIdx][contiWEqIdx] = volVars.massFraction(nPhaseIdx, wCompIdx);

        // set temperature as part of the coupling
        if (this->bcTypes_(scvIdx).isCouplingOutflow(energyEqIdx))
            this->residual_[scvIdx][energyEqIdx] = volVars.temperature();
    }

    /*!
     * \brief Check if one of the boundary conditions is coupling.
     */
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
    Scalar elementFluxes(const int scvIdx)
    {
        return elementFluxes_[scvIdx];
    }

private:
    ElementFluxVector elementFluxes_;
};

} // namespace Dumux

#endif // DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH
