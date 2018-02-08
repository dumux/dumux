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
 * \brief Element-wise calculation of the residual for the two-phase fully implicit model.
 */
#ifndef DUMUX_DECOUPLED_TWOP_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_DECOUPLED_TWOP_LOCAL_RESIDUAL_BASE_HH

#include "dumux/porousmediumflow/2p/implicit/properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase fully implicit model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag>
class DecoupledTwoPLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        dim = GridView::dimension

    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    /*!
     * \brief Constructor
     *
     * Sets the upwind weight.
     */
    DecoupledTwoPLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.

        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];
//
//         std::cout << "volVars.density(wPhaseIdx) = " << volVars.density(wPhaseIdx) << std::endl;
//          std::cout << "volVars.saturation(wPhaseIdx) = " << volVars.saturation(wPhaseIdx) << std::endl;

        storage = Scalar(0);

        // wetting phase mass
        storage[contiWEqIdx] = volVars.density(wPhaseIdx) * volVars.effPorosity()
                * volVars.saturation(wPhaseIdx);

        storage[contiNEqIdx] = volVars.density(nPhaseIdx) * volVars.effPorosity()
                * volVars.saturation(nPhaseIdx);

        int eIdx = this->problem_().model().elementMapper().index(this->element_());

//         if ( (this->problem_().coupled() == true) && (eIdx == 221) )
//         {
//             std::cout.precision(15);
//             std::cout << "usePrevSol = " << usePrevSol << std::endl;
//             std::cout << "volVars.saturation(wPhaseIdx) = " << volVars.saturation(wPhaseIdx) << std::endl;
//             std::cout << "volVars.density(wPhaseIdx) = " << volVars.density(wPhaseIdx) << std::endl;
//             std::cout << "volVars.effPorosity() = " << volVars.effPorosity() << std::endl;
//             std::cout << "storage[contiWEqIdx] at [" << scvIdx << "] = " << storage[contiWEqIdx] << std::endl;
//         }
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fIdx The index of the sub-control-volume face
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    void computeFlux(PrimaryVariables &flux, int fIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->curVolVars_(),
                        onBoundary);

        FluxVariables fluxVarsOld;
        fluxVarsOld.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->prevVolVars_(),
                        onBoundary);

        flux = 0;

        asImp_()->computeAdvectiveFlux(flux, fluxVars, fluxVarsOld);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);

//         if (this->problem_().coupled() == true)
//         {
//             std::cout.precision(15);
//             std::cout << "flux at [" << fIdx << "] = " << flux << std::endl;
//         }
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars, const FluxVariables &fluxVarsOld) const
    {
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

//             std::cout << "up.effPorosity() = " << up.effPorosity() << std::endl;

            // add advective flux of current phase
            int eqIdx = (phaseIdx == wPhaseIdx) ? contiWEqIdx : contiNEqIdx;
            flux[eqIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                *
                ((    massUpwindWeight_)*up.density(phaseIdx)
                 +
                 (1 - massUpwindWeight_)*dn.density(phaseIdx));

//             flux[eqIdx] += up.effPorosity() * up.saturation(phaseIdx)
//                         * up.density(phaseIdx) * fluxVars.timeDerivUNormal();

            int eIdx = this->problem_().model().elementMapper().index(this->element_());

//             if (eIdx == 13)
//             {
//                 if ((this->problem_().coupled() == true) && (phaseIdx == 0))
//                 {
//                     std::cout.precision(15);
//                     std::cout << "potentialGrad_ = " << fluxVars.potentialGrad(phaseIdx) << std::endl;
//                     std::cout << "Keff = " << fluxVars.Keff() << std::endl;
//                     std::cout << "volumeFlux = " << fluxVars.volumeFlux(phaseIdx) << std::endl;
//                     std::cout << std::endl;
//                     std::cout << "flux = " << flux << std::endl;
//                 }
//             }

//                 if (eIdx == 210)
//                 {
//                     std::cout << "fluxVars.volumeFlux(0) = " << fluxVars.volumeFlux(0) << std::endl;
//                 }
//             std::cout << "up.density(phaseIdx) = " << up.density(phaseIdx) << std::endl;
//             std::cout << "dn.density(phaseIdx) = " << dn.density(phaseIdx) << std::endl;
//             std::cout.precision(20);

            if (this->problem_().coupled() == true && phaseIdx == 0) {
                if (!fluxVars.onBoundary())
                {
                    Dune::FieldVector<Scalar,4> lameParams = this->problem_().spatialParams().lameParams(this->element_(), this->fvGeometry_(), 0);

                    Scalar E = lameParams[0];
                    //bulk modulus for both the pure elastic and the viscoelastic model
                    Scalar B = lameParams[1];

                    Scalar lambda = 3.0*B*((3.0*B-E)/(9.0*B-E));
                    Scalar mu = ((3.0*B*E)/(9.0*B-E));

                    Scalar dt = this->problem_().timeManager().timeStepSize();
                    Scalar beta = GET_RUNTIME_PARAM(TypeTag, Scalar, Damping.Beta);

                    Scalar currentGradP = fluxVars.potentialGrad(phaseIdx)* fluxVars.face().normal;
                    Scalar oldGradP = fluxVarsOld.potentialGrad(phaseIdx)* fluxVarsOld.face().normal;
                    Scalar dGradP = currentGradP - oldGradP;

                    Scalar stabilizationTerm = dGradP;
                    stabilizationTerm *= fluxVars.face().area * fluxVars.face().area;
                    stabilizationTerm /= beta * (lambda + 2.0 * mu);

                    stabilizationTerm /= dt;

//                     std::cout << "potentialGrad_(" << phaseIdx << ")      = " << fluxVars.potentialGrad(phaseIdx) << std::endl;
//                     std::cout << "stabilizationTerm of " << phaseIdx << " = " << fluxVars.potentialGrad(phaseIdx) << std::endl;


//                     for (int coordDir = 0; coordDir < dim; ++coordDir)
//                             tmp[coordDir] = fluxVars.potentialGrad(phaseIdx)[coordDir] / dt;

//                     flux[eqIdx] -= stabilizationTerm;
                }
            }
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current sub-control-volume
     *
     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        flux += 0.0;
    }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }
    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

private:
    Scalar massUpwindWeight_;

};

}

#endif
