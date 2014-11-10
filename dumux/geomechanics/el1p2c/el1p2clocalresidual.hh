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
 * \brief Element-wise calculation the local Jacobian for the linear elastic,
 * single-phase, two-component model in the fully implicit scheme.
 */
#ifndef DUMUX_ELASTIC1P2C_LOCAL_RESIDUAL_HH
#define DUMUX_ELASTIC1P2C_LOCAL_RESIDUAL_HH

#include "el1p2cproperties.hh"

namespace Dumux
{
    /*!
     * \ingroup ElOnePTwoCModel
     * \ingroup ImplicitLocalResidual
     * \brief Calculate the local Jacobian for a one-phase two-component
     * flow in a linear-elastic porous medium.
     *
     * This class is used to fill the gaps in BoxLocalResidual for the
     * one-phase two-component linear elasticity model.
     */
    template<class TypeTag>
    class ElOnePTwoCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
    {
        protected:
        typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

        enum { dim = GridView::dimension };
        typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
        typedef Dune::FieldVector<Scalar, dim> DimVector;

        typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
        typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
        typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
        typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

        enum {
            //phase index
            phaseIdx = Indices::phaseIdx,
            transportCompIdx = Indices::transportCompIdx
            };
        // indices of the equations
        enum {
            conti0EqIdx = Indices::conti0EqIdx,
            transportEqIdx = Indices::transportEqIdx
        };

        //! property that defines whether mole or mass fractions are used
        static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);


        public:
            /*!
             * \brief Constructor. Sets the upwind weight.
             */
            ElOnePTwoCLocalResidual()
            {
                // retrieve the upwind weight for the mass conservation equations. Use the value
                // specified via the property system as default, and overwrite
                // it by the run-time parameter from the Dune::ParameterTree
                upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
                // retrieve the property which defines if the stabilization terms in the mass balance
                // equations are switched on. Use the value specified via the property system as default,
                // and overwrite it by the run-time parameter from the Dune::ParameterTree
                withStabilization_ =  GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, WithStabilization);
            };

            /*!
             * \brief Evaluate the amount of all conservation quantities
             *        (e.g. phase mass) within a finite volume.
             *
             *        \param storage The mass of the component within the sub-control volume
             *        \param scvIdx The index of the considered face of the sub-control volume
             *        \param usePrevSol Evaluate function with solution of current or previous time step
             */
            void computeStorage(PrimaryVariables &storage, int scvIdx,
                            bool usePrevSol) const
            {
                // if flag usePrevSol is set, the solution from the previous
                // time step is used, otherwise the current solution is
                // used. The secondary variables are used accordingly.  This
                // is required to compute the derivative of the storage term
                // using the implicit euler method.
                const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[scvIdx];

                storage = 0;
                 // this model assumes incompressible fluids and solids only the bulk
                 // density i.e. the ratio of solid phase and pore fluid can vary
                 // storage term of continuity equation
                storage[conti0EqIdx] += volVars.divU;

                if(!useMoles)
                {
                    //storage term of the transport equation - massfractions
                    storage[transportEqIdx] += volVars.fluidState().massFraction(phaseIdx, transportCompIdx)
                            * volVars.effPorosity;
                }
                else
                {
                    // storage term of the transport equation - molefractions
                    storage[transportEqIdx] += volVars.fluidState().moleFraction(phaseIdx, transportCompIdx)
                        * volVars.effPorosity;
                }
            }

            /*!
             * \brief Evaluate the mass flux over a face of a sub-control
             *        volume.
             *
             *        \param flux The flux over the SCV (sub-control-volume) face for each component
             *        \param faceIdx The index of the considered face of the sub control volume
             *        \param onBoundary A boolean variable to specify whether the flux variables
             *               are calculated for interior SCV faces or boundary faces, default=false
             */
            void computeFlux(PrimaryVariables &flux, int faceIdx, const bool onBoundary=false) const
            {
                flux = 0;
                FluxVariables fluxVars(this->problem_(),
                                this->element_(),
                                this->fvGeometry_(),
                                faceIdx,
                                this->curVolVars_());

                this->computeAdvectiveFlux(flux, fluxVars);
                this->computeDiffusiveFlux(flux, fluxVars);
                this->computeStresses(flux, fluxVars, faceIdx);
            }

            /*!
             * \brief Evaluates the advective mass flux of all phases over
             *        a face of a subcontrol volume.
             *
             * \param flux The advective flux over the sub-control-volume face for each component
             * \param fluxVars The flux variables at the current SCV
             */
            void computeAdvectiveFlux(PrimaryVariables &flux,
                                      const FluxVariables &fluxVars) const
            {
                ////////
                // advective fluxes of all components in all phases
                ////////

                // data attached to upstream and the downstream vertices
                // of the current phase
                const VolumeVariables &up =
                    this->curVolVars_(fluxVars.upstreamIdx());
                const VolumeVariables &dn =
                    this->curVolVars_(fluxVars.downstreamIdx());

                // calculate the stabilization term which helps in case of stability problems
                // e.g. observed for small time steps (according to G.Aguilar, F.Gaspar, F.Lisbona
                // and C.Rodrigo (2008))
                Scalar stabilizationTerm(0.0);
                if(withStabilization_){
                // calculate distance h between nodes i and j
                DimVector hVec = this->element_().geometry().corner(fluxVars.face().j)
                                      - this->element_().geometry().corner(fluxVars.face().i);
                Scalar h = hVec.two_norm();
                stabilizationTerm = (h * h) /
                                    (4 * (fluxVars.lambda()
                                    + 2 * fluxVars.mu()));

                stabilizationTerm *= fluxVars.timeDerivGradPNormal();
                }


                // flux for mass balance of the solid-fluid mixture
                // KmvpNormal is the Darcy velocity multiplied with the normal vector,
                // calculated in 1p2cfluxvariables.hh
                flux[conti0EqIdx] +=
                   fluxVars.KmvpNormal() *
                   ((     upwindWeight_)/up.viscosity()
                    +
                   ((1 - upwindWeight_)/dn.viscosity()));


                // stabilization term
                if(withStabilization_)
                flux[conti0EqIdx] -= stabilizationTerm;

                if(!useMoles)
                {
                // mass flux of the dissolved second component - massfraction
                // advective flux of the component
                flux[transportEqIdx] +=
                    fluxVars.KmvpNormal() *
                    ((    upwindWeight_)* up.fluidState().massFraction(phaseIdx, transportCompIdx)/up.viscosity()
                     +
                     (1 - upwindWeight_)* dn.fluidState().massFraction(phaseIdx, transportCompIdx)/dn.viscosity());

                // flux of the dissolved second component due to solid displacement
                flux[transportEqIdx] +=
                     fluxVars.timeDerivUNormal() *
                     ((    upwindWeight_)* up.fluidState().massFraction(phaseIdx, transportCompIdx)
                     * up.effPorosity
                     +
                     (1 - upwindWeight_)*dn.fluidState().massFraction(phaseIdx, transportCompIdx)
                     * up.effPorosity);

                // stabilization term
                if(withStabilization_)
                flux[transportEqIdx] -=
                    stabilizationTerm *
                    ((    upwindWeight_)* up.fluidState().massFraction(phaseIdx, transportCompIdx)
                     +
                     (1 - upwindWeight_)*dn.fluidState().massFraction(phaseIdx, transportCompIdx));
                }
                else
                {
                                // mass flux of the dissolved second component - massfraction
                // advective flux of the component
                flux[transportEqIdx] +=
                    fluxVars.KmvpNormal() *
                    ((    upwindWeight_)* up.fluidState().moleFraction(phaseIdx, transportCompIdx)/up.viscosity()
                     +
                     (1 - upwindWeight_)* dn.fluidState().moleFraction(phaseIdx, transportCompIdx)/dn.viscosity());

                // flux of the dissolved second component due to solid displacement
                flux[transportEqIdx] +=
                     fluxVars.timeDerivUNormal() *
                     ((    upwindWeight_)* up.fluidState().moleFraction(phaseIdx, transportCompIdx)
                     * up.effPorosity
                     +
                     (1 - upwindWeight_)*dn.fluidState().moleFraction(phaseIdx, transportCompIdx)
                     * up.effPorosity);

                // stabilization term
                if(withStabilization_)
                flux[transportEqIdx] -=
                    stabilizationTerm *
                    ((    upwindWeight_)* up.fluidState().moleFraction(phaseIdx, transportCompIdx)
                     +
                     (1 - upwindWeight_)*dn.fluidState().moleFraction(phaseIdx, transportCompIdx));
                }

            }

            /*!
             * \brief Adds the diffusive mass flux of all components over
             *        a face of a sub-control volume.
             *
             * \param flux The diffusive flux over the sub-control-volume face for each component
             * \param fluxVars The flux variables at the current sub-control-volume face
             */
            void computeDiffusiveFlux(PrimaryVariables &flux,
                                      const FluxVariables &fluxVars) const
            {
                Scalar tmp(0);

                // diffusive flux of second component
                if(!useMoles)
                {
                    // diffusive flux of the second component - massfraction
                    tmp = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal);
                    tmp *= fluxVars.diffCoeffPM();

                    // dispersive flux of second component - massfraction
                               DimVector normalDisp;
                               fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
                               tmp -= (normalDisp * fluxVars.moleFractionGrad(transportCompIdx));

                    // convert it to a mass flux and add it
                    flux[transportEqIdx] += tmp * FluidSystem::molarMass(transportCompIdx);
                }
                else
                {
                                    // diffusive flux of the second component - molefraction
                    tmp = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal);
                    tmp *= fluxVars.diffCoeffPM();

                    // dispersive flux of second component - molefraction
                                DimVector normalDisp;
                                fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
                                tmp -= (normalDisp * fluxVars.moleFractionGrad(transportCompIdx));

                    flux[transportEqIdx] += tmp;
                }
            }

            /*!
             * \brief Evaluates the total stress induced by effective stresses and fluid
             * pressure in the solid fluid mixture.
             * \param stress The stress over the sub-control-volume face for each component
             * \param fluxVars The variables at the current sub-control-volume face
             * \param faceIdx The index of the current sub-control-volume face
             */
            void computeStresses(PrimaryVariables &stress,
                    const FluxVariables &fluxVars, const int faceIdx) const
            {
                DimMatrix pressure(0.0), sigma(0.0);
                // the normal vector corresponding to the current sub-control-volume face
                const DimVector &normal(this->fvGeometry_().subContVolFace[faceIdx].normal);

                // the pressure term of the momentum balance
                for (int i = 0; i < dim; ++i)
                    pressure[i][i] += 1.0;

                pressure *= fluxVars.pressure();
                // effective stresses
                sigma = fluxVars.sigma();
                // calculate total stresses by subtracting the pressure
                sigma -= pressure;

                DimVector tmp(0.0);
                // multiply total stress tensor with normal vector of current face
                sigma.mv(normal, tmp);

                // set the stress term equal to the calculated vector
                for (int i = 0; i < dim; ++i)
                stress[Indices::momentum(i)] = tmp[i];
            }

            /*!
             * \brief Calculate the source term of the equation
             *        \param source The source/sink in the SCV for each component
             *        \param scvIdx The index of the vertex of the sub control volume
             *
             */
            void computeSource(PrimaryVariables &source, const int scvIdx)
            {
                source = 0;

                const ElementVolumeVariables &elemVolVars = this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[scvIdx];

                DimVector tmp1(0.0), tmp2(0.0);

                this->problem_().solDependentSource(source,
                                             this->element_(),
                                             this->fvGeometry_(),
                                             scvIdx,
                                             this->curVolVars_());

                // the gravity term of the momentum balance equation is treated as a source term
                // gravity contribution of solid matrix
                tmp1 = this->problem_().gravity();
                tmp1 *= volVars.rockDensity();
                tmp1 *= (1. - volVars.effPorosity);

                // gravity contribution of the fluids
                tmp2 = this->problem_().gravity();
                tmp2 *= volVars.density();
                tmp2 *= volVars.effPorosity;

                tmp1 += tmp2;

                for (int i = 0; i < dim; ++i)
                    source[Indices::momentum(i)] += tmp1[i];
            }

            Implementation *asImp_()
            { return static_cast<Implementation *> (this); }
            const Implementation *asImp_() const
            { return static_cast<const Implementation *> (this); }

        private:
            Scalar upwindWeight_;
            bool withStabilization_;
    };

}

#endif
