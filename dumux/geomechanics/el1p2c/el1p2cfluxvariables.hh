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
 * \brief This file contains the calculation of all the fluxes over the surface of the
 * finite volume that make up the volume, the mass and the momentum balance
 * for the one-phase two-component linear-elastic model.
 *
 * This means pressure, concentration and solid-displacement gradients, phase densities at
 * the integration point, etc.
 *
 * This class inherits from the one-phase two-component model FluxVariables and from the
 * linear elasticity model FluxVariables
 */
#ifndef DUMUX_ELASTIC1P2C_FLUX_VARIABLES_HH
#define DUMUX_ELASTIC1P2C_FLUX_VARIABLES_HH

#include <dumux/geomechanics/elastic/elasticfluxvariables.hh>
#include <dumux/porousmediumflow/1p2c/implicit/fluxvariables.hh>

namespace Dumux
{
/*!
 * \ingroup ElOnePTwoCBoxModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes over the surface of the
 *           finite volume that make up the volume, the mass and the momentum balance
 *           for the one-phase two-component linear-elastic model.
 *
 * This means pressure, concentration and solid-displacement gradients, phase densities at
 * the integration point, etc.
 *
 */
    template<class TypeTag>
    class ElOnePTwoCFluxVariables: public ElasticFluxVariablesBase<TypeTag> ,
                                   public OnePTwoCFluxVariables<TypeTag>
    {
            typedef ElasticFluxVariablesBase<TypeTag> ElasticBase;
            typedef OnePTwoCFluxVariables<TypeTag> OnePTwoCBase;

            typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
            typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
            typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
            typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;

            typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
            typedef typename GridView::template Codim<0>::Entity Element;
            enum
            {
                dim = GridView::dimension,
                dimWorld = GridView::dimensionworld
            };

            typedef typename GridView::ctype CoordScalar;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
            typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
            typedef Dune::FieldVector<CoordScalar, dim> DimVector;

            typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
            typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

        public:
        /*
         * \brief The constructor
         *
         * \param problem The problem
         * \param element The finite element
         * \param fvGeometry The finite-volume geometry in the fully implicit scheme
         * \param fIdx The local index of the SCV (sub-control-volume) face
         * \param elemVolVars The volume variables of the current element
         * \param onBoundary A boolean variable to specify whether the flux variables
         * are calculated for interior SCV faces or boundary faces, default=false
         */
            ElOnePTwoCFluxVariables(const Problem &problem,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            int fIdx,
                            const ElementVolumeVariables &elemVolVars,
                            const bool onBoundary = false)
                : ElasticBase(problem, element, fvGeometry, fIdx, elemVolVars),
                  OnePTwoCBase(problem, element, fvGeometry, fIdx, elemVolVars),
                  fvGeometry_(fvGeometry), faceIdx_(fIdx), onBoundary_(onBoundary)
            {
                dU_ = Scalar(0);
                dGradP_ = Scalar(0);
                porosity_ = 0;
                effPorosity_ = 0;
                pressure_ = 0;
                timeDerivUNormal_ = 0;

                elOnePTwoCGradients_(problem, element, elemVolVars);
                calculateEffectiveValues_(problem, element, elemVolVars);
                calculateDiffCoeffPM_(problem, element, elemVolVars);
                calculateDDt_(problem, element, elemVolVars);

            }
            ;

        public:
            /*!
             * \brief Return porosity [-] at the integration point.
             */
            Scalar porosity() const
            {
                return porosity_;
            }

            /*!
             * \brief Return effective porosity [-] at the integration point.
             */
            Scalar effPorosity() const
            {
                return effPorosity_;
            }

            /*!
             * \brief Return pressure [Pa] at the integration
             *        point.
             */
            Scalar pressure() const
            {
                return pressure_;
            }

            /*!
             * \brief Return change of pressure gradient with time [Pa/m] at
             * integration point.
             */
            Scalar dGradP(int dimIdx) const
            {
                return dGradP_[dimIdx];
            }

            /*!
             * \brief Return gradient of time derivative of pressure [Pa].
             */
            Scalar timeDerivGradPNormal() const
            {
                return timeDerivGradPNormal_;
            }

            /*!
             * \brief Return change of u [m] with time at integration point
             *        point.
             */
            Scalar dU(int dimIdx) const
            {
                return dU_[dimIdx];
            }

            /*!
             * \brief Return time derivative of u [m/s] in normal direction at integration point
             */
            Scalar timeDerivUNormal() const
            {
                return timeDerivUNormal_;
            }

            /*!
             * \brief Return porous medium diffusion coefficient [m^2]
             */
            Scalar diffCoeffPM() const
            {
                return diffCoeffPM_;
            }

            const SCVFace &face() const
            {
                return fvGeometry_.subContVolFace[faceIdx_];
            }

        protected:
            /*!
             * \brief Calculation of the solid displacement and pressure gradients.
             *
             *        \param problem The considered problem file
             *        \param element The considered element of the grid
             *        \param elemVolVars The parameters stored in the considered element
             */
            void elOnePTwoCGradients_(const Problem &problem,
                            const Element &element,
                            const ElementVolumeVariables &elemVolVars)
            {
                // calculate gradients
                GlobalPosition tmp(0.0);
                for (int idx = 0; idx < fvGeometry_.numScv; idx++) // loop over adjacent vertices
                {
                    // FE gradient at vertex idx
                    const DimVector &feGrad = face().grad[idx];

                    // the gradient of the temporal pressure change (needed for stabilization term)
                    tmp = feGrad;
                    tmp *= elemVolVars[idx].dPressure();
                    dGradP_ += tmp;

                    // average the pressure at integration point
                    pressure_ += elemVolVars[idx].pressure()
                                    * face().shapeValue[idx];
                    // average temporal displacement change at integration point (for calculation of solid displacement velocity)
                    for (int i = 0; i < dim; ++i)
                    dU_[i] += elemVolVars[idx].dU(i)
                                        * face().shapeValue[idx];
                    // average porosity at integration point
                    porosity_ += elemVolVars[idx].porosity()
                                    * face().shapeValue[idx];
                }
            }

            /*!
             * \brief Calculation of the effective porosity.
             *
             *        \param problem The considered problem file
             *        \param element The considered element of the grid
             *        \param elemVolVars The parameters stored in the considered element
             */
            void calculateEffectiveValues_(const Problem &problem,
                            const Element &element,
                            const ElementVolumeVariables &elemVolVars)
            {

                // the effective porosity is calculated as a function of solid displacement and initial porosity
                // according to Han & Dusseault (2003)

                // calculate effective porosity as a function of solid displacement and initial porosity
                effPorosity_ = (porosity_ + this->divU())
                                  / (1 + this->divU());
            }

            /*!
             * \brief Calculation of the effective porous media diffusion coefficient.
             *
             *        \param problem The considered problem file
             *        \param element The considered element of the grid
             *        \param elemVolVars The parameters stored in the considered element
             */
            void calculateDiffCoeffPM_(const Problem &problem,
                            const Element &element,
                            const ElementVolumeVariables &elemVolVars)
            {
                const VolumeVariables &volVarsI = elemVolVars[face().i];
                const VolumeVariables &volVarsJ = elemVolVars[face().j];

                const Scalar diffCoeffI =
                    EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                    /*sat=*/1.0,
                                                                    volVarsI.diffCoeff());

                const Scalar diffCoeffJ =
                    EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                    /*sat=*/1.0,
                                                                    volVarsJ.diffCoeff());

                diffCoeffPM_ = harmonicMean(diffCoeffI, diffCoeffJ);
            }

            /*!
             * \brief Calculation of the time derivative of solid displacement and pressure gradient
             *        \param problem The considered problem file
             *        \param element The considered element of the grid
             *        \param elemVolVars The parameters stored in the considered element
             */
            void calculateDDt_(const Problem &problem,
                    const Element &element,
                    const ElementVolumeVariables &elemVolVars)
            {
                Scalar dt= problem.timeManager().timeStepSize();
                DimVector tmp(0.0);

                //time derivative of solid displacement times normal vector
                for (int i = 0; i < dim; ++i)
                    tmp[i] = dU_[i] / dt;
                       timeDerivUNormal_ = tmp * face().normal;
                    //time derivative of pressure gradient times normal vector
                for (int i = 0; i < dim; ++i)
                    tmp[i] = dGradP_[i] / dt;
                      timeDerivGradPNormal_ = tmp * face().normal;
            }

            const FVElementGeometry &fvGeometry_;
            int faceIdx_;
            const bool onBoundary_;

            //! change of solid displacement with time at integration point
            GlobalPosition dU_;
            //! change of pressure gradient with time at integration point
            GlobalPosition dGradP_;
            //! porosity at integration point
            Scalar porosity_;
            //! effective porosity at integration point
            Scalar effPorosity_;
            //! pressure at integration point
            Scalar pressure_;
            //! time derivative of solid displacement times normal vector at integration point
            Scalar timeDerivUNormal_;
            //! time derivative of pressure gradient times normal vector at integration point
            Scalar timeDerivGradPNormal_;
            //! Parameters
            Scalar diffCoeffPM_;
    };

} // end namespace

#endif
