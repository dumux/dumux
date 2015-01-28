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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume in the
 *        two phase discrete fracture-matrix model.
 */
#ifndef DUMUX_MODELS_2PDFM_FLUX_VARIABLES_HH
#define DUMUX_MODELS_2PDFM_FLUX_VARIABLES_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>

#include "2pdfmproperties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPDFMModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate the fluxes of 
 *        the fluid phases over a face of a finite volume for the two-phase
 *        discrete fracture-matrix model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoPDFMFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

public:
    /*!
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
    TwoPDFMFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        : ImplicitDarcyFluxVariables<TypeTag>(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary), 
          fvGeometry_(fvGeometry), faceIdx_(fIdx), onBoundary_(onBoundary)
    {
        faceSCV_ = &this->face();

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            potentialGradFracture_[phaseIdx] = 0.0;
        }

        calculateGradientsInFractures_(problem, element, elemVolVars, fIdx);
        calculateVelocitiesFracture_(problem, element, elemVolVars, fIdx);
    };

public:
    /*!
     * \brief Calculates the velocities in the lower-dimenstional fracture.
     * 
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     * \param fIdx The local index of the SCV (sub-control-volume) face
     */
    void calculateVelocitiesFracture_(const Problem &problem,
                                      const Element &element,
                                      const ElementVolumeVariables &elemVolVars,
                                      int fIdx)
    {
        isFracture_ = problem.spatialParams().isEdgeFracture(element, fIdx);
        fractureWidth_ = problem.spatialParams().fractureWidth(element, fIdx);

        Scalar KFracture, KFi, KFj; //permeabilities of the fracture
        if (isFracture_)
        {
            KFi = problem.spatialParams().
                intrinsicPermeabilityFracture(element, this->fvGeometry_, this->face().i);
            KFj = problem.spatialParams().
                intrinsicPermeabilityFracture(element, this->fvGeometry_, this->face().j);
        }
        else
        {
            KFi = 0;
            KFj = 0;
        }

        KFracture = Dumux::harmonicMean(KFi, KFj);

        // temporary vector for the Darcy velocity
        Scalar vDarcyFracture = 0;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            vDarcyFracture = - (KFracture * potentialGradFracture_[phaseIdx]);

            Valgrind::CheckDefined(KFracture);
            Valgrind::CheckDefined(fractureWidth_);
            vDarcyFracture_[phaseIdx] = (vDarcyFracture * fractureWidth_);
            Valgrind::CheckDefined(vDarcyFracture_[phaseIdx]);
        }

        // set the upstream and downstream vertices
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            upstreamFractureIdx[phaseIdx] = faceSCV_->i;
            downstreamFractureIdx[phaseIdx] = faceSCV_->j;

            if (vDarcyFracture_[phaseIdx] < 0)
            {
                std::swap(upstreamFractureIdx[phaseIdx],
                          downstreamFractureIdx[phaseIdx]);
            }
        }
    }

    /*!
     * \brief Return the pressure potential gradient in the lower dimensional fracture.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Scalar &potentialGradFracture(int phaseIdx) const
    {
        return potentialGradFracture_[phaseIdx];
    }


    PhasesVector vDarcyFracture_;

    int upstreamFractureIdx[numPhases];
    int downstreamFractureIdx[numPhases];
protected:
    // gradients
    Scalar potentialGradFracture_[numPhases];
    const FVElementGeometry &fvGeometry_;
    int faceIdx_;
    const bool onBoundary_;
    bool isFracture_;
    Scalar fractureWidth_;
    const SCVFace *faceSCV_;

private:
    /*!
     * \brief Calculates the gradients in the lower-dimenstional fracture.
     * 
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     * \param fIdx The local index of the SCV (sub-control-volume) face
     */
    void calculateGradientsInFractures_(const Problem &problem,
                                        const Element &element,
                                        const ElementVolumeVariables &elemVolVars,
                                        int fIdx)
    {
        // calculate gradients, loop over adjacent vertices
        for (unsigned int idx = 0; idx < this->face().numFap; idx++)
        {
            int i = this->face().i;
            int j = this->face().j;

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                const GlobalPosition localIdx_i = element.geometry().corner(i);
                const GlobalPosition localIdx_j = element.geometry().corner(j);

                isFracture_ = problem.spatialParams().isEdgeFracture(element, fIdx);

                if (isFracture_)
                {
                    GlobalPosition diff_ij = localIdx_j;
                    diff_ij -= localIdx_i;
                    potentialGradFracture_[phaseIdx] =
                        (elemVolVars[j].pressure(phaseIdx) - elemVolVars[i].pressure(phaseIdx))
                        / diff_ij.two_norm();

                    // correct the pressure gradient by the gravitational acceleration
                    if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
                    {
                        // ask for the gravitational acceleration at the given SCV face
                        GlobalPosition g(problem.gravityAtPos(this->face().ipGlobal));

                        // calculate the phase density at the integration point. we
                        // only do this if the wetting phase is present in both cells
                        Scalar SI = elemVolVars[this->face().i].fluidState().saturation(phaseIdx);
                        Scalar SJ = elemVolVars[this->face().j].fluidState().saturation(phaseIdx);
                        Scalar rhoI = elemVolVars[this->face().i].fluidState().density(phaseIdx);
                        Scalar rhoJ = elemVolVars[this->face().j].fluidState().density(phaseIdx);
                        Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                        Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                        if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(fI + fJ, 0.0, 1.0e-30))
                            // doesn't matter because no wetting phase is present in
                            // both cells!
                            fI = fJ = 0.5;
                        const Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

                        // make gravity acceleration a force
                        GlobalPosition f(g);
                        f *= density;

                        //transform into fracture local coordinates
                        diff_ij/=diff_ij.two_norm();
                        Scalar fractureLocalGravityForce=diff_ij*f;

                        // calculate the final potential gradient
                        potentialGradFracture_[phaseIdx] -= fractureLocalGravityForce;
                    }// gravity
                }
                else
                {
                    potentialGradFracture_[phaseIdx] = 0;
                }
            }
        }
    }
};

} // end namespace

#endif // DUMUX_MODELS_2PDFM_FLUX_VARIABLES_HH
