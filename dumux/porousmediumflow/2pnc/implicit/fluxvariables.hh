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
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase two-component model fully implicit model.
 */
#ifndef DUMUX_2PNC_FLUX_VARIABLES_HH
#define DUMUX_2PNC_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase n-component fully implicit model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */

template <class TypeTag>
class TwoPNCFluxVariables : public GET_PROP_TYPE(TypeTag, BaseFluxVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseFluxVariables) BaseFluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld,
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
          };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
//     typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<CoordScalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
            wPhaseIdx = FluidSystem::wPhaseIdx,
            nPhaseIdx = FluidSystem::nPhaseIdx,
            wCompIdx  = FluidSystem::wCompIdx,
         };

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the sub-control-volume face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    TwoPNCFluxVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int fIdx,
                     const ElementVolumeVariables &elemVolVars,
                     const bool onBoundary = false)
    : BaseFluxVariables(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Scalar(0);
            molarDensity_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                massFractionGrad_[phaseIdx][compIdx] = Scalar(0);
                moleFractionGrad_[phaseIdx][compIdx] = Scalar(0);
            }
        }
        calculateGradients_(problem, element, elemVolVars);
        calculateVelocities_(problem, element, elemVolVars);
        calculateporousDiffCoeff_(problem, element, elemVolVars);
    };

protected:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // calculate gradients
        DimVector tmp(0.0);
        for (int idx = 0;
             idx < this->fvGeometry_.numScv;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector &feGrad = face().grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemVolVars[idx].pressure(phaseIdx); //FE grad times phase pressure
                potentialGrad_[phaseIdx] += tmp;
            }

            // the concentration gradient of the non-wetting
            // component in the wetting phase

            for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if(compIdx != phaseIdx) //No grad is needed for this case
                    {
                        tmp = feGrad;
                        tmp *= elemVolVars[idx].massFraction(phaseIdx, compIdx);
                        massFractionGrad_[phaseIdx][compIdx] += tmp;

                        tmp = feGrad;
                        tmp *= elemVolVars[idx].moleFraction(phaseIdx, compIdx);
                        moleFractionGrad_[phaseIdx][compIdx] += tmp;
                    }
                }
            }
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            int i = face().i;
            int j = face().j;
            Scalar fI = rhoFactor_(phaseIdx, i, elemVolVars);
            Scalar fJ = rhoFactor_(phaseIdx, j, elemVolVars);
            if (fI + fJ <= 0)
                fI = fJ = 0.5; // doesn't matter because no phase is
                               // present in both cells!
            density_[phaseIdx] =
                (fI*elemVolVars[i].density(phaseIdx) +
                 fJ*elemVolVars[j].density(phaseIdx))
                /
                (fI + fJ);
            // phase density
            molarDensity_[phaseIdx]
                =
                (fI*elemVolVars[i].molarDensity(phaseIdx) +
                 fJ*elemVolVars[j].molarDensity(phaseIdx))
                /
                (fI + fJ); //arithmetic averaging

            tmp = problem.gravity();
            tmp *= density_[phaseIdx];

            potentialGrad_[phaseIdx] -= tmp;
        }
    }

    DUNE_DEPRECATED_MSG("This method will be removed without replacement!")
    Scalar rhoFactor_(int phaseIdx, int scvIdx, const ElementVolumeVariables &vDat)
    {

        static const Scalar eps = 1e-2;
        const Scalar sat = vDat[scvIdx].density(phaseIdx);
        if (sat > eps)
            return 0.5;
        if (sat <= 0)
            return 0;

        static const Dumux::Spline<Scalar> sp(0, eps, // x0, x1
                                              0, 0.5, // y0, y1
                                              0, 0); // m0, m1
        return sp.eval(sat);
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        const SpatialParams &spatialParams = problem.spatialParams();
        // multiply the pressure potential with the intrinsic
        // permeability
        DimMatrix K(0.0);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            auto K_i = spatialParams.intrinsicPermeability(element,this->fvGeometry_,face().i);
            //K_i *= volVarsI.permFactor();

            auto K_j = spatialParams.intrinsicPermeability(element,this->fvGeometry_,face().j);
            //K_j *= volVarsJ.permFactor();

            spatialParams.meanK(K,K_i,K_j);

            K.mv(potentialGrad_[phaseIdx], Kmvp_[phaseIdx]);
            KmvpNormal_[phaseIdx] = - (Kmvp_[phaseIdx] * face().normal);
        }

        // set the upstream and downstream vertices
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            upstreamIdx_[phaseIdx] = face().i;
            downstreamIdx_[phaseIdx] = face().j;

            if (KmvpNormal_[phaseIdx] < 0) {
                std::swap(upstreamIdx_[phaseIdx],
                          downstreamIdx_[phaseIdx]);
            }
        }
    }

    void calculateporousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            /* If there is no phase saturation on either side of the face
                * no diffusion takes place */

            if (volVarsI.saturation(phaseIdx) <= 0 ||
                volVarsJ.saturation(phaseIdx) <= 0)
                {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        porousDiffCoeff_[phaseIdx][compIdx] = 0.0;
                    }
                }

            else
            {
            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient
            Scalar tauI =  1.0/(volVarsI.porosity() * volVarsI.porosity()) *
                            pow(volVarsI.porosity() * volVarsI.saturation(phaseIdx), 7.0/3);

            Scalar tauJ =   1.0/(volVarsJ.porosity() * volVarsJ.porosity()) *
                            pow(volVarsJ.porosity() * volVarsJ.saturation(phaseIdx), 7.0/3);
            // Diffusion coefficient in the porous medium

            // -> harmonic mean
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if(phaseIdx==compIdx)
                        porousDiffCoeff_[phaseIdx][compIdx] = 0.0;
                    else
                    {
                        porousDiffCoeff_[phaseIdx][compIdx] = harmonicMean(volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * volVarsI.diffCoeff(phaseIdx, compIdx),
                                                                            volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * volVarsJ.diffCoeff(phaseIdx, compIdx));
                    }
                }
            }
        }
    }

public:
    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability which goes from vertex i to
     *        vertex j.
     *
     * Note that the length of the face's normal is the area of the
     * face, so this is not the actual velocity by the integral of
     * the velocity over the face's area. Also note that the phase
     * mobility is not yet included here since this would require a
     * decision on the upwinding approach (which is done in the
     * model and/or local residual file).
     *
     *   \param phaseIdx The phase index
     */
    Scalar KmvpNormal(int phaseIdx) const
    { return KmvpNormal_[phaseIdx]; }

    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability as vector (for velocity output)
     *
     *   \param phaseIdx The phase index
     */
    DimVector Kmvp(int phaseIdx) const
    { return Kmvp_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     *   \param phaseIdx The phase index
     */
    int upstreamIdx(int phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     *   \param phaseIdx The phase index
     */
    int downstreamIdx(int phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     *
     *   \param phaseIdx The phase index
     *   \param compIdx The component index
     */
    Scalar porousDiffCoeff(int phaseIdx, int compIdx) const
    { return porousDiffCoeff_[phaseIdx][compIdx];}

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief The concentration gradient of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    const DimVector &massFractionGrad(int phaseIdx, int compIdx) const
    { return massFractionGrad_[phaseIdx][compIdx]; }

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    const DimVector &moleFractionGrad(int phaseIdx, int compIdx) const
    { return moleFractionGrad_[phaseIdx][compIdx]; }

    const SCVFace &face() const
    {
    if (this->onBoundary_)
        return this->fvGeometry_.boundaryFace[this->faceIdx_];
    else
        return this->fvGeometry_.subContVolFace[this->faceIdx_];
    }

protected:

    // gradients
    DimVector potentialGrad_[numPhases];
    DimVector massFractionGrad_[numPhases][numComponents];
    DimVector moleFractionGrad_[numPhases][numComponents];

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // intrinsic permeability times pressure potential gradient
    DimVector Kmvp_[numPhases];
    // projected on the face normal
    Scalar KmvpNormal_[numPhases];

    // local index of the upwind vertex for each phase
    int upstreamIdx_[numPhases];
    // local index of the downwind vertex for each phase
    int downstreamIdx_[numPhases];

    // the diffusion coefficient for the porous medium
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff_;
};

} // end namespace

#endif
