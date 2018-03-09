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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the three-phase three-component model.
 */
#ifndef DUMUX_3P3C_FLUX_VARIABLES_HH
#define DUMUX_3P3C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the three-phase three-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class ThreePThreeCFluxVariables : public GET_PROP_TYPE(TypeTag, BaseFluxVariables)
{
    friend typename GET_PROP_TYPE(TypeTag, BaseFluxVariables); // be friends with base class
    typedef typename GET_PROP_TYPE(TypeTag, BaseFluxVariables) BaseFluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim>  DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        gCompIdx = Indices::gCompIdx
    };

public:

    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        BaseFluxVariables::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);
        calculateValues_(problem, element, elemVolVars);
    }

private:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Scalar(0);
            molarDensity_[phaseIdx] = Scalar(0);
        }

        // calculate densities at the integration points of the face
        GlobalPosition tmp(0.0);
        for (unsigned int idx = 0; idx < this->face().numFap; idx++) // loop over adjacent vertices
        {
            // index for the element volume variables
            int volVarsIdx = this->face().fapIndices[idx];

            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                density_[phaseIdx] += elemVolVars[volVarsIdx].density(phaseIdx)*
                    this->face().shapeValue[idx];
                molarDensity_[phaseIdx] += elemVolVars[volVarsIdx].molarDensity(phaseIdx)*
                    this->face().shapeValue[idx];
            }
        }

        calculatePorousDiffCoeff_(problem, element, elemVolVars);
    }

    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // initialize to zero
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            massFractionCompWGrad_[phaseIdx] = Scalar(0);
            massFractionCompNGrad_[phaseIdx] = Scalar(0);
            massFractionCompGGrad_[phaseIdx] = Scalar(0);
            moleFractionCompWGrad_[phaseIdx] = Scalar(0);
            moleFractionCompNGrad_[phaseIdx] = Scalar(0);
            moleFractionCompGGrad_[phaseIdx] = Scalar(0);
        }

        BaseFluxVariables::calculateGradients_(problem, element, elemVolVars);

        // calculate gradients
        GlobalPosition tmp(0.0);
        for (unsigned int idx = 0;
             idx < this->face().numFap;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const GlobalPosition &feGrad = this->face().grad[idx];

            // index for the element volume variables
            int volVarsIdx = this->face().fapIndices[idx];

            // the concentration gradient of the components
            // component in the phases
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(wPhaseIdx, wCompIdx);
            massFractionCompWGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(nPhaseIdx, wCompIdx);
            massFractionCompWGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(gPhaseIdx, wCompIdx);
            massFractionCompWGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(wPhaseIdx, nCompIdx);
            massFractionCompNGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(nPhaseIdx, nCompIdx);
            massFractionCompNGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(gPhaseIdx, nCompIdx);
            massFractionCompNGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(wPhaseIdx, gCompIdx);
            massFractionCompGGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(nPhaseIdx, gCompIdx);
            massFractionCompGGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].massFraction(gPhaseIdx, gCompIdx);
            massFractionCompGGrad_[gPhaseIdx] += tmp;

            // the molar concentration gradients of the components
            // in the phases
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, wCompIdx);
            moleFractionCompWGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(nPhaseIdx, wCompIdx);
            moleFractionCompWGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(gPhaseIdx, wCompIdx);
            moleFractionCompWGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, nCompIdx);
            moleFractionCompNGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(nPhaseIdx, nCompIdx);
            moleFractionCompNGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(gPhaseIdx, nCompIdx);
            moleFractionCompNGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, gCompIdx);
            moleFractionCompGGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(nPhaseIdx, gCompIdx);
            moleFractionCompGGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(gPhaseIdx, gCompIdx);
            moleFractionCompGGrad_[gPhaseIdx] += tmp;
        }
    }

    void calculatePorousDiffCoeff_(const Problem &problem,
                               const Element &element,
                               const ElementVolumeVariables &elemVolVars)
    {

        const VolumeVariables &volVarsI = elemVolVars[this->face().i];
        const VolumeVariables &volVarsJ = elemVolVars[this->face().j];

        // the effective diffusion coefficients at vertex i and j
        Scalar diffCoeffI;
        Scalar diffCoeffJ;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to calculate only diffusion coefficients
            // for phases which exist in both finite volumes
            /* \todo take care: This should be discussed once again
             * as long as a meaningful value can be found for the required mole fraction
             * diffusion should work even without this one here */
            if (volVarsI.saturation(phaseIdx) <= 0 ||
                volVarsJ.saturation(phaseIdx) <= 0)
            {
                porousDiffCoeff_[phaseIdx][wCompIdx] = 0.0;
                porousDiffCoeff_[phaseIdx][nCompIdx] = 0.0;
                porousDiffCoeff_[phaseIdx][gCompIdx] = 0.0;
                continue;
            }

            // Diffusion coefficient in the porous medium
            diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                         volVarsI.saturation(phaseIdx),
                                                                         volVarsI.diffCoeff(phaseIdx, wCompIdx));
            diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                         volVarsJ.saturation(phaseIdx),
                                                                         volVarsJ.diffCoeff(phaseIdx, wCompIdx));

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx][wCompIdx] = harmonicMean(diffCoeffI, diffCoeffJ);

            // Diffusion coefficient in the porous medium
            diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                         volVarsI.saturation(phaseIdx),
                                                                         volVarsI.diffCoeff(phaseIdx, nCompIdx));
            diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                         volVarsJ.saturation(phaseIdx),
                                                                         volVarsJ.diffCoeff(phaseIdx, nCompIdx));

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx][nCompIdx] = harmonicMean(diffCoeffI, diffCoeffJ);

            // Diffusion coefficient in the porous medium
            diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                         volVarsI.saturation(phaseIdx),
                                                                         volVarsI.diffCoeff(phaseIdx, gCompIdx));
            diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                         volVarsJ.saturation(phaseIdx),
                                                                         volVarsJ.diffCoeff(phaseIdx, gCompIdx));

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx][gCompIdx] = harmonicMean(diffCoeffI, diffCoeffJ);
        }
    }

public:
    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     *
     *   \param phaseIdx The phase index
     *   \param compIdx The component index
     */
    Scalar porousDiffCoeff(int phaseIdx, int compIdx) const
    { return porousDiffCoeff_[phaseIdx][compIdx];}

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief The mass fraction gradient of the water in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &massFractionCompWGrad(int phaseIdx) const
    {return massFractionCompWGrad_[phaseIdx];}

    /*!
     * \brief The mass fraction gradient of the contaminant in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &massFractionCompNGrad(int phaseIdx) const
    { return massFractionCompNGrad_[phaseIdx]; };

    /*!
     * \brief The mass fraction gradient of gas in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &massFractionCompGGrad(int phaseIdx) const
    { return massFractionCompGGrad_[phaseIdx]; };

    /*!
     * \brief The mole fraction gradient of the water in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &moleFractionCompWGrad(int phaseIdx) const
    { return moleFractionCompWGrad_[phaseIdx]; };

    /*!
     * \brief The mole fraction gradient of the contaminant in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &moleFractionCompNGrad(int phaseIdx) const
    { return moleFractionCompNGrad_[phaseIdx]; };

    /*!
     * \brief The mole fraction gradient of gas in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &moleFractionCompGGrad(int phaseIdx) const
    { return moleFractionCompGGrad_[phaseIdx]; };

protected:
    // gradients
    GlobalPosition massFractionCompWGrad_[numPhases];
    GlobalPosition massFractionCompNGrad_[numPhases];
    GlobalPosition massFractionCompGGrad_[numPhases];
    GlobalPosition moleFractionCompWGrad_[numPhases];
    GlobalPosition moleFractionCompNGrad_[numPhases];
    GlobalPosition moleFractionCompGGrad_[numPhases];

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // the diffusivity matrix for the porous medium
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff_;
};

} // end namespace

#endif
