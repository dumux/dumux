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
#ifndef DUMUX_2P2C_FLUX_VARIABLES_HH
#define DUMUX_2P2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase two-component model fully implicit model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPTwoCFluxVariables : public GET_PROP_TYPE(TypeTag, BaseFluxVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseFluxVariables) BaseFluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld} ;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

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
    TwoPTwoCFluxVariables(const Problem &problem,
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
            moleFractionGrad_[phaseIdx] = Scalar(0);
        }

        calculateValues_(problem, element, elemVolVars);
    }

 protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        // calculate densities at the integration points of the face
        GlobalPosition tmp(0.0);
        for (unsigned int idx = 0;
             idx < this->face().numFap;
             idx++) // loop over adjacent vertices
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

        calculateGradients_(problem, element, elemVolVars);
        calculatePorousDiffCoeff_(problem, element, elemVolVars);
    }

    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
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

            // the mole fraction gradient of the wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, nCompIdx);
            moleFractionGrad_[wPhaseIdx] += tmp;

            // the mole fraction gradient of the non-wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(nPhaseIdx, wCompIdx);
            moleFractionGrad_[nPhaseIdx] += tmp;
        }
    }

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
            // make sure to only calculate diffusion coefficients
            // for phases which exist in both finite volumes
            if (volVarsI.saturation(phaseIdx) <= 0 || volVarsJ.saturation(phaseIdx) <= 0)
            {
                porousDiffCoeff_[phaseIdx] = 0.0;
                continue;
            }

            diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                         volVarsI.saturation(phaseIdx),
                                                                         volVarsI.diffCoeff(phaseIdx));

            diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                         volVarsJ.saturation(phaseIdx),
                                                                         volVarsJ.diffCoeff(phaseIdx));

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx] = harmonicMean(diffCoeffI, diffCoeffJ);
        }
    }

 public:
    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     *        for each fluid phase in the porous medium.
     *
     * \param phaseIdx The phase index
     */
    Scalar porousDiffCoeff(int phaseIdx) const
    { return porousDiffCoeff_[phaseIdx]; };

    /*!
     * \brief Returns the density \f$\mathrm{[kg/m^3]}\f$ of a phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ of a phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief Returns the mole fraction gradient \f$\mathrm{[1/m]}\f$
     *        of the dissolved component in a phase.
     *
     * \param phaseIdx The phase index
     */
    const GlobalPosition &moleFractionGrad(int phaseIdx) const
    { return moleFractionGrad_[phaseIdx]; };

 protected:
    // mole fraction gradients
    GlobalPosition moleFractionGrad_[numPhases];

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_[numPhases];
};

} // end namespace

#endif
