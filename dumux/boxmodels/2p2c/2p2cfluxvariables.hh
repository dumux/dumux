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
 * \brief   This file contains the data which is required to calculate
 *          all fluxes of components over a face of a finite volume for
 *          the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_FLUX_VARIABLES_HH
#define DUMUX_2P2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "2p2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
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
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
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

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
    * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Distinguishes if we are on a SCV face or on a boundary face
     */
    TwoPTwoCFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const int faceIdx,
                          const ElementVolumeVariables &elemVolVars,
                          const bool onBoundary = false)
    : BaseFluxVariables(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Scalar(0);
            molarDensity_[phaseIdx] = Scalar(0);
            massFractionGrad_[phaseIdx] = Scalar(0);
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
        DimVector tmp(0.0);
        for (int idx = 0;
             idx < this->fvGeometry_.numFAP;
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
        DimVector tmp(0.0);
        for (int idx = 0;
             idx < this->fvGeometry_.numFAP;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector &feGrad = this->face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = this->face().fapIndices[idx];

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(wPhaseIdx, nCompIdx);
            massFractionGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(wPhaseIdx, nCompIdx);
            moleFractionGrad_[wPhaseIdx] += tmp;

            //            // the concentration gradient of the wetting component
            //            // in the non-wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().massFraction(nPhaseIdx, wCompIdx);
            massFractionGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(nPhaseIdx, wCompIdx);
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

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to only calculate diffusion coefficients
            // for phases which exist in both finite volumes
            if (volVarsI.saturation(phaseIdx) <= 0 ||
                volVarsJ.saturation(phaseIdx) <= 0)
            {
                porousDiffCoeff_[phaseIdx] = 0.0;
                continue;
            }

            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient
            Scalar tauI =
                1.0/(volVarsI.porosity() * volVarsI.porosity()) *
                pow(volVarsI.porosity() * volVarsI.saturation(phaseIdx), 7.0/3);
            Scalar tauJ =
                1.0/(volVarsJ.porosity() * volVarsJ.porosity()) *
                pow(volVarsJ.porosity() * volVarsJ.saturation(phaseIdx), 7.0/3);
            // Diffusion coefficient in the porous medium

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx] = harmonicMean(volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * volVarsI.diffCoeff(phaseIdx),
                                                      volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * volVarsJ.diffCoeff(phaseIdx));
        }
    }

public:
    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability which goes from vertex i to
     *        vertex j.
     *
     * Note that the length of the face's normal is the area of the
     * phase, so this is not the actual velocity by the integral of
     * the velocity over the face's area. Also note that the phase
     * mobility is not yet included here since this would require a
     * decision on the upwinding approach (which is done in the
     * actual model).
     */
    DUNE_DEPRECATED_MSG("use the methods of the base flux variables instead")
    Scalar KmvpNormal(int phaseIdx) const
    { return -this->kGradPNormal(phaseIdx); }

    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability as vector (for velocity output)
     */
    DUNE_DEPRECATED_MSG("use the methods of the base flux variables instead")
    DimVector Kmvp(int phaseIdx) const
    { return this->kGradP_[phaseIdx]; }

    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     */
    Scalar porousDiffCoeff(int phaseIdx) const
    { return porousDiffCoeff_[phaseIdx]; };

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    DUNE_DEPRECATED_MSG("use density instead")
    Scalar densityAtIP(int phaseIdx) const
    { return density(phaseIdx); }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase.
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    DUNE_DEPRECATED_MSG("use molarDensity instead")
    Scalar molarDensityAtIP(int phaseIdx) const
    { return molarDensity(phaseIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase.
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief The concentration gradient of a component in a phase.
     */
    DUNE_DEPRECATED_MSG("use massFractionGrad instead")
    const DimVector &concentrationGrad(int phaseIdx) const
    { return massFractionGrad(phaseIdx); };

    /*!
     * \brief The mass fraction gradient of the dissolved component in a phase.
     */
    DUNE_DEPRECATED_MSG("use moleFractionGrad instead")
    const DimVector &massFractionGrad(int phaseIdx) const
    { return massFractionGrad_[phaseIdx]; };

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     */
    DUNE_DEPRECATED_MSG("use moleFractionGrad instead")
    const DimVector &molarConcGrad(int phaseIdx) const
    { return moleFractionGrad(phaseIdx); };

    /*!
     * \brief The mole fraction gradient of the dissolved component in a phase.
     */
    const DimVector &moleFractionGrad(int phaseIdx) const
    { return moleFractionGrad_[phaseIdx]; };

protected:
    // gradients
    DimVector massFractionGrad_[numPhases];
    DimVector moleFractionGrad_[numPhases];

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // intrinsic permeability times pressure potential gradient
    DimVector Kmvp_[numPhases];

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_[numPhases];
};

} // end namepace

#endif
