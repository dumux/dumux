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
    friend typename GET_PROP_TYPE(TypeTag, BaseFluxVariables); // be friends with base class
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
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
            wPhaseIdx = FluidSystem::wPhaseIdx,
            nPhaseIdx = FluidSystem::nPhaseIdx,
            wCompIdx  = FluidSystem::wCompIdx,
         };

public:
    /*!
     * \brief The old constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the sub-control-volume face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    DUNE_DEPRECATED_MSG("FluxVariables now have to be default constructed and updated.")
    TwoPNCFluxVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int fIdx,
                     const ElementVolumeVariables &elemVolVars,
                     const bool onBoundary = false)
    : BaseFluxVariables(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary) {}

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    TwoPNCFluxVariables() = default;

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
        calculatePorousDiffCoeff_(problem, element, elemVolVars);
    }

protected:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        BaseFluxVariables::calculateGradients_(problem, element, elemVolVars);

        // initialize to mole/mass fraction gradients to zero
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                moleFractionGrad_[phaseIdx][compIdx] = 0.0; // deprecated
                massFractionGrad_[phaseIdx][compIdx] = 0.0; // deprecated
            }
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            concentrationGrad_[phaseIdx].fill(GlobalPosition(0.0));

         // loop over number of flux approximation points
        for (unsigned int idx = 0; idx < this->face().numFap; ++idx)
        {
            // FE gradient at vertex idx
            const GlobalPosition &feGrad = this->face().grad[idx];

            // index for the element volume variables
            auto volVarsIdx = this->face().fapIndices[idx];

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if(compIdx != phaseIdx) //No grad is needed for this case
                    {
                        GlobalPosition tmp(feGrad); // deprecated
                        tmp *= elemVolVars[volVarsIdx].massFraction(phaseIdx, compIdx);
                        massFractionGrad_[phaseIdx][compIdx] += tmp;

                        tmp = feGrad; // deprecated
                        tmp *= elemVolVars[volVarsIdx].moleFraction(phaseIdx, compIdx);
                        moleFractionGrad_[phaseIdx][compIdx] += tmp;

                        tmp = feGrad;
                        tmp *= elemVolVars[volVarsIdx].moleFraction(phaseIdx, compIdx)*elemVolVars[volVarsIdx].molarDensity(phaseIdx);
                        concentrationGrad_[phaseIdx][compIdx] += tmp;
                    }
                }
            }
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

        static const Spline<Scalar> sp(0, eps, // x0, x1
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
            /* If there is no phase saturation on either side of the face
                * no diffusion takes place */
            if (volVarsI.saturation(phaseIdx) <= 0 || volVarsJ.saturation(phaseIdx) <= 0)
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
                                std::pow(volVarsI.porosity() * volVarsI.saturation(phaseIdx), 7.0/3);

                Scalar tauJ =   1.0/(volVarsJ.porosity() * volVarsJ.porosity()) *
                                std::pow(volVarsJ.porosity() * volVarsJ.saturation(phaseIdx), 7.0/3);
                // Diffusion coefficient in the porous medium

                // -> harmonic mean
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if(phaseIdx==compIdx)
                    {
                        porousDiffCoeff_[phaseIdx][compIdx] = 0.0;
                    }
                    else
                    {
                        auto porousDiffI = volVarsI.porosity() * volVarsI.saturation(phaseIdx) * tauI * volVarsI.diffCoeff(phaseIdx, compIdx);
                        auto porousDiffJ = volVarsJ.porosity() * volVarsJ.saturation(phaseIdx) * tauJ * volVarsJ.diffCoeff(phaseIdx, compIdx);
                        porousDiffCoeff_[phaseIdx][compIdx] = harmonicMean(porousDiffI, porousDiffJ);
                    }
                }
            }
        }
    }

    DUNE_DEPRECATED_MSG("Use calculatePorousDiffCoeff_ (captial P)")
    void calculateporousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    { calculatePorousDiffCoeff_(problem, element, elemVolVars); }

public:
    DUNE_DEPRECATED_MSG("Use darcy flux variables interface.")
    Scalar KmvpNormal(int phaseIdx) const
    { return this->kGradPNormal_[phaseIdx]; }

    DUNE_DEPRECATED_MSG("Will be removed without replacement. Use darcy flux variables interface.")
    GlobalPosition Kmvp(int phaseIdx) const
    { return this->kGradP_[phaseIdx]; }

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
    DUNE_DEPRECATED_MSG("Use concentrationGrad!")
    const GlobalPosition &massFractionGrad(int phaseIdx, int compIdx) const
    { return massFractionGrad_[phaseIdx][compIdx]; }

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    DUNE_DEPRECATED_MSG("Use concentrationGrad!")
    const GlobalPosition &moleFractionGrad(int phaseIdx, int compIdx) const
    { return moleFractionGrad_[phaseIdx][compIdx]; }

    /*!
     * \brief The concentration gradient of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    const GlobalPosition &concentrationGrad(int phaseIdx, int compIdx) const
    { return concentrationGrad_[phaseIdx][compIdx]; }

protected:

    // gradients
    GlobalPosition massFractionGrad_[numPhases][numComponents]; // deprecated
    GlobalPosition moleFractionGrad_[numPhases][numComponents]; // deprecated
    std::array<std::array<GlobalPosition, numComponents>, numPhases> concentrationGrad_;

    // density of each face at the integration point
    Scalar density_[numPhases], molarDensity_[numPhases];

    // the diffusion coefficient for the porous medium
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff_;
};

} // end namespace

#endif
