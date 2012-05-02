// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
class TwoPTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::lPhaseIdx,
        nPhaseIdx = Indices::gPhaseIdx,
        wCompIdx = Indices::lCompIdx,
        nCompIdx = Indices::gCompIdx
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
        : fvGeometry_(fvGeometry), faceIdx_(faceIdx), onBoundary_(onBoundary)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Scalar(0);
            molarDensity_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
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
             idx < fvGeometry_.numVertices;
             idx++) // loop over adjacent vertices
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                density_[phaseIdx] += elemVolVars[idx].density(phaseIdx)*
                    face().shapeValue[idx];
                molarDensity_[phaseIdx] += elemVolVars[idx].molarDensity(phaseIdx)*
                    face().shapeValue[idx];
            }
        }

        calculateGradients_(problem, element, elemVolVars);
        calculateVelocities_(problem, element, elemVolVars);
        calculatePorousDiffCoeff_(problem, element, elemVolVars);
    }

    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // calculate gradients
        DimVector tmp(0.0);
        for (int idx = 0;
             idx < fvGeometry_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector &feGrad = face().grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemVolVars[idx].pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[idx].fluidState().massFraction(wPhaseIdx, nCompIdx);
            massFractionGrad_[wPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[idx].fluidState().moleFraction(wPhaseIdx, nCompIdx);
            moleFractionGrad_[wPhaseIdx] += tmp;

            //            // the concentration gradient of the wetting component
            //            // in the non-wetting phase
            tmp = feGrad;
            tmp *= elemVolVars[idx].fluidState().massFraction(nPhaseIdx, wCompIdx);
            massFractionGrad_[nPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemVolVars[idx].fluidState().moleFraction(nPhaseIdx, wCompIdx);
            moleFractionGrad_[nPhaseIdx] += tmp;
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            DimVector g(problem.boxGravity(element, fvGeometry_, face().i));
            g += problem.boxGravity(element, fvGeometry_, face().j);
            g /= 2;

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = elemVolVars[face().i].saturation(phaseIdx);
                Scalar SJ = elemVolVars[face().j].saturation(phaseIdx);
                Scalar rhoI = elemVolVars[face().i].density(phaseIdx);
                Scalar rhoJ = elemVolVars[face().j].density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (fI + fJ == 0)
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

                // make gravity acceleration a force
                DimVector f(g);
                f *= density;

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            }
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

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        const SpatialParams &spatialParams = problem.spatialParams();
        // multiply the pressure potential by the intrinsic permeability
        DimMatrix K;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_,
                                                                    face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_,
                                                                    face().j));
            K.mv(potentialGrad_[phaseIdx], Kmvp_[phaseIdx]);
            KmvpNormal_[phaseIdx] = -(Kmvp_[phaseIdx]*face().normal);
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

    void calculatePorousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

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
    Scalar KmvpNormal(int phaseIdx) const
    { return KmvpNormal_[phaseIdx]; }

    /*!
     * \brief Return the pressure potential multiplied with the
     *        intrinsic permeability as vector (for velocity output)
     */
    DimVector Kmvp(int phaseIdx) const
    { return Kmvp_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     */
    int upstreamIdx(int phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     */
    int downstreamIdx(int phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     */
    Scalar porousDiffCoeff(int phaseIdx) const
    { return porousDiffCoeff_[phaseIdx]; };

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    DUMUX_DEPRECATED_MSG("use density instead")
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
    DUMUX_DEPRECATED_MSG("use molarDensity instead")
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
    DUMUX_DEPRECATED_MSG("use massFractionGrad instead")
    const DimVector &concentrationGrad(int phaseIdx) const
    { return massFractionGrad(phaseIdx); };

    /*!
     * \brief The mass fraction gradient of the dissolved component in a phase.
     */
    DUMUX_DEPRECATED_MSG("use moleFractionGrad instead")
    const DimVector &massFractionGrad(int phaseIdx) const
    { return massFractionGrad_[phaseIdx]; };

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     */
    DUMUX_DEPRECATED_MSG("use moleFractionGrad instead")
    const DimVector &molarConcGrad(int phaseIdx) const
    { return moleFractionGrad(phaseIdx); };

    /*!
     * \brief The mole fraction gradient of the dissolved component in a phase.
     */
    const DimVector &moleFractionGrad(int phaseIdx) const
    { return moleFractionGrad_[phaseIdx]; };

    /*!
     * \brief The face of the current sub-control volume. This may be either
     *        an inner sub-control-volume face or a face on the boundary.
     */
    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_.boundaryFace[faceIdx_];
        else
            return fvGeometry_.subContVolFace[faceIdx_];
    }

protected:
    const FVElementGeometry &fvGeometry_;
    const int faceIdx_;
    const bool onBoundary_;

    // gradients
    DimVector potentialGrad_[numPhases];
    DimVector massFractionGrad_[numPhases];
    DimVector moleFractionGrad_[numPhases];

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
    Scalar porousDiffCoeff_[numPhases];
};

} // end namepace

#endif
