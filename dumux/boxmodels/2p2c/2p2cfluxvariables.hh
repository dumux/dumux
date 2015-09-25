// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
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
/*!
 * \ingroup TwoPTwoCModel
 */
#ifndef DUMUX_2P2C_FLUX_VARIABLES_HH
#define DUMUX_2P2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

/*!
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,
    };

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemDat The volume variables of the current element
     */
    TwoPTwoCFluxVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int faceIdx,
                     const ElementVolumeVariables &elemDat)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = faceIdx;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            densityAtIP_[phaseIdx] = Scalar(0);
            molarDensityAtIP_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
            concentrationGrad_[phaseIdx] = Scalar(0);
            molarConcGrad_[phaseIdx] = Scalar(0);
        }

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
    };

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        // calculate gradients
        Vector tmp(0.0);
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = face().grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().massFrac(lPhaseIdx, gCompIdx);
            concentrationGrad_[lPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().moleFrac(lPhaseIdx, gCompIdx);
            molarConcGrad_[lPhaseIdx] += tmp;

//            // the concentration gradient of the wetting component
//            // in the non-wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().massFrac(gPhaseIdx, lCompIdx);
            concentrationGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().moleFrac(gPhaseIdx, lCompIdx);
            molarConcGrad_[gPhaseIdx] += tmp;
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            int i = face().i;
            int j = face().j;
            Scalar fI = rhoFactor_(phaseIdx, i, elemDat);
            Scalar fJ = rhoFactor_(phaseIdx, j, elemDat);
            if (fI + fJ <= 0)
                fI = fJ = 0.5; // doesn't matter because no phase is
                               // present in both cells!
            densityAtIP_[phaseIdx] =
                (fI*elemDat[i].density(phaseIdx) +
                 fJ*elemDat[j].density(phaseIdx))
                /
                (fI + fJ);
            // phase density
            molarDensityAtIP_[phaseIdx]
                =
                (fI*elemDat[i].molarDensity(phaseIdx) +
                 fJ*elemDat[j].molarDensity(phaseIdx))
                /
                (fI + fJ);

            tmp = problem.gravity();
            tmp *= densityAtIP_[phaseIdx];

            potentialGrad_[phaseIdx] -= tmp;
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
                              const ElementVolumeVariables &elemDat)
    {
        const SpatialParameters &spatialParams = problem.spatialParameters();
        // multiply the pressure potential with the intrinsic
        // permeability
        Tensor K;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvElemGeom_,
                                                                    face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvElemGeom_,
                                                                    face().j));
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

    void calculateDiffCoeffPM_(const Problem &problem,
                               const Element &element,
                               const ElementVolumeVariables &elemDat)
    {
        const VolumeVariables &vDat_i = elemDat[face().i];
        const VolumeVariables &vDat_j = elemDat[face().j];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to only calculate diffusion coefficents
            // for phases which exist in both finite volumes
            if (vDat_i.saturation(phaseIdx) <= 0 ||
                vDat_j.saturation(phaseIdx) <= 0)
            {
                porousDiffCoeff_[phaseIdx] = 0.0;
                continue;
            }

            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient

            Scalar tau_i =
                1.0/(vDat_i.porosity() * vDat_i.porosity()) *
                pow(vDat_i.porosity() * vDat_i.saturation(phaseIdx), 7.0/3);
            Scalar tau_j =
                1.0/(vDat_j.porosity() * vDat_j.porosity()) *
                pow(vDat_j.porosity() * vDat_j.saturation(phaseIdx), 7.0/3);
            // Diffusion coefficient in the porous medium

            // -> harmonic mean
            porousDiffCoeff_[phaseIdx] = harmonicMean(vDat_i.porosity() * vDat_i.saturation(phaseIdx) * tau_i * vDat_i.diffCoeff(phaseIdx),
                                         vDat_j.porosity() * vDat_j.saturation(phaseIdx) * tau_j * vDat_j.diffCoeff(phaseIdx));

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
    Vector Kmvp(int phaseIdx) const
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
    Scalar densityAtIP(int phaseIdx) const
    { return densityAtIP_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensityAtIP(int phaseIdx) const
    { return molarDensityAtIP_[phaseIdx]; }

    /*!
     * \brief The concentration gradient of a component in a phase.
     */
    const Vector &concentrationGrad(int phaseIdx) const
    { return concentrationGrad_[phaseIdx]; };

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     */
    const Vector &molarConcGrad(int phaseIdx) const
    { return molarConcGrad_[phaseIdx]; };

    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

protected:
    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_[numPhases];
    Vector concentrationGrad_[numPhases];
    Vector molarConcGrad_[numPhases];

    // density of each face at the integration point
    Scalar densityAtIP_[numPhases], molarDensityAtIP_[numPhases];

    // intrinsic permeability times pressure potential gradient
    Vector Kmvp_[numPhases];
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
