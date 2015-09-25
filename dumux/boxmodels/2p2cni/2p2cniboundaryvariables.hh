// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of the fluid phase over the boundary of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point of the boundary, etc.
 */
#ifndef DUMUX_2P2CNI_BOUNDARY_VARIABLES_HH
#define DUMUX_2P2CNI_BOUNDARY_VARIABLES_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over the boundary of a
 *        finite volume for the 2p2cni model.
 *
 * This means pressure and velocity gradients, phase density and viscosity at
 * the integration point of the boundary, etc.
 */
template <class TypeTag>
class TwoPTwoCNIBoundaryVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum {
        dim = GridView::dimension,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> VelocityVector;
    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;
    typedef Dune::FieldMatrix<Scalar, dim, dim> VectorGradient;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::BoundaryFace BoundaryFace;

    enum {
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

public:
    TwoPTwoCNIBoundaryVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int boundaryFaceIdx,
                     const ElementVolumeVariables &elemDat,
                     int scvIdx)
        : fvElemGeom_(elemGeom), scvIdx_(scvIdx)
    {
        boundaryFace_ = &fvElemGeom_.boundaryFace[boundaryFaceIdx];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            densityAtIP_[phaseIdx] = Scalar(0);
            molarDensityAtIP_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
            concentrationGrad_[phaseIdx] = Scalar(0);
            molarConcGrad_[phaseIdx] = Scalar(0);
        }

        calculateBoundaryValues_(problem, element, elemDat);
    };

    Scalar KmvpNormal(int phaseIdx) const
    { return KmvpNormal_[phaseIdx]; }

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
    const ScalarGradient &concentrationGrad(int phaseIdx) const
    { return concentrationGrad_[phaseIdx]; };

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     */
    const ScalarGradient &molarConcGrad(int phaseIdx) const
    { return molarConcGrad_[phaseIdx]; };

    /*!
     * \brief The total heat flux \f$\mathrm{[J/s]}\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

    const BoundaryFace& boundaryFace() const
    { return *boundaryFace_; }

protected:
    void calculateBoundaryValues_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        ScalarGradient tmp(0.0);

        // calculate gradients and secondary variables at IPs of the boundary
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const ScalarGradient& feGrad = boundaryFace_->grad[idx];

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

            // the concentration gradient of the wetting component
            // in the non-wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().massFrac(gPhaseIdx, lCompIdx);
            concentrationGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().moleFrac(gPhaseIdx, lCompIdx);
            molarConcGrad_[gPhaseIdx] += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].temperature();
            temperatureGrad_ += tmp;

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                densityAtIP_[phaseIdx] += elemDat[idx].density(phaseIdx)*boundaryFace_->shapeValue[idx];
                molarDensityAtIP_[phaseIdx] += elemDat[idx].molarDensity(phaseIdx)*boundaryFace_->shapeValue[idx];
            }
        }

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            tmp = problem.gravity();
            tmp *= densityAtIP_[phaseIdx];

            potentialGrad_[phaseIdx] -= tmp;

            Scalar k = problem.spatialParameters().intrinsicPermeability(element, fvElemGeom_, scvIdx_);
            VectorGradient K(0);
            K[0][0] = K[1][1] = k;
            ScalarGradient Kmvp;
            K.mv(potentialGrad_[phaseIdx], Kmvp);
            KmvpNormal_[phaseIdx] = - (Kmvp*boundaryFace_->normal);

            const VolumeVariables &vertDat = elemDat[scvIdx_];

            if (vertDat.saturation(phaseIdx) <= 0)
                porousDiffCoeff_[phaseIdx] = 0.0;
            else
            {
                Scalar tau = 1.0/(vertDat.porosity()*vertDat.porosity())*
                    pow(vertDat.porosity()*vertDat.saturation(phaseIdx), 7.0/3);

                porousDiffCoeff_[phaseIdx] = vertDat.porosity()*vertDat.saturation(phaseIdx)*tau*vertDat.diffCoeff(phaseIdx);
            }
        }

        // The spatial parameters calculates the actual heat flux vector
        problem.spatialParameters().boundaryMatrixHeatFlux(tmp,
                                                   elemDat,
                                                   temperatureGrad_,
                                                   element,
                                                   fvElemGeom_,
                                                   scvIdx_);
        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = tmp*boundaryFace_->normal;

    }

    const FVElementGeometry &fvElemGeom_;
    const BoundaryFace *boundaryFace_;

    // gradients
    ScalarGradient potentialGrad_[numPhases];
    ScalarGradient concentrationGrad_[numPhases];
    ScalarGradient molarConcGrad_[numPhases];

    // density of each face at the integration point
    Scalar densityAtIP_[numPhases];
    Scalar molarDensityAtIP_[numPhases];

    // intrinsic permeability times pressure potential gradient
    // projected on the face normal
    Scalar KmvpNormal_[numPhases];

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_[numPhases];

    ScalarGradient temperatureGrad_;
    Scalar normalMatrixHeatFlux_;

    int scvIdx_;
};

} // end namespace

#endif
