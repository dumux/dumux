// $Id: 2pfluxdata.hh 3736 2010-06-15 09:52:10Z lauser $
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_2P_FLUX_DATA_HH
#define DUMUX_2P_FLUX_DATA_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPBoxModel
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the two-phase model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoPFluxData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))    Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;

    typedef typename GridView::ctype                     CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef std::vector<VertexData>                      VertexDataArray;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef Dune::FieldVector<CoordScalar, dimWorld>  Vector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename FVElementGeometry::SubControlVolume             SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace         SCVFace;

    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

public:
    TwoPFluxData(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &elemGeom,
                 int faceIdx,
                 const VertexDataArray &elemDat)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = faceIdx;

        for (int phase = 0; phase < numPhases; ++phase) {
            densityAtIP_[phase] = Scalar(0);
            potentialGrad_[phase] = Scalar(0);
        }

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
    };

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
     * \brief Return density [kg/m^3] of a phase at the integration
     *        point.
     */
    Scalar densityAtIP(int phaseIdx) const
    { return densityAtIP_[phaseIdx]; }

    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

protected:
    const FVElementGeometry &fvElemGeom_;
    int                      scvfIdx_;

    // gradients
    Vector potentialGrad_[numPhases];
    Vector concentrationGrad_[numPhases];

    // density of each face at the integration point
    Scalar densityAtIP_[numPhases];

    // intrinsic permeability times pressure potential gradient
    // projected on the face normal
    Scalar KmvpNormal_[numPhases];

    // local index of the upwind vertex for each phase
    int upstreamIdx_[numPhases];
    // local index of the downwind vertex for each phase
    int downstreamIdx_[numPhases];

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_[numPhases];


private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const VertexDataArray &elemDat)
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
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure(phase);
                potentialGrad_[phase] += tmp;

                // phase density
                densityAtIP_[phase]
                    +=
                    elemDat[idx].density(phase) *
                    face().shapeValue[idx];
            }
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
        {
            tmp = problem.gravity();
            tmp *= densityAtIP_[phase];

            potentialGrad_[phase] -= tmp;
        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const VertexDataArray &elemDat)
    {
        const SpatialParameters &spatialParams = problem.spatialParameters();
        // multiply the pressure potential with the intrinsic
        // permeability
        Tensor K;
        Vector Kmvp;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            spatialParams.meanK(K,
                   spatialParams.intrinsicPermeability(element,
                                                       fvElemGeom_,
                                                       face().i),
                    spatialParams.intrinsicPermeability(element,
                                                        fvElemGeom_,
                                                        face().j));
            K.mv(potentialGrad_[phaseIdx], Kmvp);
            KmvpNormal_[phaseIdx] = - (Kmvp * face().normal);
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
};

} // end namepace

#endif
