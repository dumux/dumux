// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 *        the flux of fluid over a face of a finite volume.
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \brief This template class contains the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        dimWorld = GridView::dimensionworld,
        
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename GridView::template Codim<0>::Entity Element;

public:
    RichardsFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int faceIdx,
                          const ElementVolumeVariables &elemVolVars)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = faceIdx;

        calculateGradients_(problem, element, elemVolVars);
        calculateK_(problem, element, elemVolVars);
    };

    /*
     * \brief Return the intrinsic permeability.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient.
     */
    const Vector &potentialGradW() const
    { return potentialGrad_; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().j:face().i; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux > 0)?face().i:face().j; }

    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

protected:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        potentialGrad_ = 0.0;
        // calculate gradients
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex index
            const Vector &feGrad = face().grad[idx];

            // the pressure gradient
            Vector tmp(feGrad);
            tmp *= elemVolVars[idx].pressure(wPhaseIdx);
            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        ///////////////

        // calculate the phase density at the integration point. we
        // only do this if the wetting phase is present in both cells
        Scalar SI = elemVolVars[face().i].saturation(wPhaseIdx);
        Scalar SJ = elemVolVars[face().j].saturation(wPhaseIdx);
        Scalar rhoI = elemVolVars[face().i].density(wPhaseIdx);
        Scalar rhoJ = elemVolVars[face().j].density(wPhaseIdx);
        Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
        Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
        if (fI + fJ == 0)
            // doesn't matter because no wetting phase is present in
            // both cells!
            fI = fJ = 0.5; 
        Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

        Vector tmp(problem.gravity());
        tmp *= density;
        potentialGrad_ -= tmp;
    }

    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemDat)
    {
        const SpatialParameters &sp = problem.spatialParameters();
        // calculate the intrinsic permeability
        sp.meanK(K_,
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().i),
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().j));
    }

    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_;

    // intrinsic permeability
    Tensor K_;
};

} // end namepace

#endif
