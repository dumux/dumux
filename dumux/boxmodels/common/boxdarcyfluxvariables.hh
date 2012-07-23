// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_BOX_DARCY_FLUX_VARIABLES_HH
#define DUMUX_BOX_DARCY_FLUX_VARIABLES_HH

#include "boxproperties.hh"

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

namespace Dumux
{
    
namespace Properties
{
// forward declaration of properties 
NEW_PROP_TAG(MobilityUpwindWeight);
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(EnableGravity);
}   

/*!
 * \ingroup BoxModel
 * \ingroup BoxFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity 
 * on a (sub)control volume face.
 */
template <class TypeTag>
class BoxDarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

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
     * \param elemDat The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    BoxDarcyFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int faceIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        : fvGeometry_(fvGeometry), faceIdx_(faceIdx), onBoundary_(onBoundary)
    {
        mobilityUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MobilityUpwindWeight);
        
        calculateNormalVelocity_(problem, element, elemVolVars);
    }

public:
    /*!
     * \brief Return the normal velocity for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    Scalar normalVelocity(int phaseIdx) const
    { return normalVelocity_[phaseIdx]; }
    
    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    int downstreamIdx(int phaseIdx) const
    { return (normalVelocity_[phaseIdx] >= 0) ? face().j : face().i; }
    
    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    int upstreamIdx(int phaseIdx) const
    { return (normalVelocity_[phaseIdx] > 0) ? face().i : face().j; }

    /*!
     * \brief Return the SCV (sub-control-volume) face
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
    int faceIdx_;
    const bool onBoundary_;
    Scalar normalVelocity_[numPhases];
    Scalar mobilityUpwindWeight_;

private:
    void calculateNormalVelocity_(const Problem &problem,
                                  const Element &element,
                                  const ElementVolumeVariables &elemVolVars)
    {
        // calculate the mean intrinsic permeability
        const SpatialParams &spatialParams = problem.spatialParams();
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(element,
                                                                fvGeometry_,
                                                                face().i),
                            spatialParams.intrinsicPermeability(element,
                                                                fvGeometry_,
                                                                face().j));
        
        // loop over all phases 
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the phase pressure gradient
            Vector gradPotential(0.0);         
            for (int idx = 0;
                 idx < fvGeometry_.numFAP;
                 idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const Vector &feGrad = face().grad[idx];

                // index for the element volume variables 
                int volVarsIdx = face().fapIndices[idx];
	    
                // the pressure gradient
                Vector tmp(feGrad);
                tmp *= elemVolVars[volVarsIdx].pressure(phaseIdx);
                gradPotential += tmp;
            }

            // correct the pressure gradient by the gravitational acceleration
            if (GET_PARAM(TypeTag, bool, EnableGravity))
            {
                // estimate the gravitational acceleration at a given SCV face
                // using the arithmetic mean
                Vector g(problem.boxGravity(element, fvGeometry_, face().i));
                g += problem.boxGravity(element, fvGeometry_, face().j);
                g /= 2;

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
                Vector f(g);
                f *= density;

                // calculate the final potential gradient
                gradPotential -= f;
            } // gravity

            // calculate the flux in the normal direction of the
            // current sub control volume face:
            //
            // v = - (K grad phi) * n
            //
            // (the minus comes from the Darcy law which states that
            // the flux is from high to low potentials.)
            Vector kGradPotential;
            K.mv(gradPotential, kGradPotential);
            Scalar normalFlux = -(kGradPotential*face().normal);

            // determine the upwind direction
            int upstreamIdx, downstreamIdx;
            if (normalFlux > 0) 
            {
                upstreamIdx = face().i;
                downstreamIdx = face().j;
            }
            else 
            {
                upstreamIdx = face().j;
                downstreamIdx = face().i;
            }
            
            // obtain the upwind volume variables
            const VolumeVariables& upVolVars = elemVolVars[upstreamIdx];
            const VolumeVariables& downVolVars = elemVolVars[downstreamIdx];
                
            // set the normal velocity
            normalVelocity_[phaseIdx] = normalFlux
                  *( mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx));
        } // loop over all phases
    }
};

} // end namepace

#endif
