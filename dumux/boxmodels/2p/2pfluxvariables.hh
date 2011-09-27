/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_2P_FLUX_VARIABLES_HH
#define DUMUX_2P_FLUX_VARIABLES_HH

#include "2pproperties.hh"

#include <dumux/common/parameters.hh>
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
class TwoPFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVariables)) ElementVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*
     * \brief The constructor
     */
    TwoPFluxVariables()
    {}

#warning Docme
    void update(const ElementVariables &elemVars, int scvfIdx)
    {
        insideScvIdx_ = elemVars.fvElemGeom().subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemVars.fvElemGeom().subContVolFace[scvfIdx].j;

        calculateGradients_(elemVars, scvfIdx);
        calculateNormalFlux_(elemVars, scvfIdx);
    };

    /*!
     * \brief Return the extrusion factor of the SCVF.
     */
    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \brief Return a phase's pressure potential gradient.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return a phase's pressure potential gradient times
     *        intrinsic permeability times the normal of the sub
     *        control volume face times the area of the SCVF.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar normalFlux(int phaseIdx) const
    { return normalFlux_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    int downstreamIdx(int phaseIdx) const
    { return (normalFlux_[phaseIdx] > 0)?outsideScvIdx_:insideScvIdx_; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    int upstreamIdx(int phaseIdx) const
    { return (normalFlux_[phaseIdx] > 0)?insideScvIdx_:outsideScvIdx_; }

protected:
    void calculateGradients_(const ElementVariables &elemVars,
                             int scvfIdx)
    {
        // reset all gradients to 0
        for (int phase = 0; phase < numPhases; ++phase) {
            potentialGrad_[phase] = Scalar(0);
        }
        
        typedef typename FVElementGeometry::SubControlVolumeFace Scvf;
        const Scvf &scvf = elemVars.fvElemGeom().subContVolFace[scvfIdx];

        // calculate gradients
        for (int scvIdx = 0;
             scvIdx < elemVars.numScv();
             scvIdx ++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = scvf.grad[scvIdx];

            // compute sum of pressure gradients for each phase
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
                Vector tmp(feGrad);
                tmp *= elemVars.volVars(scvIdx, /*historyIdx=*/0).pressure(phase);
                potentialGrad_[phase] += tmp;
            }
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(elemVars.problem().boxGravity(elemVars.element(), 
                                                   elemVars.fvElemGeom(), 
                                                   insideScvIdx_));
            g += elemVars.problem().boxGravity(elemVars.element(), 
                                               elemVars.fvElemGeom(), 
                                               outsideScvIdx_);
            g /= 2;
            
            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = elemVars.volVars(insideScvIdx_, /*historyIdx=*/0).saturation(phaseIdx);
                Scalar SJ = elemVars.volVars(outsideScvIdx_, /*historyIdx=*/0).saturation(phaseIdx);
                Scalar rhoI = elemVars.volVars(insideScvIdx_, /*historyIdx=*/0).density(phaseIdx);
                Scalar rhoJ = elemVars.volVars(outsideScvIdx_, /*historyIdx=*/0).density(phaseIdx);
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
                potentialGrad_[phaseIdx] -= f;
            }
        }
    }

    void calculateNormalFlux_(const ElementVariables &elemVars, 
                              int scvfIdx)
    {
        const SpatialParameters &spatialParams = elemVars.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemVars,
                                                                insideScvIdx_),
                            spatialParams.intrinsicPermeability(elemVars,
                                                                outsideScvIdx_));

        const Vector &normal = elemVars.fvElemGeom().subContVolFace[scvfIdx].normal;

        // calculate the flux in the normal direction of the
        // current sub control volume face:
        //
        // v = - (K grad p) * n
        //
        // (the minus comes from the Darcy law which states that
        // the flux is from high to low pressure potentials.)
        Vector tmpVec;
                            
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            K.mv(potentialGrad(phaseIdx), tmpVec);

            // scalar product with the face normal
            normalFlux_[phaseIdx] = 0.0;
            for (int i = 0; i < Vector::size; ++i) 
                normalFlux_[phaseIdx] += tmpVec[i]*normal[i];

            // flux is along negative pressure gradients
            normalFlux_[phaseIdx] *= -1;
        }
    }

    // local indices of the inside and the outside sub-control volumes
    int insideScvIdx_;
    int outsideScvIdx_;

    // gradients
    Vector potentialGrad_[numPhases];

    // normal fluxes
    Scalar normalFlux_[numPhases];
};

} // end namepace

#endif
