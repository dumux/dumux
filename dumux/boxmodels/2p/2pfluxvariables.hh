// $Id$
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

    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;

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
    TwoPFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &elemGeom,
                 int faceIdx,
                 const ElementVolumeVariables &elemDat)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = faceIdx;

        for (int phase = 0; phase < numPhases; ++phase) {
            potentialGrad_[phase] = Scalar(0);
        }

        calculateGradients_(problem, element, elemDat);
        calculateK_(problem, element, elemDat);
    };

public:
    /*
     * \brief Return the intrinsic permeability.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param normalFlux The normal flux i.e. the given intrinsic permeability
     *                   times the pressure potential gradient and SCV face normal.
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().j:face().i; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param normalFlux The normal flux i.e. the given intrinsic permeability
     *                   times the pressure potential gradient and SCV face normal.
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux > 0)?face().i:face().j; }

    /*!
     * \brief Return the SCV (sub-control-volume) face
    */
    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

protected:
    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_[numPhases];

    // intrinsic permeability
    Tensor K_;

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        // calculate gradients
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
                Vector tmp(feGrad);
                tmp *= elemDat[idx].pressure(phase);
                potentialGrad_[phase] += tmp;
            }
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the phase density at the integration
            // point. for this we make sure only existing phases
            // matter
            Scalar SI = elemDat[face().i].saturation(phaseIdx);
            Scalar SJ = elemDat[face().j].saturation(phaseIdx);
            Scalar rhoI = elemDat[face().i].density(phaseIdx);
            Scalar rhoJ = elemDat[face().j].density(phaseIdx);
            Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
            Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
            if (fI + fJ <= 0.0)
                fI = fJ = 0.5; // doesn't matter because no phase is
                               // present in both cells!
            Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

            Vector tmp(problem.gravity());
            tmp *= density;

            potentialGrad_[phaseIdx] -= tmp;
        }
    }

    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemDat)
    {
        const SpatialParameters &spatialParams = problem.spatialParameters();
        // calculate the intrinsic permeability
        spatialParams.meanK(K_,
                            spatialParams.intrinsicPermeability(element,
                                                                fvElemGeom_,
                                                                face().i),
                            spatialParams.intrinsicPermeability(element,
                                                                fvElemGeom_,
                                                                face().j));
    }
};

} // end namepace

#endif
