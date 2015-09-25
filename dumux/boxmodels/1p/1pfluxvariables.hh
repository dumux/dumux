// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Onur Dogan                                   *
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
 *        the flux of the fluid over a face of a finite volume for the one-phase model.
 *
 *        This means pressure and temperature gradients, phase densities at
 *           the integration point, etc.
 */
#ifndef DUMUX_1P_FLUX_VARIABLES_HH
#define DUMUX_1P_FLUX_VARIABLES_HH

#include "1pproperties.hh"

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \brief This template class contains the data which is required to
 *        calculate the flux of the fluid over a face of a
 *        finite volume for the one-phase model.
 */
template <class TypeTag>
class OnePFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

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
    OnePFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &elemGeom,
                 int faceIdx,
                 const ElementVolumeVariables &elemDat)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = faceIdx;

        calculateGradients_(problem, element, elemDat);
        calculateK_(problem, element);
    };


    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

    /*!
     * \brief Return the intrinsic permeability.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient.
     */
    const Vector &potentialGrad() const
    { return potentialGrad_; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The normal flux i.e. the given intrinsic permeability
     *                   times the pressure potential gradient and SCV face normal.
     *
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().i:face().j; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The normal flux i.e. the given intrinsic permeability
     *                   times the pressure potential gradient and SCV face normal.
     *
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().j:face().i; }

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        potentialGrad_ = 0.0;
        Scalar densityAtIP = 0.0;

        // calculate potential gradient
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const LocalPosition &feGrad = face().grad[idx];

            // the pressure gradient
            Vector tmp(feGrad);
            tmp *= elemDat[idx].pressure();
            potentialGrad_ += tmp;

            // fluid density
            densityAtIP +=
                elemDat[idx].density()*face().shapeValue[idx];
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        Vector tmp(problem.gravity());
        tmp *= densityAtIP;

        potentialGrad_ -= tmp;
    }

    void calculateK_(const Problem &problem,
                     const Element &element)
    {
        const SpatialParameters &spatialParams = problem.spatialParameters();
        spatialParams.meanK(K_,
                            spatialParams.intrinsicPermeability(element,
                                                                fvElemGeom_,
                                                                face().i),
                            spatialParams.intrinsicPermeability(element,
                                                                fvElemGeom_,
                                                                face().j));
    }

protected:
    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_;

    // intrinsic permeability
    Tensor K_;

    // local index of the upwind vertex
    int upstreamIdx_;
    // local index of the downwind vertex
    int downstreamIdx_;
};

} // end namepace

#endif
