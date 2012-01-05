// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Katherina Baber
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
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the flux of the fluid over a face of a
 *        finite volume for the one-phase model.
 */
template <class TypeTag>
class OnePFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

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

        calculateK_(problem, element);
        calculateGradients_(problem, element, elemDat);
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
                             const ElementVolumeVariables &elemVolVars)
    {
        potentialGrad_ = 0.0;

        // calculate potential gradient
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = face().grad[idx];

            // the pressure gradient
            Vector tmp(feGrad);
            tmp *= elemVolVars[idx].pressure();
            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(problem.boxGravity(element, fvElemGeom_, face().i));
            g += problem.boxGravity(element, fvElemGeom_, face().j);
            g /= 2;

            // calculate the phase density at the integration point. we
            // only do this if the wetting phase is present in both cells
            Scalar rhoI = elemVolVars[face().i].density();
            Scalar rhoJ = elemVolVars[face().j].density();
            Scalar density = (rhoI + rhoJ)/2;

            // make it a force
            Vector f(g);
            f *= density;

            // calculate the final potential gradient
            potentialGrad_ -= f;
        }
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
