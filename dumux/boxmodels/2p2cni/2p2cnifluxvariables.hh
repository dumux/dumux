// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 *        all fluxes (mass of components and energy) over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_2P2CNI_FLUX_VARIABLES_HH
#define DUMUX_2P2CNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/boxmodels/2p2c/2p2cfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes (mass of components and energy) over a face of a finite
 *        volume for the non-isothermal two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPTwoCNIFluxVariables : public TwoPTwoCFluxVariables<TypeTag>
{
    typedef TwoPTwoCFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     */
    TwoPTwoCNIFluxVariables(const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int faceIdx,
                       const ElementVolumeVariables &elemVolVars,
                       bool onBoundary = false)
        : ParentType(problem, element, elemGeom, faceIdx, elemVolVars, onBoundary)
    {
        scvfIdx_ = faceIdx;

        if (!onBoundary)
            calculateValues_(problem, element, this->face(), elemVolVars);
    }

    /*!
     * \brief The total heat flux \f$\mathrm{[J/s]}\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

    Vector temperatureGradient() const
    { return temperatureGrad_; }

protected:
    template<class FaceType>
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const FaceType &face,
                          const ElementVolumeVariables &elemVolVars,
                          bool onBoundary = false)
    {
        if (onBoundary)
            ParentType::calculateValues_(problem, element, face, elemVolVars, onBoundary);

        // calculate temperature gradient using finite element
        // gradients
        temperatureGrad_ = 0;
        Vector tmp(0.0);
        for (int vertIdx = 0; vertIdx < this->fvGeom_.numVertices; vertIdx++)
        {
            tmp = face.grad[vertIdx];
            tmp *= elemVolVars[vertIdx].temperature();
            temperatureGrad_ += tmp;
        }

        // The spatial parameters calculates the actual heat flux vector
        if (!onBoundary)
            problem.spatialParameters().matrixHeatFlux(tmp,
                                                       *this,
                                                       elemVolVars,
                                                       temperatureGrad_,
                                                       element,
                                                       this->fvGeom_,
                                                       scvfIdx_);
        else
            problem.spatialParameters().matrixHeatFlux(tmp,
                                                       *this,
                                                       elemVolVars,
                                                       temperatureGrad_,
                                                       element,
                                                       this->fvGeom_,
                                                       scvfIdx_);

        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = tmp*face.normal;
    }

private:
    Scalar normalMatrixHeatFlux_;
    Vector temperatureGrad_;
    int scvfIdx_;
};

} // end namespace

#endif
