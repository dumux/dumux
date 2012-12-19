// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes (mass and energy) of all phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_2PNI_FLUX_VARIABLES_HH
#define DUMUX_2PNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes (mass and energy) of all phases over a
 *        face of a finite volume for the non-isothermal two-phase model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPNIFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
{
    typedef ImplicitDarcyFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;

public:

    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */

    TwoPNIFluxVariables(const Problem &problem,
                   const Element &element,
                   const FVElementGeometry &fvGeometry,
                   int faceIdx,
                   const ElementVolumeVariables &elemVolVars,
                   const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {
        // calculate temperature gradient using finite element
        // gradients
        Vector temperatureGrad(0);
        for (int idx = 0; idx < fvGeometry.numFAP; idx++)
        {
            Vector feGrad = this->face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = this->face().fapIndices[idx];
            
            feGrad *= elemVolVars[volVarsIdx].temperature();
            temperatureGrad += feGrad;
        }

        // The spatial parameters calculates the actual heat flux vector
        Vector heatFlux;
        problem.spatialParams().matrixHeatFlux(heatFlux,
                                                   *this,
                                                   elemVolVars,
                                                   temperatureGrad,
                                                   element,
                                                   fvGeometry,
                                                   faceIdx);
        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = heatFlux * this->face().normal;
    }

    /*!
     * \brief The total heat flux \f$\mathrm{[J/s]}\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

private:
    Scalar normalMatrixHeatFlux_;
};

} // end namepace

#endif
