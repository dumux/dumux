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
 * \brief This file contains data which is required to calculate
 *        the heat fluxes over a face of a finite volume.
 *
 * This means temperature gradients and the normal matrix
 * heat flux.
 */
#ifndef DUMUX_3P3CNI_FLUX_VARIABLES_HH
#define DUMUX_3P3CNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/3p3c/3p3cfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePThreeCNIModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the heat fluxes over a face of a finite
 *        volume for the non-isothermal three-phase three-component model.
 *        The mass fluxes are computed in the parent class.
 *
 * This means temperature gradients and the normal matrix
 * heat flux.
 */
template <class TypeTag>
class ThreePThreeCNIFluxVariables : public ThreePThreeCFluxVariables<TypeTag>
{
    typedef ThreePThreeCFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<CoordScalar, dim> DimVector;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Distinguishes if we are on a sub-control-volume face or on a boundary face
     */
    ThreePThreeCNIFluxVariables(const Problem &problem,
                                const Element &element,
                                const FVElementGeometry &fvGeometry,
                                const int faceIdx,
                                const ElementVolumeVariables &elemVolVars,
                                const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {
        // calculate temperature gradient using finite element
        // gradients
        GlobalPosition temperatureGrad(0);
        GlobalPosition tmp(0.0);
        for (unsigned int idx = 0; idx < this->face().numFap; idx++)
        {
            tmp = this->face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = this->face().fapIndices[idx];

            tmp *= elemVolVars[volVarsIdx].temperature();
            temperatureGrad += tmp;
        }

        // The spatial parameters calculates the actual heat flux vector
        problem.spatialParams().matrixHeatFlux(tmp,
                                                   *this,
                                                   elemVolVars,
                                                   temperatureGrad,
                                                   element,
                                                   fvGeometry,
                                                   faceIdx);
        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = tmp*this->face().normal;
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

} // end namespace

#endif
