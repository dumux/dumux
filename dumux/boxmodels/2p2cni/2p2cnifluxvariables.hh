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
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

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
                            const FVElementGeometry &fvGeometry,
                            const int faceIdx,
                            const ElementVolumeVariables &elemVolVars,
                            bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {
        faceIdx_ = faceIdx;

        calculateValues_(problem, element, elemVolVars);
    }

    /*!
     * \brief The total heat flux \f$\mathrm{[J/s]}\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

    DimVector temperatureGradient() const
    { return temperatureGrad_; }

protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        temperatureGrad_ = 0;
        DimVector tmp(0.0);
        for (int idx = 0; idx < this->fvGeometry_.numFAP; idx++)
        {
            tmp = this->face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = this->face().fapIndices[idx];

            tmp *= elemVolVars[volVarsIdx].temperature();
            temperatureGrad_ += tmp;
        }

        // The spatial parameters calculates the actual heat flux vector
        if (this->face().i != this->face().j)
            problem.spatialParams().matrixHeatFlux(tmp,
                                                   *this,
                                                   elemVolVars,
                                                   temperatureGrad_,
                                                   element,
                                                   this->fvGeometry_,
                                                   faceIdx_);
        else // heat flux at outflow boundaries
            problem.spatialParams().boundaryMatrixHeatFlux(tmp,
                                                           *this,
                                                           elemVolVars,
                                                           this->face(),
                                                           element,
                                                           this->fvGeometry_);

        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = tmp*this->face().normal;
    }

private:
    Scalar normalMatrixHeatFlux_;
    DimVector temperatureGrad_;
    int faceIdx_;
};

} // end namespace

#endif
