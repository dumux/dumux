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
 * \brief flux variables for the implicit non-isothermal models
 */
#ifndef DUMUX_NI_FLUX_VARIABLES_HH
#define DUMUX_NI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup NIModel
 * \ingroup ImplicitFluxVariables
 * \brief This class contains data which is required to
 *        calculate the heat fluxes over a face of a finite
 *        volume for the implicit non-isothermal models.
 *        The mass fluxes are computed in the parent class.
 */
template <class TypeTag>
class NIFluxVariables : public GET_PROP_TYPE(TypeTag, IsothermalFluxVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, IsothermalFluxVariables) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dimWorld = GridView::dimensionworld };
    enum { dim = GridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the sub-control-volume face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Distinguishes if we are on a sub-control-volume face or on a boundary face
     */
    NIFluxVariables(const Problem &problem,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const int fIdx,
                            const ElementVolumeVariables &elemVolVars,
                            bool onBoundary = false)
    : ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary)
    {
        calculateValues_(problem, element, elemVolVars);
    }

    /*!
     * \brief The total heat flux \f$\mathrm{[J/s]}\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

    /*!
     * \brief The local temperature gradient at the IP of the considered scv face.
     */
    GlobalPosition temperatureGradient() const
    { return temperatureGrad_; }

    /*!
     * \brief The harmonically averaged effective thermal conductivity.
     */
    Scalar effThermalConductivity() const
    { return lambdaEff_; }

protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        temperatureGrad_ = 0;
        GlobalPosition tmp(0.0);
        for (unsigned int idx = 0; idx < this->face().numFap; idx++)
        {
            tmp = this->face().grad[idx];

            // index for the element volume variables
            int volVarsIdx = this->face().fapIndices[idx];

            tmp *= elemVolVars[volVarsIdx].temperature();
            temperatureGrad_ += tmp;
        }

        lambdaEff_ = 0;
        calculateEffThermalConductivity_(problem, element, elemVolVars);

        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = temperatureGrad_* this->face().normal;
        normalMatrixHeatFlux_ *= -lambdaEff_;
    }

    void calculateEffThermalConductivity_(const Problem &problem,
                                          const Element &element,
                                          const ElementVolumeVariables &elemVolVars)
    {
        const unsigned i = this->face().i;
        const unsigned j = this->face().j;
        Scalar lambdaI, lambdaJ;

        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            lambdaI =
              ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i],
                                                                     problem.spatialParams(),
                                                                     element, this->fvGeometry_, i);

            lambdaJ =
              ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j],
                                                                     problem.spatialParams(),
                                                                     element, this->fvGeometry_, j);
        }
        else
        {
            const Element& elementI = this->fvGeometry_.neighbors[i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            lambdaI =
              ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i],
                                                                     problem.spatialParams(),
                                                                     elementI, fvGeometryI, 0);

            const Element& elementJ = this->fvGeometry_.neighbors[j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            lambdaJ =
              ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j],
                                                                     problem.spatialParams(),
                                                                     elementJ, fvGeometryJ, 0);
        }

        lambdaEff_ = harmonicMean(lambdaI, lambdaJ);
    }

private:
    Scalar lambdaEff_;
    Scalar normalMatrixHeatFlux_;
    GlobalPosition temperatureGrad_;
};

} // end namespace

#endif
