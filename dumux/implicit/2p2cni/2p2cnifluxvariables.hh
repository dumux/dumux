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
 * \brief This file contains data which is required to calculate
 *        the heat fluxes over a face of a finite volume.
 *
 * This means temperature gradients and the normal matrix
 * heat flux.
 */
#ifndef DUMUX_2P2CNI_FLUX_VARIABLES_HH
#define DUMUX_2P2CNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/2p2c/2p2cfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the heat fluxes over a face of a finite
 *        volume for the non-isothermal two-phase, two-component model.
 *        The mass fluxes are computed in the parent class.
 *
 * This means temperature gradients and the normal matrix
 * heat flux.
 */
template <class TypeTag>
class TwoPTwoCNIFluxVariables : public TwoPTwoCFluxVariables<TypeTag>
{
    typedef TwoPTwoCFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;

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

    /*!
     * \brief The local temperature gradient at the IP of the considered scv face.
     */
    DimVector temperatureGradient() const
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
        DimVector tmp(0.0);
        for (int idx = 0; idx < this->fvGeometry_.numFap; idx++)
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
        normalMatrixHeatFlux_ = temperatureGrad_*
                                this->face().normal;
        normalMatrixHeatFlux_ *= -lambdaEff_;
    }

    void calculateEffThermalConductivity_(const Problem &problem,
                                        const Element &element,
                                        const ElementVolumeVariables &elemVolVars)
    {
        const Scalar lambdaI =
                ThermalConductivityModel::effectiveThermalConductivity(element,
                                                                       elemVolVars,
                                                                       this->fvGeometry_,
                                                                       problem.spatialParams(),
                                                                       this->face().i);
        const Scalar lambdaJ =
                ThermalConductivityModel::effectiveThermalConductivity(element,
                                                                       elemVolVars,
                                                                       this->fvGeometry_,
                                                                       problem.spatialParams(),
                                                                       this->face().j);
        // -> harmonic mean
        lambdaEff_ = harmonicMean(lambdaI, lambdaJ);
    }

private:
    Scalar lambdaEff_;
    Scalar normalMatrixHeatFlux_;
    DimVector temperatureGrad_;
    int faceIdx_;
};

} // end namespace

#endif
