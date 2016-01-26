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
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc fully implicit model.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitFluxVariables
 * \brief Variables for the enthalpy fluxes in the MpNc model
 */
template <class TypeTag, bool enableEnergy/*=false*/, int numEnergyEquations/*=0*/>
class MPNCFluxVariablesEnergy
{
    static_assert(!(numEnergyEquations && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!numEnergyEquations,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     */
    MPNCFluxVariablesEnergy()
    {
    }
    /*!
     * \brief update
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables of the current element
     */
    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const FluxVariables & fluxVars,
                const ElementVolumeVariables & elemVolVars)
    {}
};

template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true,  /*numEnergyEquations=*/1>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum{dimWorld = GridView::dimensionworld};
    enum{dim = GridView::dimension};
    enum{nPhaseIdx = FluidSystem::nPhaseIdx};
    enum{wPhaseIdx = FluidSystem::wPhaseIdx};

    typedef Dune::FieldVector<CoordScalar, dim>  DimVector;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;

public:
    /*!
     * \brief The constructor
     */
    MPNCFluxVariablesEnergy()
    {}
    /*!
     * \brief update
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables of the current element
     */
    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const FluxVariables & fluxVars,
                const ElementVolumeVariables & elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        GlobalPosition tmp(0.0);
        GlobalPosition temperatureGradient(0.);
        for (int idx = 0; idx < face.numFap; idx++)
        {
            tmp = face.grad[idx];

            // index for the element volume variables
            int volVarsIdx = face.fapIndices[idx];

            tmp *= elemVolVars[volVarsIdx].fluidState().temperature(/*phaseIdx=*/0);
            temperatureGradient += tmp;
        }

        // project the heat flux vector on the face's normal vector
        temperatureGradientNormal_ = temperatureGradient * face.normal;

        lambdaEff_ = 0;
        calculateEffThermalConductivity_(problem,
                                         element,
                                         fvGeometry,
                                         face,
                                         elemVolVars);
    }

    /*!
     * \brief The lumped / average conductivity of solid plus phases \f$[W/mK]\f$.
     */
    Scalar lambdaEff() const
    { return lambdaEff_; }

    /*!
     * \brief The normal of the gradient of temperature .
     */
    Scalar temperatureGradientNormal() const
    {
        return temperatureGradientNormal_;
    }

protected:
    /*!
     * \brief Calculate the effective thermal conductivity of
     *        the porous medium plus residing phases \f$[W/mK]\f$.
     *        This basically means to access the model for averaging
     *        the individual conductivities, set by the property ThermalConductivityModel.
     *        Except the adapted arguments, this is the same function
     *        as used in the implicit TwoPTwoCNIFluxVariables.
     */
    void calculateEffThermalConductivity_(const Problem &problem,
                                          const Element &element,
                                          const FVElementGeometry & fvGeometry,
                                          const SCVFace & face,
                                          const ElementVolumeVariables &elemVolVars)
    {
        const unsigned i = face.i;
        const unsigned j = face.j;
        Scalar lambdaI, lambdaJ;

        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            lambdaI =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i].saturation(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(element, fvGeometry, i),
                                                                   problem.spatialParams().porosity(element, fvGeometry, i));
            lambdaJ =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j].saturation(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(element, fvGeometry, j),
                                                                   problem.spatialParams().porosity(element, fvGeometry, j));
        }
        else
        {
            const Element & elementI = fvGeometry.neighbors[i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            lambdaI =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i].saturation(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(elementI, fvGeometryI, 0),
                                                                   problem.spatialParams().porosity(elementI, fvGeometryI, 0));

            const Element & elementJ = fvGeometry.neighbors[j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            lambdaJ =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j].saturation(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(elementJ, fvGeometryJ, 0),
                                                                   problem.spatialParams().porosity(elementJ, fvGeometryJ, 0));
        }

        // -> harmonic mean
        lambdaEff_ = harmonicMean(lambdaI, lambdaJ);
    }

private:
    Scalar lambdaEff_ ;
    Scalar temperatureGradientNormal_;
};

} // end namespace

#endif
