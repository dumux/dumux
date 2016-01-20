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
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc box model with kinetic energy transfer enabled.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_KINETIC_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_KINETIC_HH

#include <dune/common/fvector.hh>

#include <dumux/common/spline.hh>
#include <dumux/porousmediumflow/mpnc/implicit/fluxvariables.hh>

namespace Dumux
{

/*!
 * \brief Specialization for the case of *3* energy balance equations.
 */
template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/3>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices)  Indices;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {dim = GridView::dimension};
    enum {dimWorld = GridView::dimensionworld};
    enum {numEnergyEqs             = Indices::numPrimaryEnergyVars};

    typedef Dune::FieldVector<CoordScalar, dim>  DimVector;
    typedef Dune::FieldVector<CoordScalar, dimWorld>  GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

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
        GlobalPosition tmp ;

        for(int energyEqIdx=0; energyEqIdx<numEnergyEqs; energyEqIdx++)
            temperatureGradient_[energyEqIdx] = 0.;

        for (unsigned int idx = 0;
                idx < face.numFap;
                idx++){
            // FE gradient at vertex idx
            const GlobalPosition & feGrad = face.grad[idx];

            for (int energyEqIdx =0; energyEqIdx < numEnergyEqs; ++energyEqIdx){
                // index for the element volume variables
                int volVarsIdx = face.fapIndices[idx];

                tmp = feGrad;
                tmp   *= elemVolVars[volVarsIdx].temperature(energyEqIdx);
                temperatureGradient_[energyEqIdx] += tmp;
            }
        }
    }

    /*!
     * \brief The total heat flux \f$[J/s]\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face.
     *
     * \param energyEqIdx The index of the energy equation
     */
    GlobalPosition temperatureGradient(const unsigned int energyEqIdx) const
    {
        return temperatureGradient_[energyEqIdx];
    }

private:
    GlobalPosition temperatureGradient_[numEnergyEqs];
};


/*!
 * \brief Specialization for the case of *2* energy balance equations.
 */
template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/2>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices)  Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)  Scalar;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {dim = GridView::dimension};
    enum {dimWorld = GridView::dimensionworld};
    enum {numEnergyEqs          = Indices::numPrimaryEnergyVars};
    enum {wPhaseIdx             = FluidSystem::wPhaseIdx};
    enum {nPhaseIdx             = FluidSystem::nPhaseIdx};


    typedef Dune::FieldVector<CoordScalar, dimWorld>  GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
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
        GlobalPosition tmp ;

        for(int energyEqIdx=0; energyEqIdx<numEnergyEqs; energyEqIdx++)
            temperatureGradient_[energyEqIdx] = 0.;

        for (unsigned int idx = 0;
                idx < face.numFap;
                idx++){
            // FE gradient at vertex idx
            const GlobalPosition & feGrad = face.grad[idx];

            for (int energyEqIdx =0; energyEqIdx < numEnergyEqs; ++energyEqIdx){
                // index for the element volume variables
                int volVarsIdx = face.fapIndices[idx];

                tmp = feGrad;
                tmp   *= elemVolVars[volVarsIdx].temperature(energyEqIdx);
                temperatureGradient_[energyEqIdx] += tmp;
            }
        }

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
     * \brief The total heat flux \f$[J/s]\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face.
     *
     * \param energyEqIdx The index of the energy equation
     */
    GlobalPosition temperatureGradient(const unsigned int energyEqIdx) const
    {
        return temperatureGradient_[energyEqIdx];
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
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i].fluidState().saturation(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(element, fvGeometry, i),
                                                                   problem.spatialParams().porosity(element, fvGeometry, i));
            lambdaJ =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j].fluidState().saturation(wPhaseIdx),
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
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[i].fluidState().saturation(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[i].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(elementI, fvGeometryI, 0),
                                                                   problem.spatialParams().porosity(elementI, fvGeometryI, 0));

            const Element & elementJ = fvGeometry.neighbors[j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            lambdaJ =
                ThermalConductivityModel::effectiveThermalConductivity(elemVolVars[j].fluidState().saturation(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(wPhaseIdx),
                                                                   elemVolVars[j].thermalConductivity(nPhaseIdx),
                                                                   problem.spatialParams().solidThermalConductivity(elementJ, fvGeometryJ, 0),
                                                                   problem.spatialParams().porosity(elementJ, fvGeometryJ, 0));
        }

        // -> arithmetic mean, open to discussion
        lambdaEff_ = 0.5 * (lambdaI+lambdaJ);
    }


private:
    Scalar lambdaEff_ ;
    GlobalPosition temperatureGradient_[numEnergyEqs];
};

} // end namespace

#endif
