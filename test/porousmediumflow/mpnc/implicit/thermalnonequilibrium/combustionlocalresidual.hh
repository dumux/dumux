// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup MPNCTests
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium model.
 */
#ifndef DUMUX_ENERGY_COMBUSTION_LOCAL_RESIDUAL_HH
#define DUMUX_ENERGY_COMBUSTION_LOCAL_RESIDUAL_HH

#include <cmath>
#include <dune/common/math.hh>
#include <dumux/common/spline.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup MPNCTests
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium  model.
 */
template<class TypeTag>
class CombustionEnergyLocalResidual
: public  EnergyLocalResidualNonEquilibrium<TypeTag, 1/*numEnergyEqFluid*/>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr auto numEnergyEqFluid = ModelTraits::numEnergyEqFluid();
    static constexpr auto numEnergyEqSolid = ModelTraits::numEnergyEqSolid();
    static constexpr auto energyEq0Idx = Indices::energyEq0Idx;

public:
    /*!
     * \brief Calculates the source term of the equation
     *
     * \param source The source term
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {
        //specialization for 2 fluid phases
        const auto& volVars = elemVolVars[scv];
        const auto& fs = volVars.fluidState() ;
        const Scalar characteristicLength = volVars.characteristicLength()  ;

        //interfacial area
        // Shi & Wang, Transport in porous media (2011)
        const Scalar as = volVars.fluidSolidInterfacialArea();

        //temperature fluid is the same for both fluids
        const Scalar TFluid = volVars.temperatureFluid(0);
        const Scalar TSolid = volVars.temperatureSolid();

        const Scalar satW = fs.saturation(0) ;
        const Scalar satN = fs.saturation(1) ;

        const Scalar eps = 1e-6 ;
        Scalar solidToFluidEnergyExchange ;

        Scalar fluidConductivity ;
        if (satW < 1.0 - eps)
            fluidConductivity = volVars.fluidThermalConductivity(1) ;
        else if (satW >= 1.0 - eps)
            fluidConductivity = volVars.fluidThermalConductivity(0) ;
        else
            DUNE_THROW(Dune::NotImplemented, "wrong range");

        const Scalar factorEnergyTransfer   = volVars.factorEnergyTransfer()  ;

        solidToFluidEnergyExchange = factorEnergyTransfer * (TSolid - TFluid) / characteristicLength * as * fluidConductivity ;
        const Scalar epsRegul = 1e-3 ;

        if (satW < (0 + eps) )
        {
            solidToFluidEnergyExchange *=  volVars.nusseltNumber(1) ;
        }
        else if ( (satW >= 0 + eps) and (satW < 1.0-eps) )
        {
            solidToFluidEnergyExchange *=  (volVars.nusseltNumber(1) * satN );
            Scalar qBoil ;
            if (satW<=epsRegul)
            {// regularize
                typedef Dumux::Spline<Scalar> Spline;
                    Spline sp(0.0, epsRegul,                           // x1, x2
                        QBoilFunc(volVars, 0.0), QBoilFunc(volVars, epsRegul),       // y1, y2
                        0.0,dQBoil_dSw(volVars, epsRegul));    // m1, m2

                qBoil = sp.eval(satW);
            }

            else if (satW>= (1.0-epsRegul) )
            {   // regularize
                typedef Dumux::Spline<Scalar> Spline;
                Spline sp(1.0-epsRegul, 1.0,    // x1, x2
                        QBoilFunc(volVars, 1.0-epsRegul), 0.0,    // y1, y2
                        dQBoil_dSw(volVars, 1.0-epsRegul), 0.0 );      // m1, m2

                qBoil = sp.eval(satW) ;
            }
            else
            {
                qBoil = QBoilFunc(volVars, satW)  ;
            }

            solidToFluidEnergyExchange += qBoil;
        }
        else if (satW >= 1.0-eps)
        {
            solidToFluidEnergyExchange *=  volVars.nusseltNumber(0) ;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "wrong range");

        using std::isfinite;
        if (!isfinite(solidToFluidEnergyExchange))
            DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "TFluid="<< TFluid << " TSolid="<< TSolid);

        for(int energyEqIdx =0; energyEqIdx<numEnergyEqFluid+numEnergyEqSolid; ++energyEqIdx)
        {
            switch (energyEqIdx)
            {
            case 0 :
                source[energyEq0Idx + energyEqIdx] += solidToFluidEnergyExchange;
                break;
            case 1 :
                source[energyEq0Idx + energyEqIdx] -=   solidToFluidEnergyExchange;
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                        "wrong index");
            } // end switch
        }// end energyEqIdx
    }// end source


  /*! \brief Calculates the energy transfer during boiling, i.e. latent heat
   *
   * \param volVars The volume variables
   * \param satW The wetting phase saturation. Not taken from volVars, because we regularize.
   */
    static Scalar QBoilFunc(const VolumeVariables & volVars,
                            const  Scalar satW)
    {
        // using saturation as input (instead of from volVars)
        // in order to make regularization (evaluation at different points) easyer
        const auto& fs = volVars.fluidState();
        const Scalar g( 9.81 );
        const Scalar gamma(0.0589);
        const Scalar TSolid = volVars.temperatureSolid();
        using std::pow;
        using Dune::power;
        const Scalar as = volVars.fluidSolidInterfacialArea();
        const Scalar mul = fs.viscosity(0);
        const Scalar deltahv = fs.enthalpy(1) - fs.enthalpy(0);
        const Scalar deltaRho = fs.density(0) - fs.density(1);
        const Scalar firstBracket = pow(g * deltaRho / gamma, 0.5);
        const Scalar cp = FluidSystem::heatCapacity(fs, 0);
        // This use of Tsat is only justified if the fluid is always boiling (tsat equals boiling conditions)
        // If a different state is to be simulated, please use the actual fluid temperature instead.
        const Scalar Tsat = FluidSystem::vaporTemperature(fs, 1 ) ;
        const Scalar deltaT = TSolid - Tsat;
        const Scalar secondBracket = power( (cp *deltaT / (0.006 * deltahv)  ) , 3);
        const Scalar Prl = volVars.prandtlNumber(0);
        const Scalar thirdBracket = pow( 1/Prl , (1.7/0.33));
        const Scalar QBoil = satW * as * mul * deltahv * firstBracket * secondBracket * thirdBracket;
        return QBoil;
    }

    /*! \brief Calculates the derivative of the energy transfer function during boiling. Needed for regularization.
   *
   * \param volVars The volume variables
   * \param satW The wetting phase saturation. Not taken from volVars, because we regularize.
   */
    static Scalar dQBoil_dSw(const VolumeVariables & volVars,
                                const Scalar satW)
    {
        // on the fly derive w.r.t. Sw.
        // Only linearly depending on it (directly)
        return (QBoilFunc(volVars, satW) / satW ) ;
    }
};
} // end namespace Dumux

#endif
