// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
* \ingroup CCTpfaFlux
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
*/
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_NONEQUILIBRIUM_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_NONEQUILIBRIUM_HH

#include <dumux/common/deprecated.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FouriersLawNonEquilibriumImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag>
class FouriersLawNonEquilibriumImplementation<TypeTag, DiscretizationMethod::cctpfa>
{
    using Implementation = FouriersLawNonEquilibriumImplementation<TypeTag, DiscretizationMethod::cctpfa>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVarsCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto numEnergyEqFluid = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr auto sPhaseIdx = ModelTraits::numFluidPhases();

public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    using Cache = FluxVariablesCaching::EmptyHeatConductionCache;

    //! Compute the heat condution flux assuming thermal equilibrium
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        Scalar tInside = 0.0;
        Scalar tOutside = 0.0;
        // get the inside/outside temperatures
        if (phaseIdx < numEnergyEqFluid)
        {
            tInside += elemVolVars[scvf.insideScvIdx()].temperatureFluid(phaseIdx);
            tOutside += elemVolVars[scvf.outsideScvIdx()].temperatureFluid(phaseIdx);
        }
        else //temp solid
        {
            tInside += elemVolVars[scvf.insideScvIdx()].temperatureSolid();
            tOutside += elemVolVars[scvf.outsideScvIdx()].temperatureSolid();
        }

        Scalar tij = calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx);
        return tij*(tInside - tOutside);
    }

    //! Compute transmissibilities
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf,
                                            const int phaseIdx)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        Scalar insideLambda = 0.0;
        Scalar outsideLambda = 0.0;

        // effective diffusion tensors
        if (phaseIdx != sPhaseIdx)
        {
            //when number of energyEq for the fluid are smaller than numPhases that means that we need an effecitve law
           if (numEnergyEqFluid < ModelTraits::numFluidPhases())
            {
                insideLambda += insideVolVars.effectiveThermalConductivity();
            }
            else //numEnergyEqFluid >1
            {
                insideLambda += insideVolVars.fluidThermalConductivity(phaseIdx)*insideVolVars.saturation(phaseIdx)*insideVolVars.porosity();
            }
        }
        //solid phase
        else
        {
            insideLambda += insideVolVars.solidThermalConductivity()*(1.0-insideVolVars.porosity());
        }

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, insideLambda, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = scvf.area()*ti;
        }
        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

       // effective diffusion tensors
        if (phaseIdx != sPhaseIdx)
        {
            //when number of energyEq for the fluid are smaller than numPhases that means that we need an effecitve law
           if (numEnergyEqFluid < ModelTraits::numFluidPhases())
            {
                outsideLambda += outsideVolVars.effectiveThermalConductivity();
            }
            else
            {
                outsideLambda += outsideVolVars.fluidThermalConductivity(phaseIdx)*outsideVolVars.saturation(phaseIdx)*outsideVolVars.porosity();
            }
        }
        //solid phase
        else
        {
            outsideLambda +=outsideVolVars.solidThermalConductivity()*(1.0-outsideVolVars.porosity());
        }
            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*computeTpfaTransmissibility(scvf, outsideScv, outsideLambda, outsideVolVars.extrusionFactor());
            else
                tj = computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv, outsideLambda, outsideVolVars.extrusionFactor());

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }
        return tij;
    }
};

} // end namespace Dumux

#endif
