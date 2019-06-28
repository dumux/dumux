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
 * \ingroup BoxFlux
 * \brief This file contains the data which is required to calculate
 *        energy fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH
#define DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/deprecated.hh>

#include <dumux/discretization/method.hh>


namespace Dumux {

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod>
class FouriersLawNonEquilibriumImplementation;

/*!
 * \ingroup BoxFouriersLaw
 * \brief Specialization of Fourier's Law for the box method for thermal nonequilibrium models.
 */
template <class TypeTag>
class FouriersLawNonEquilibriumImplementation<TypeTag, DiscretizationMethod::box>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto numEnergyEqFluid = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr auto sPhaseIdx = ModelTraits::numFluidPhases();

public:
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        Scalar insideLambda = 0.0;
        Scalar outsideLambda = 0.0;
       // effective diffusion tensors
        if (phaseIdx != sPhaseIdx)
        {
            //when number of energyEq for the fluid are smaller than numPhases that means that we need an effecitve law
            if (numEnergyEqFluid < ModelTraits::numFluidPhases())
            {
                insideLambda += insideVolVars.effectiveThermalConductivity();
                outsideLambda += outsideVolVars.effectiveThermalConductivity();
            }
            else
            {
                insideLambda += insideVolVars.fluidThermalConductivity(phaseIdx)*insideVolVars.saturation(phaseIdx)*insideVolVars.porosity();
                outsideLambda += outsideVolVars.fluidThermalConductivity(phaseIdx)*outsideVolVars.saturation(phaseIdx)*outsideVolVars.porosity();
            }
        }
        //solid phase
        else
        {
            insideLambda += insideVolVars.solidThermalConductivity()*(1.0-insideVolVars.porosity());
            outsideLambda += outsideVolVars.solidThermalConductivity()*(1.0-outsideVolVars.porosity());
        }

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged diffusion tensor
        const auto lambda = problem.spatialParams().harmonicMean(insideLambda, outsideLambda, scvf.unitOuterNormal());

        // evaluate gradTemp at integration point
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        Dune::FieldVector<Scalar, GridView::dimensionworld> gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            // compute the temperature gradient with the shape functions
            if (phaseIdx < numEnergyEqFluid)
                gradTemp.axpy(elemVolVars[scv].temperatureFluid(phaseIdx), fluxVarsCache.gradN(scv.indexInElement()));
            else
               gradTemp.axpy(elemVolVars[scv].temperatureSolid(), fluxVarsCache.gradN(scv.indexInElement()));
        }

        // comute the heat conduction flux
        return -1.0*vtmv(scvf.unitOuterNormal(), lambda, gradTemp)*scvf.area();
    }
};

} // end namespace Dumux

#endif
