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
 * \brief This file contains the data which is required to calculate
 *        energy fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH
#define DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>

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
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using IndexType = typename GridView::IndexSet::IndexType;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;

    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numPhases()} ;
    enum { numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid) };
    enum {sPhaseIdx = FluidSystem::numPhases};

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

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
           if (numEnergyEqFluid == 1)
            {   //when only one energy equation for fluids is used, we need an effective law for that
                insideLambda += ThermalConductivityModel::effectiveThermalConductivity(insideVolVars, problem.spatialParams(), element, fvGeometry, insideScv);
                outsideLambda += ThermalConductivityModel::effectiveThermalConductivity(outsideVolVars, problem.spatialParams(), element, fvGeometry, outsideScv);
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
            insideLambda += insideVolVars.solidThermalConductivity()*(1-insideVolVars.porosity());
            outsideLambda +=outsideVolVars.solidThermalConductivity()*(1-outsideVolVars.porosity());
        }

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged diffusion tensor
        const auto lambda = problem.spatialParams().harmonicMean(insideLambda, outsideLambda, scvf.unitOuterNormal());

        // evaluate gradTemp at integration point
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        Dune::FieldVector<Scalar, dimWorld> gradTemp(0.0);
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
