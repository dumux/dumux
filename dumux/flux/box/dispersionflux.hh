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
 *        dispersive fluxes.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_DISPERSION_FLUX_HH
#define DUMUX_DISCRETIZATION_BOX_DISPERSION_FLUX_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/traits.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of a dispersion flux for the box method
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation<TypeTag, DiscretizationMethods::Box, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVarCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluxTraits = typename Dumux::FluxTraits<FluxVariables>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = ModelTraits::numFluidPhases(),
        numComponents = ModelTraits::numFluidComponents()
    };

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using HeatFluxScalar = Scalar;

    static constexpr bool stationaryVelocityField = FluxTraits::hasStationaryVelocityField();

public:

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    /*!
     * \brief Returns the dispersive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
    static ComponentFluxVector compositionalDispersionFlux(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace& scvf,
                                                           const int phaseIdx,
                                                           const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        ComponentFluxVector componentFlux(0.0);

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

        // density interpolation
        Scalar rhoMassOrMole(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto rho = massOrMolarDensity(elemVolVars[scv], referenceSystem, phaseIdx);
            rhoMassOrMole += rho * shapeValues[scv.indexInElement()][0];
        }

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            // collect the dispersion tensor, the fluxVarsCache and the shape values
            const auto& dispersionTensor =
                VolumeVariables::CompositionalDispersionModel::compositionalDispersionTensor(problem, scvf, fvGeometry,
                                                                                             elemVolVars, elemFluxVarsCache,
                                                                                             phaseIdx, compIdx);

            // the mole/mass fraction gradient
            Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto x = massOrMoleFraction(elemVolVars[scv], referenceSystem, phaseIdx, compIdx);
                gradX.axpy(x, fluxVarsCache.gradN(scv.indexInElement()));
            }

            // compute the dispersion flux
            componentFlux[compIdx] = -1.0 * rhoMassOrMole * vtmv(scvf.unitOuterNormal(), dispersionTensor, gradX) * scvf.area();
            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[phaseIdx] -= componentFlux[compIdx];
        }
        return componentFlux;
    }

    /*!
     * \brief Returns the thermal dispersive flux
     *        across the given sub-control volume face.
     */
    static HeatFluxScalar thermalDispersionFlux(const Problem& problem,
                                                const Element& element,
                                                const FVElementGeometry& fvGeometry,
                                                const ElementVolumeVariables& elemVolVars,
                                                const SubControlVolumeFace& scvf,
                                                const int phaseIdx,
                                                const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // collect the dispersion tensor
        const auto& dispersionTensor =
            ModelTraits::ThermalDispersionModel::thermalDispersionTensor(problem, scvf, fvGeometry,
                                                                         elemVolVars, elemFluxVarsCache,
                                                                         phaseIdx);
        // compute the temperature gradient with the shape functions
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, GridView::dimensionworld> gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
            gradTemp.axpy(elemVolVars[scv].temperature(), fluxVarsCache.gradN(scv.indexInElement()));

        // compute the heat conduction flux
        return -1.0*vtmv(scvf.unitOuterNormal(), dispersionTensor, gradTemp)*Extrusion::area(scvf);
    }

};

} // end namespace Dumux

#endif
