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
#ifndef DUMUX_NETWORK_TISSUE_TRANSPORT_PROPERTIES_HH
#define DUMUX_NETWORK_TISSUE_TRANSPORT_PROPERTIES_HH

// # Model (`TypeTag`) and property specializations (`properties.hh`)
//
// This file contains the __property definitions__ which configure
// the model and discretization method (at compile time). This file
// implements three models: a common base model `NetworkTissueModel`,
// and the sub-domain models `NetworkTransportModel` and `TissueTransportModel`
// that inherit some of the basic properties from the common ancestor `NetworkTissueModel`.
// By specializing properties for a model tag, every class that knows
// the model tag can extract these properties at compile time. (In C++ speak, the
// template structs that are specialized for different types are called "traits".)
//
// [[content]]
//
// ### Include headers
// [[codeblock]]
// 3d Cartesian grid implementation
#include <dune/grid/yaspgrid.hh>
// 1d in 3d embedded network grid implementation
#include <dune/foamgrid/foamgrid.hh>
// properties
#include <dumux/common/properties.hh>
// dof mapper that reorders dofs for optimal matrix pattern
#include <dumux/common/reorderingdofmapper.hh>
// discretization method
#include <dumux/discretization/cctpfa.hh>
// model equations
#include <dumux/porousmediumflow/tracer/model.hh>
// effective diffusion coefficient model
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
// 1d-3d coupling manager with surface average operator
#include <dumux/multidomain/embedded/couplingmanager1d3d_average.hh>

#include "spatialparams.hh"
#include "problem.hh"
#include "tracerfluidsystem.hh"
// [[/codeblock]]
//
// ## Common properties of both models
// [[codeblock]]
namespace Dumux::Properties {

// inheriting all properties from the tracer model (advection-diffusion equation)
// and inheriting all properties from the cell-centered finite volume discretization with two-point flux approximation
namespace TTag { struct NetworkTissueModel { using InheritsFrom = std::tuple<Tracer, CCTpfaModel>; }; }

// caching options, enabling caching uses more memory but is faster in terms of runtime
template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = true; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::NetworkTissueModel> { static constexpr bool value = false; };

// Set the fluid system (see `tracerfluidsystem.hh`)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::NetworkTissueModel>
{ using type = TracerFluidSystem<double>; };

// Use molar balances and a mole fraction as primary variable
template<class TypeTag>
struct UseMoles<TypeTag, TTag::NetworkTissueModel>
{ static constexpr bool value = true; };

} // end namespace Dumux::Properties
// [[/codeblock]]
//
// ## Network transport model
// [[codeblock]]
namespace Dumux::Properties {

namespace TTag { struct NetworkTransportModel { using InheritsFrom = std::tuple<NetworkTissueModel>; }; }

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::NetworkTransportModel>
{ using type = Dune::FoamGrid<1, 3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::NetworkTransportModel>
{ using type = NetworkTransportProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::NetworkTransportModel>
{ using type = NetworkSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (cc)
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::NetworkTransportModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using ElementMapper = ReorderingDofMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

template<class Scalar>
class EffectiveDiffusionTissuePVS
{
public:
    template<class VolumeVariables>
    static Scalar effectiveDiffusionCoefficient(const VolumeVariables& volVars,
                                                const int phaseIdx,
                                                const int compIdxI,
                                                const int compIdxJ)
    {
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        return volVars.porosity() * diffCoeff;
    }
};

template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::NetworkTransportModel>
{ using type = EffectiveDiffusionTissuePVS<double>; };

} // end namespace Dumux::Properties
// [[/codeblock]]
//
// ## Tissue transport model
// [[codeblock]]
namespace Dumux::Properties {

namespace TTag { struct TissueTransportModel { using InheritsFrom = std::tuple<NetworkTissueModel>; }; }

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TissueTransportModel>
{ using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TissueTransportModel>
{ using type = TissueTransportProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TissueTransportModel>
{ using type = TissueSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

template<class Scalar>
class EffectiveDiffusionTissue
{
public:
    template<class VolumeVariables>
    static Scalar effectiveDiffusionCoefficient(const VolumeVariables& volVars,
                                                const int phaseIdx,
                                                const int compIdxI,
                                                const int compIdxJ)
    {
        static const Scalar tau = getParam<Scalar>("Tissue.SpatialParams.Tortuosity");
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        return tau * volVars.porosity() * diffCoeff;
    }
};

// effective diffusivity model
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::TissueTransportModel>
{ using type = EffectiveDiffusionTissue<double>; };

} // end namespace Dumux::Properties
// [[/codeblock]]
//
// ## Multi-domain coupling
// [[codeblock]]
namespace Dumux::Properties {

using CouplingTransport = Embedded1d3dCouplingManager<MultiDomainTraits<
    Properties::TTag::TissueTransportModel, Properties::TTag::NetworkTransportModel>,
    Embedded1d3dCouplingMode::Average
>;

// tell the tissue sub-model about the coupling
template<class TypeTag> struct CouplingManager<TypeTag, TTag::TissueTransportModel> { using type = CouplingTransport; };
template<class TypeTag> struct PointSource<TypeTag, TTag::TissueTransportModel> { using type = CouplingTransport::PointSourceTraits::template PointSource<0>; };
template<class TypeTag> struct PointSourceHelper<TypeTag, TTag::TissueTransportModel> { using type = CouplingTransport::PointSourceTraits::template PointSourceHelper<0>; };

// tell the network sub-model about the coupling
template<class TypeTag> struct CouplingManager<TypeTag, TTag::NetworkTransportModel> { using type = CouplingTransport; };
template<class TypeTag> struct PointSource<TypeTag, TTag::NetworkTransportModel> { using type = CouplingTransport::PointSourceTraits::template PointSource<1>; };
template<class TypeTag> struct PointSourceHelper<TypeTag, TTag::NetworkTransportModel> { using type = CouplingTransport::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties
// [[/codeblock]]
// [[/content]]
#endif
