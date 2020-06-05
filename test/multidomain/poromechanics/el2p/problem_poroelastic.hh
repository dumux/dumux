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
 * \ingroup PoromechanicsTests
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
#ifndef DUMUX_POROELASTIC_SUBPROBLEM_HH
#define DUMUX_POROELASTIC_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>

#include "spatialparams_poroelastic.hh"
#include "co2tables_el2p.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class PoroElasticSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct PoroElasticSub { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = Dune::YaspGrid<3>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticSub> { using type = Dumux::PoroElasticSubProblem<TypeTag>; };

// Set the fluid system for TwoPSubProblem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineCO2<Scalar, El2P::CO2Tables>;
};

// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticSub>
{
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry> >;
};
} // end namespace Properties

/*!
 * \ingroup PoromechanicsTests
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
template<class TypeTag>
class PoroElasticSubProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    PoroElasticSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<CouplingManager> couplingManagerPtr,
                          const std::string& paramGroup = "PoroElastic")
    : ParentType(gridGeometry, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Returns the temperature in the domain
    static constexpr Scalar temperature()
    { return 273.15; }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Neumann boundary segment.
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the effective fluid density.
     */
    Scalar effectiveFluidDensity(const Element& element, const SubControlVolume& scv) const
    {
        // get porous medium flow volume variables from coupling manager
        const auto pmFlowVolVars = couplingManager().getPMFlowVolVars(element);

        Scalar wPhaseDensity = pmFlowVolVars.density(FluidSystem::phase0Idx);
        Scalar nPhaseDensity = pmFlowVolVars.density(FluidSystem::phase1Idx);
        Scalar Sw = pmFlowVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = pmFlowVolVars.saturation(FluidSystem::phase1Idx);
        return (wPhaseDensity * Sw + nPhaseDensity * Sn);
    }

    /*!
     * \brief Returns the effective pore pressure.
     */
    template< class FluxVarsCache >
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    {
        // get porous medium flow volume variables from coupling manager
        const auto pmFlowVolVars = couplingManager().getPMFlowVolVars(element);

        Scalar pw = pmFlowVolVars.pressure(FluidSystem::phase0Idx);
        Scalar pn = pmFlowVolVars.pressure(FluidSystem::phase1Idx);
        Scalar Sw = pmFlowVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = pmFlowVolVars.saturation(FluidSystem::phase1Idx);
        return (pw * Sw + pn * Sn);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
