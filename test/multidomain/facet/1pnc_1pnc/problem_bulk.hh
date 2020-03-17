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
 * \ingroup FacetTests
 * \brief The problem for the bulk domain in the 1pnc facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_BULKPROBLEM_HH

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePNCBulkProblem;

namespace Properties {

// create the type tag nodes
namespace TTag {
struct BaseBulk {};
struct OnePNCBulk { using InheritsFrom = std::tuple<OnePNC, BaseBulk>; };
struct OnePNCNIBulk { using InheritsFrom = std::tuple<OnePNCNI, BaseBulk>; };

struct OnePNCBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCNIBulk>; };
} // end namespace TTag

// Set the grid type (DIMWORLD is defined in CMakeLists.txt)
template<class TypeTag>
struct Grid<TypeTag, TTag::BaseBulk> { using type = Dune::ALUGrid<2, DIMWORLD, Dune::cube, Dune::nonconforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::BaseBulk> { using type = OnePNCBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BaseBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BaseBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
public:
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief Test problem for the 1pnc model with
 *        coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePNCBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr bool enableHeatConduction = ModelTraits::enableEnergyBalance();

    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx)
    };

public:
    OnePNCBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManager,
                      const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if (!getPropValue<TypeTag, Properties::UseMoles>())
            DUNE_THROW(Dune::InvalidStateException, "Problem is implemented for molar formulation!");

        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[1] < 1e-6 || globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
            values.setAllDirichlet();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);

        if (globalPos[1] < 1e-6)
        {
            values[N2Idx] = 1e-3;
            if constexpr (enableHeatConduction)
            {
                values[Indices::temperatureIdx] += 100; // DeltaT = 100 K
            }
        }

        return values;
    }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = 1.0e5;
        values[N2Idx] = 0.0;

        if constexpr (enableHeatConduction)
        {
            values[Indices::temperatureIdx] = temperature();
        }

        return values;
    }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
