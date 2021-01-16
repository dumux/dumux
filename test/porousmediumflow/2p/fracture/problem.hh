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
 * \ingroup TwoPTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *
 * The fracture is a 2D network embedded in 3D.
 */
#ifndef DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class FractureProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Fracture { using InheritsFrom = std::tuple<TwoP>; };
struct FractureBox { using InheritsFrom = std::tuple<Fracture, BoxModel>; };
struct FractureCCTpfa { using InheritsFrom = std::tuple<Fracture, CCTpfaModel>; };
struct FractureCCMpfa { using InheritsFrom = std::tuple<Fracture, CCMpfaModel>; };
} // end namespace TTag

// set the grid property
#if HAVE_DUNE_FOAMGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::Fracture> { using type = Dune::FoamGrid<2, 3>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Fracture> { using type = Dumux::FractureProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Fracture>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Fracture>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FractureSpatialParams<GridGeometry, Scalar>;
};

// Use global caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };

// permeablility is solution-independent
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Fracture> { static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief Trichloroethene (DNAPL) transport through a fracture network (2d in 3d).
 */
template <class TypeTag>
class FractureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        // equation indices
        contiTCEEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FractureProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position where to set the BC types
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    {
        return 0.1;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position where to set the BC types
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();
        if (onInlet_(globalPos))
            values.setAllNeumann();
        if (globalPos[2] > 1.0 - eps_ || globalPos[2] < eps_)
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        const auto depth = this->gridGeometry().bBoxMax()[dimWorld-1] - globalPos[dimWorld-1];
        const auto g = this->spatialParams().gravity(globalPos)[dimWorld-1];

        PrimaryVariables values;
        values[pressureIdx] = 1e5 + 1000*g*depth;
        values[saturationIdx] = 0.0;
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos)) {
            values[contiTCEEqIdx] = -0.04; // kg / (m * s)
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos); }
    // \}

private:

    bool onInlet_(const GlobalPosition &globalPos) const
    { return globalPos[0] < eps_ && globalPos[1] > -0.5 - eps_; }

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} // end namespace Dumux

#endif
