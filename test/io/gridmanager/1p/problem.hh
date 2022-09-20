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
 * \ingroup OnePTests
 * \brief
 */

#ifndef DUMUX_LANDFILLGRIDMANAGER_PROBLEM_HH
#define DUMUX_LANDFILLGRIDMANAGER_PROBLEM_HH

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/common/timeloop.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class LandfillGridmanagerTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct LandfillGridmanagerTest { using InheritsFrom = std::tuple<OneP>; };
struct LandfillGridmanagerTestBox { using InheritsFrom = std::tuple<LandfillGridmanagerTest, BoxModel>; };
struct LandfillGridmanagerTestCCTpfa { using InheritsFrom = std::tuple<LandfillGridmanagerTest, CCTpfaModel>; };
struct LandfillGridmanagerTestCCMpfa { using InheritsFrom = std::tuple<LandfillGridmanagerTest, CCMpfaModel>; };
} // end namespace TTag

// Specialize the fluid system type for this type tag
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LandfillGridmanagerTest>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Specialize the grid type for this type tag
template<class TypeTag>
struct Grid<TypeTag, TTag::LandfillGridmanagerTest>
{ using type = Dune::UGGrid<3>; };

// Specialize the problem type for this type tag
template<class TypeTag>
struct Problem<TypeTag, TTag::LandfillGridmanagerTest>
{ using type = LandfillGridmanagerTestProblem<TypeTag>; };

// Specialize the spatial params type for this type tag
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LandfillGridmanagerTest>
{
    using GridGeometry = GetPropType<TypeTag, GridGeometry>;
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = LandfillGridmanagerTestSpatialParams<GridGeometry, Scalar>;
};

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::LandfillGridmanagerTest> { static constexpr bool value = true; };

} // end namespace Properties

/*!
 * \ingroup OnePNCTests
 * \brief
 */
template <class TypeTag>
class LandfillGridmanagerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy some indices for convenience
    enum
    {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx,
    };

    // world dimension to access gravity vector
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    LandfillGridmanagerTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return initial_(globalPos); }

    // \}

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1e5+ 9.81*1000*(this->gridGeometry().bBoxMax()[dimWorld-1] - globalPos[dimWorld-1]); // initial condition for the pressure
        return priVars;
    }

    std::string name_;

    static constexpr Scalar eps_ = 1.0e-6;

};

} // end namespace Dumux

#endif
