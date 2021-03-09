// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later vesion.                                      *
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
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal injection problem where water is injected into a
 *        sand column with a NAPL contamination.
 */

#ifndef DUMUX_COLUMNXYLOLPROBLEM_HH
#define DUMUX_COLUMNXYLOLPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux {
/*!
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal injection problem where water is injected into a
 *        sand column with a NAPL contamination.
 */
template <class TypeTag>
class ColumnProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Column { using InheritsFrom = std::tuple<ThreePThreeCNI>; };
struct ColumnBox { using InheritsFrom = std::tuple<Column, BoxModel>; };
struct ColumnCCTpfa { using InheritsFrom = std::tuple<Column, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Column> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Column> { using type = ColumnProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Column>
{ using type = FluidSystems::H2OAirXylene<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Column>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Component = Dumux::Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, Component>;
};

//! The two-phase model uses the immiscible fluid state
template<class TypeTag>
struct SolidState<TypeTag, TTag::Column>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = InertSolidState<Scalar, SolidSystem>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Column>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ColumnSpatialParams<GridGeometry, Scalar>;
};
} // end namespace Properties

/*!
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal injection problem where a water is injected into a
 *        sand column with a NAPL contamination.
 *
 * The 2D domain of this test problem is 0.1m long and 1.2m high.
 * Initially the column is filled with NAPL, gas and water, the NAPL saturation
 * increases to the bottom of the columns, the water saturation is constant.
 * Then water is injected from the top with a rate of 0.395710 mol/(s m) and
 * an enthalpy of 17452.97 [J/(s m)]. *
 *
 * Left, right and top boundaries are Neumann boundaries. Left and right are
 * no-flow boundaries, on top the injection takes places.
 * The bottom is a Dirichlet boundary.
 *
 * This problem uses the \ref ThreePThreeCModel and \ref NIModel model.
 *
 * This problem should typically be simulated for 200 days.
 * A good choice for the initial time step size is 1 s.
 * To adjust the simulation time it is necessary to edit the input file.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box3p3cnicolumnxylol test_columnxylol_fv.input</tt> or
 * <tt>./test_cc3p3cnicolumnxylol test_columnxylol_fv.input</tt>
 */
template <class TypeTag >
class ColumnProblem : public PorousMediumFlowProblem<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,

        // equation indices
        energyEqIdx = Indices::energyEqIdx,
        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wCompIdx, //!< Index of the mass conservation equation for the water component,
        contiGEqIdx = Indices::conti0EqIdx + FluidSystem::gCompIdx, //!< Index of the mass conservation equation for the gas component,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx, //!< Index of the mass conservation equation for the contaminant component

        // Phase State
        threePhases = Indices::threePhases
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    ColumnProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
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
    const std::string name() const
    { return name_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[1] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
    * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // negative values for injection
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
        {
            values[contiWEqIdx] = -0.395710;
            values[contiGEqIdx] = -0.000001;
            values[contiNEqIdx] = -0.00;
            values[energyEqIdx] = -17452.97;
        }
        return values;
    }

    // \}


    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }


    /*!
     * \brief Appends all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * Adjust this in case of anisotropic permeabilities.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        const auto& gg = this->gridGeometry();
        Kxx_.resize(gg.numDofs());
        vtk.addField(Kxx_, "permeability");

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(gg);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
                Kxx_[scv.dofIndex()] = this->spatialParams().intrinsicPermeabilityAtPos(scv.dofPosition());
        }
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values.setState(threePhases);
        const auto y = globalPos[1];
        const auto yMax = this->gridGeometry().bBoxMax()[1];

        values[temperatureIdx] = 296.15;
        values[pressureIdx] = 1.e5;
        values[switch1Idx] = 0.005;

        if (y > yMax - eps_)
            values[switch2Idx] = 0.112;
        else if (y > yMax - 0.0148 - eps_)
            values[switch2Idx] = 0 + ((yMax - y)/0.0148)*0.112;
        else if (y > yMax - 0.0296 - eps_)
            values[switch2Idx] = 0.112 + (((yMax - y) - 0.0148)/0.0148)*(0.120 - 0.112);
        else if (y > yMax - 0.0444 - eps_)
            values[switch2Idx] = 0.120 + (((yMax - y) - 0.0296)/0.0148)*(0.125 - 0.120);
        else if (y > yMax - 0.0592 - eps_)
            values[switch2Idx] = 0.125 + (((yMax - y) - 0.0444)/0.0148)*(0.137 - 0.125);
        else if (y > yMax - 0.0740 - eps_)
            values[switch2Idx] = 0.137 + (((yMax - y) - 0.0592)/0.0148)*(0.150 - 0.137);
        else if (y > yMax - 0.0888 - eps_)
            values[switch2Idx] = 0.150 + (((yMax - y) - 0.0740)/0.0148)*(0.165 - 0.150);
        else if (y > yMax - 0.1036 - eps_)
            values[switch2Idx] = 0.165 + (((yMax - y) - 0.0888)/0.0148)*(0.182 - 0.165);
        else if (y > yMax - 0.1184 - eps_)
            values[switch2Idx] = 0.182 + (((yMax - y) - 0.1036)/0.0148)*(0.202 - 0.182);
        else if (y > yMax - 0.1332 - eps_)
            values[switch2Idx] = 0.202 + (((yMax - y) - 0.1184)/0.0148)*(0.226 - 0.202);
        else if (y > yMax - 0.1480 - eps_)
            values[switch2Idx] = 0.226 + (((yMax - y) - 0.1332)/0.0148)*(0.257 - 0.226);
        else if (y > yMax - 0.1628 - eps_)
            values[switch2Idx] = 0.257 + (((yMax - y) - 0.1480)/0.0148)*(0.297 - 0.257);
        else if (y > yMax - 0.1776 - eps_)
            values[switch2Idx] = 0.297 + (((yMax - y) - 0.1628)/0.0148)*(0.352 - 0.297);
        else if (y > yMax - 0.1924 - eps_)
            values[switch2Idx] = 0.352 + (((yMax - y) - 0.1776)/0.0148)*(0.426 - 0.352);
        else if (y > yMax - 0.2072 - eps_)
            values[switch2Idx] = 0.426 + (((yMax - y) - 0.1924)/0.0148)*(0.522 - 0.426);
        else if (y > yMax - 0.2220 - eps_)
            values[switch2Idx] = 0.522 + (((yMax - y) - 0.2072)/0.0148)*(0.640 - 0.522);
        else if (y > yMax - 0.2368 - eps_)
            values[switch2Idx] = 0.640 + (((yMax - y) - 0.2220)/0.0148)*(0.767 - 0.640);
        else if (y > yMax - 0.2516 - eps_)
            values[switch2Idx] = 0.767 + (((yMax - y) - 0.2368)/0.0148)*(0.878 - 0.767);
        else if (y > yMax - 0.2664 - eps_)
            values[switch2Idx] = 0.878 + (((yMax - y) - 0.2516)/0.0148)*(0.953 - 0.878);
        else if (y > yMax - 0.2812 - eps_)
            values[switch2Idx] = 0.953 + (((yMax - y) - 0.2664)/0.0148)*(0.988 - 0.953);
        else if (y > yMax - 0.3000 - eps_)
            values[switch2Idx] = 0.988;
        else
            values[switch2Idx] = 1.e-4;
        return values;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    std::vector<Scalar> Kxx_;
};

} // end namespace Dumux

#endif
