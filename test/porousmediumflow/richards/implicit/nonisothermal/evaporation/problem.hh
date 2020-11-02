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
/**
 * \file
 * \ingroup RichardsTests
 * \brief Test for the extended richards problem:
 * The simulation domain is a tube a constant evaporation rate is set at the top and the soil gradually dries out.
 */

#ifndef DUMUX_RICHARDS_EVAPORATION_PROBLEM_HH
#define DUMUX_RICHARDS_EVAPORATION_PROBLEM_HH

#include <cmath>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include "../spatialparams.hh"

namespace Dumux {

/**
 * \ingroup RichardsTests
 * \brief Test for the RichardsModel in combination with the NI model for evaporation.
 */
template <class TypeTag>
class RichardsNIEvaporationProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct RichardsNIEvaporation { using InheritsFrom = std::tuple<RichardsNI>; };
struct RichardsNIEvaporationBox { using InheritsFrom = std::tuple<RichardsNIEvaporation, BoxModel>; };
struct RichardsNIEvaporationCC { using InheritsFrom = std::tuple<RichardsNIEvaporation, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsNIEvaporation> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsNIEvaporation> { using type = RichardsNIEvaporationProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RichardsNIEvaporation> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsNIEvaporation>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsNISpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableWaterDiffusionInAir<TypeTag, TTag::RichardsNIEvaporation> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup RichardsTests
 *
 * \brief Test for the RichardsModel in combination with the NI model for evaporation

 * The result of the analytical solution is written into the vtu files.
 * This problem uses the \ref RichardsModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_boxRichardsnievaporation -ParameterFile ./test_boxRichardsnievaporation.input</tt> or <br>
 * <tt>./test_ccRichardsnievaporation -ParameterFile ./test_ccRichardsnievaporation.input</tt>
 */
template <class TypeTag>
class RichardsNIEvaporationProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using IapwsH2O = Components::H2O<Scalar>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum { dimWorld = GridView::dimensionworld };

    enum {
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    RichardsNIEvaporationProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // initialize fluid system
        FluidSystem::init();

        name_ =  getParam<std::string>("Problem.Name");
        pressure_ = 9.9e4;
        temperatureInitial_ = 291;
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
    const std::string& name() const
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
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if(globalPos[1] < eps_)
        {
            values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }
        return values;
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
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The subcontrolvolume face
     *  Negative values mean influx.
     */
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        Scalar boundaryLayerThickness = 0.0016;

        if(globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
        {
             values[conti0EqIdx] = 1e-3;
             values[energyEqIdx] = FluidSystem::enthalpy( volVars.fluidState(), FluidSystem::gasPhaseIdx) * values[conti0EqIdx];
             values[energyEqIdx] += FluidSystem::thermalConductivity(volVars.fluidState(), FluidSystem::gasPhaseIdx)
                                    * (volVars.temperature() - temperatureInitial_)/boundaryLayerThickness;
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Returns the reference pressure [Pa] of the nonwetting
     *        fluid phase within a finite volume.
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonwettingReferencePressure() const
    { return 1e5; };

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

    // \}

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::bothPhases);
        priVars[pressureIdx] = pressure_; // initial condition for the pressure
        priVars[temperatureIdx] = temperatureInitial_;
        return priVars;
    }

    Scalar temperatureInitial_;
    Scalar pressure_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};

} // end namespace Dumux

#endif
