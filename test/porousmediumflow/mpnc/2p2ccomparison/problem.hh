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
 * \ingroup MPNCTests
 * \brief Problem where air is injected in a unsaturated porous medium.
 *
 * This test compares a mpnc problem with a 2p2c problem.
 */
#ifndef DUMUX_MPNC_TWOPTWOC_COMPARISON_OBSTACLEPROBLEM_HH
#define DUMUX_MPNC_TWOPTWOC_COMPARISON_OBSTACLEPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/mpnc/initialconditionhelper.hh>

namespace Dumux {

/*!
 * \ingroup MPNCTests
 * \brief Problem where air is injected in a unsaturated porous medium.
 *
 * This test compares a mpnc problem with a 2p2c problem.
 *
 */
template <class TypeTag>
class MPNCComparisonProblem
    : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NeumannFluxes = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using ParameterCache = typename FluidSystem::ParameterCache;

    static constexpr auto numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
    static constexpr auto numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();
    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr auto liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    static constexpr auto wCompIdx = FluidSystem::H2OIdx;
    static constexpr auto nCompIdx = FluidSystem::N2Idx;
    static constexpr auto fug0Idx = Indices::fug0Idx;
    static constexpr auto s0Idx = Indices::s0Idx;
    static constexpr auto p0Idx = Indices::p0Idx;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    static constexpr bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box;

public:
    MPNCComparisonProblem(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry)
    {
        const auto temperature = getParam<Scalar>("SpatialParams.Temperature");

        // initialize the tables of the fluid system
        Scalar Tmin = temperature - 1.0;
        Scalar Tmax = temperature + 1.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
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
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onOutlet_(globalPos))
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     * \param globalPos The global position
     */
    NeumannFluxes neumannAtPos(const GlobalPosition& globalPos) const
    {
        NeumannFluxes values(0.0);
        Scalar injectedAirMass = -1e-3;
        Scalar injectedAirMolarMass = injectedAirMass/FluidSystem::molarMass(FluidSystem::N2Idx);
        if (onInlet_(globalPos))
            values[Indices::conti0EqIdx + FluidSystem::N2Idx] = injectedAirMolarMass;
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
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
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        FluidState fs;

       // set the fluid temperatures
        fs.setTemperature(this->spatialParams().temperatureAtPos(globalPos));

        // set water saturation
        fs.setSaturation(liquidPhaseIdx, 0.8);
        fs.setSaturation(gasPhaseIdx, 1.0 - fs.saturation(liquidPhaseIdx));
        // set pressure of the gas phase
        fs.setPressure(gasPhaseIdx, 1e5);
        // calulate the capillary pressure
        const auto& fm =
            this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
        const int wPhaseIdx = this->spatialParams().template wettingPhaseAtPos<FluidSystem>(globalPos);
        const auto pc = fm.capillaryPressures(fs, wPhaseIdx);
        fs.setPressure(liquidPhaseIdx,
                       fs.pressure(gasPhaseIdx) + pc[liquidPhaseIdx] - pc[gasPhaseIdx]);

        ParameterCache paramCache;

        using InitialHelper = MPNCInitialConditionHelper<Scalar, PrimaryVariables, FluidSystem, GetPropType<TypeTag, Properties::ModelTraits>>;
        values = InitialHelper::solveForPrimaryVariables(fs, paramCache, MPNCInitialConditions::AllPhasesPresent{.refPhaseIdx = 0});

        return values;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10 + eps_;
    }

    bool onOutlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10 + eps_;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};
} // end namespace Dumux

#endif
