// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem for the 1pnc problem:
 * Component transport of nitrogen dissolved in the water phase.
 */

#ifndef DUMUX_1P2C_TEST_PROBLEM_HH
#define DUMUX_1P2C_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>


namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem for the 1pnc problem:
 *  Component transport of nitrogen dissolved in the water phase.
 *
 * Nitrogen is dissolved in the water phase and is transported with the
 * water flow from the left side to the right.
 *
 * The model domain is specified in the input file and
 * we use homogeneous soil properties.
 * Initially, the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction.
 * The water phase flows from the left side to the right if the applied pressure
 * gradient is >0. The nitrogen is transported with the water flow
 * and leaves the domain at theboundary, where again Dirichlet boundary
 * conditions are applied.
 *
 * This problem uses the \ref OnePNCModel model.
 */
template <class TypeTag>
class OnePNCDispersionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
#if NONISOTHERMAL
        energyEqIdx = Indices::energyEqIdx,
#endif
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),
        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx
    };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int numFluidComps = FluidSystem::numComponents;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePNCDispersionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if constexpr (useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';

        boundaryConcentration_ = getParam<Scalar>("Problem.BoundaryConcentration");
        temperatureDifference_ = getParam<Scalar>("Problem.BoundaryTemperatureDifference");
        counterFlowRate_ = getParam<Scalar>("Problem.CounterFlowRate");
        pressure_ = getParam<Scalar>("Problem.Pressure");
        problemName_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();

        if (isRightBoundary_(globalPos) ||
            isLowerBoundary_(globalPos) ||
            isUpperBoundary_(globalPos) )
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        if (isLeftBoundary_(globalPos))
        {
            values[contiN2EqIdx] =  boundaryConcentration_;
#if NONISOTHERMAL
            values[energyEqIdx] = this->spatialParams().temperatureAtPos(globalPos) + temperatureDifference_;
#endif
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the box scheme).
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = element.geometry().corner(scvf.insideScvIdx());
        NumEqVector values(0.0);
        if (isRightBoundary_(globalPos))
            values[contiH2OEqIdx] = -1.0 * counterFlowRate_;
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilated per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions (mole/(m^3*s) or kg/(m^3*s)).
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    //#### initial value for a control volume
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        GetPropType<TypeTag, Properties::FluidState> fluidState;
        fluidState.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
        fluidState.setPressure(0, pressure_);
        Scalar density = FluidSystem::density(fluidState, 0);
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar gravity = this->spatialParams().gravity(globalPos)[1];

        values[Indices::pressureIdx] = pressure_ - (density * gravity * depth); //initial condition for the pressure
        values[contiN2EqIdx] = 0.0; //initial condition for the molefraction
#if NONISOTHERMAL
        values[energyEqIdx] = this->spatialParams().temperatureAtPos(globalPos);
#endif
        return values;
    }

private:
    bool isLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isRightBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isLowerBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool isUpperBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    static constexpr Scalar eps_ = 1e-6;
    Scalar pressure_;
    Scalar counterFlowRate_;
    Scalar temperatureDifference_;
    Scalar boundaryConcentration_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
