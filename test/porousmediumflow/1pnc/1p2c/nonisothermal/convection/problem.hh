// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Test for the OnePNCModel in combination with the NI model for a conduction problem.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 */

#ifndef DUMUX_1P2CNI_CONVECTION_TEST_PROBLEM_HH
#define DUMUX_1P2CNI_CONVECTION_TEST_PROBLEM_HH
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/components/h2o.hh>
namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Test for the OnePTwoCModel in combination with the NI model for a convection problem.
 *
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially, the domain is fully saturated with water at a constant temperature.
 * On the left hand side water is injected at a constant rate and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The result of the analytical solution is written into the vtu files.
 *
 * This problem uses the \ref OnePModel and \ref NIModel model.
 */

template <class TypeTag>
class OnePTwoCNIConvectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IapwsH2O = Components::H2O<Scalar>;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),

        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx,
        energyEqIdx = Indices::energyEqIdx
    };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePTwoCNIConvectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        temperatureExact_.resize(gridGeometry->numDofs(), 290.0);

        darcyVelocity_ = getParam<Scalar>("Problem.DarcyVelocity");

        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
    }

    //! Get the analytical temperature
    const std::vector<Scalar>& getExactTemperature()
    {
        return temperatureExact_;
    }

    /*!
    * \brief Update the analytical temperature
    * The results are compared to an analytical solution where a retarded front velocity is calculated as follows:
     \f[
        v_{Front}=\frac{q S_{water}}{\phi S_{total}}
     \f]
    */
    void updateExactTemperature(const SolutionVector& curSol, Scalar time)
    {
        const auto someElement = *(elements(this->gridGeometry().gridView()).begin());

        auto someElemSol = elementSolution(someElement, curSol, this->gridGeometry());
        const auto someInitSol = initialAtPos(someElement.geometry().center());

        const auto someFvGeometry = localView(this->gridGeometry()).bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = IapwsH2O::liquidHeatCapacity(someInitSol[temperatureIdx], someInitSol[pressureIdx]);
        const auto storageW =  densityW*heatCapacityW*porosity;
        const auto densityS = volVars.solidDensity();
        const auto heatCapacityS = volVars.solidHeatCapacity();
        const auto storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);
        std::cout << "storage: " << storageTotal << '\n';

        using std::max;
        time = max(time, 1e-10);
        const Scalar retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;
        std::cout << "retarded velocity: " << retardedFrontVelocity << '\n';

        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                auto dofPosition = scv.dofPosition();
                temperatureExact_[dofIdxGlobal] = (dofPosition[0] < retardedFrontVelocity*time) ? temperatureHigh_ : temperatureLow_;
            }
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

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
        BoundaryTypes values;

        if(globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
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
     * \brief Evaluates the boundary conditions for a Dirichlet Sboundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache The cache related to flux computation
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        const auto& globalPos = scvf.ipGlobal();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if(globalPos[0] < eps_)
        {
             flux[contiH2OEqIdx] = -darcyVelocity_*elemVolVars[scv].molarDensity();
             flux[contiN2EqIdx] = -darcyVelocity_*elemVolVars[scv].molarDensity()*elemVolVars[scv].moleFraction(0, N2Idx);
             flux[energyEqIdx] = -darcyVelocity_
                                 *elemVolVars[scv].density()
                                 *IapwsH2O::liquidEnthalpy(temperatureHigh_, elemVolVars[scv].pressure());
        }

        return flux;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

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

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}
private:

    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[N2Idx] = 1e-10;  // initial condition for the N2 molefraction
        priVars[temperatureIdx] = temperatureLow_;
        return priVars;
    }
        static constexpr Scalar eps_ = 1e-6;
        Scalar temperatureHigh_;
        Scalar temperatureLow_;
        Scalar pressureHigh_;
        Scalar pressureLow_;
        Scalar darcyVelocity_;

        std::vector<Scalar> temperatureExact_;
    };

} // end namespace Dumux

#endif
