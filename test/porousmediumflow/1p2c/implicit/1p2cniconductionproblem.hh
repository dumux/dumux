// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Test for the OnePTwoCModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 *
 */
#ifndef DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH
#define DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1p2cnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCNIConductionProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCNIConductionProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCNIConductionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCNIConductionProblem));
NEW_TYPE_TAG(OnePTwoCNIConductionCCProblem, INHERITS_FROM(CCTpfaModel, OnePTwoCNIConductionProblem));
NEW_TYPE_TAG(OnePTwoCNIConductionCCMpfaProblem, INHERITS_FROM(CCMpfaModel, OnePTwoCNIConductionProblem));

// Set the grid type
SET_TYPE_PROP(OnePTwoCNIConductionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(OnePTwoCNIConductionProblem, Problem,
              OnePTwoCNIConductionProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCNIConductionProblem,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCNIConductionProblem,
              SpatialParams,
              OnePTwoCNISpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCNIConductionProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCNIConductionProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCNIConductionProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the OnePTwoCModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side there is a Dirichlet boundary condition with an increased temperature and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution for a diffusion process:
  \f[
     T =T_{high} + (T_{init} - T_{high})erf \left(0.5\sqrt{\frac{x^2 S_{total}}{t \lambda_{eff}}}\right)
 \f]
 *
 * This problem uses the \ref OnePTwoCModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1p2cniconduction -ParameterFile ./test_box1p2cniconduction.input</tt> or <br>
 * <tt>./test_cc1p2cniconduction -ParameterFile ./test_cc1p2cniconduction.input</tt>
 */
template <class TypeTag>
class OnePTwoCNIConductionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx,

        // index of the equations
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static const bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static const int dofCodim = isBox ? dimWorld : 0;

public:
    OnePTwoCNIConductionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                  int, Problem, OutputInterval);

        temperatureHigh_ = 300.;
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeFinished() ||
            this->timeManager().willBeFinished();
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void addVtkOutputFields(VtkOutputModule<TypeTag>& outputModule) const
    {
        auto& temperatureExact = outputModule.createScalarField("temperatureExact", dofCodim);

        const auto someElement = *(elements(this->gridView()).begin());
        const auto someElemSol = this->model().elementSolution(someElement, this->model().curSol());
        const auto someInitSol = initial_(someElement.geometry().center());

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storage = densityW*heatCapacityW*porosity + densityS*heatCapacityS*(1 - porosity);
        const auto effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars, this->spatialParams(),
                                                                                                         someElement, someFvGeometry, someScv);
        Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto globalIdx = scv.dofIndex();
                const auto& globalPos = scv.dofPosition();

                temperatureExact[globalIdx] = temperatureHigh_ + (someInitSol[temperatureIdx] - temperatureHigh_)
                                              *std::erf(0.5*std::sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));
            }
        }
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
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if(globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);

        // condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_)
            values[temperatureIdx] = temperatureHigh_;

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
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
        priVars[pressureIdx] = 1e5; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 1e-5;  // initial condition for the N2 molefraction
        priVars[temperatureIdx] = 290.;
        return priVars;
    }

    Scalar temperatureHigh_;
    const Scalar eps_;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif // DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH
