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
 * \brief Test for the OnePTwoCModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 */
#ifndef DUMUX_1P2CNI_CONVECTION_PROBLEM_HH
#define DUMUX_1P2CNI_CONVECTION_PROBLEM_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1p2cnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCNIConvectionProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCNIConvectionProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCNIConvectionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCNIConvectionProblem));
NEW_TYPE_TAG(OnePTwoCNIConvectionCCProblem, INHERITS_FROM(CCTpfaModel, OnePTwoCNIConvectionProblem));
NEW_TYPE_TAG(OnePTwoCNIConvectionCCMpfaProblem, INHERITS_FROM(CCMpfaModel, OnePTwoCNIConvectionProblem));

// Set the grid type
SET_TYPE_PROP(OnePTwoCNIConvectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(OnePTwoCNIConvectionProblem, Problem,
              OnePTwoCNIConvectionProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCNIConvectionProblem,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCNIConvectionProblem, SpatialParams, OnePTwoCNISpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCNIConvectionProblem, UseMoles, true);

}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the OnePTwoCModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side water is injected at a constant rate and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution where a retarded front velocity is calculated as follows:
  \f[
     v_{Front}=\frac{q S_{water}}{\phi S_{total}}
 \f]
 *
 * The result of the analytical solution is written into the vtu files.
 *
 * This problem uses the \ref OnePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1p2cniconvection -ParameterFile ./test_box1p2cniconvection.input</tt> or <br>
 * <tt>./test_cc1p2cniconvection -ParameterFile ./test_cc1p2cniconvection.input</tt>
 */
template <class TypeTag>
class OnePTwoCNIConvectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using IapwsH2O =  H2O<Scalar>;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum
    {
        phaseIdx = Indices::phaseIdx,
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };


    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static const bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static const int dofCodim = isBox ? dimWorld : 0;

public:
    OnePTwoCNIConvectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);
        darcyVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, DarcyVelocity);

        temperatureHigh_ = 291.0;
        temperatureLow_ = 290.0;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
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

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
        const auto storageW =  densityW*heatCapacityW*porosity;
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);
        std::cout << "storage: " << storageTotal << '\n';

        const Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);
        const Scalar retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;
        std::cout << "retarded velocity: " << retardedFrontVelocity << '\n';

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                auto dofPosition = scv.dofPosition();
                temperatureExact[dofIdxGlobal] = (dofPosition[0] < retardedFrontVelocity*time) ? temperatureHigh_ : temperatureLow_;
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

        if(globalPos[0] > this->bBoxMax()[0] - eps_)
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
        return initial_(globalPos);
    }

    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        if(scvf.ipGlobal()[0] < eps_)
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            values[pressureIdx] = -darcyVelocity_
                                   *volVars.molarDensity(phaseIdx);
            values[massOrMoleFracIdx] = -darcyVelocity_
                                         *volVars.molarDensity(phaseIdx)
                                         *volVars.moleFraction(phaseIdx, massOrMoleFracIdx);
            values[temperatureIdx] = -darcyVelocity_
                                      *volVars.density(phaseIdx)
                                      *IapwsH2O::liquidEnthalpy(temperatureHigh_, volVars.pressure(phaseIdx));
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

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
    {
        return initial_(globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 1e-10;  // initial condition for the N2 molefraction
        priVars[temperatureIdx] = temperatureLow_;
        return priVars;
    }

    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    int outputInterval_;
};

} //end namespace Dumux

#endif // DUMUX_1P2CNI_CONVECTION_PROBLEM_HH
