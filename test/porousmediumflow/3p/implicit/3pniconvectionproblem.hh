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
 * \brief Definition of a 1p2cni problem:
 *        Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_3PNI_CONVECTION_PROBLEM_HH
#define DUMUX_3PNI_CONVECTION_PROBLEM_HH

#include <math.h>

#include <dumux/porousmediumflow/3p/implicit/model.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>
#include "3pnispatialparams.hh"


namespace Dumux
{

template <class TypeTag>
class ThreePNIConvectionProblem;

namespace Properties
{
NEW_TYPE_TAG(ThreePNIConvectionProblem, INHERITS_FROM(ThreePNI, ThreePNISpatialParams));
NEW_TYPE_TAG(ThreePNIConvectionBoxProblem, INHERITS_FROM(BoxModel, ThreePNIConvectionProblem));
NEW_TYPE_TAG(ThreePNIConvectionCCProblem, INHERITS_FROM(CCTpfaModel, ThreePNIConvectionProblem));
NEW_TYPE_TAG(ThreePNIConvectionCCMpfaProblem, INHERITS_FROM(CCMpfaModel, ThreePNIConvectionProblem));

// Set the grid type
SET_TYPE_PROP(ThreePNIConvectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ThreePNIConvectionProblem, Problem, ThreePNIConvectionProblem<TypeTag>);


// Set the fluid system
SET_TYPE_PROP(ThreePNIConvectionProblem,
              FluidSystem,
              FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the spatial parameters
SET_TYPE_PROP(ThreePNIConvectionProblem,
              SpatialParams,
              ThreePNISpatialParams<TypeTag>);
}


/*!
 * \ingroup ThreePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the ThreePModel in combination with the NI model for a convection problem:
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
 * This problem uses the \ref ThreePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box3pcniconvection -ParameterFile ./test_box3pniconvection.input</tt> or <br>
 * <tt>./test_cc3pcniconvection -ParameterFile ./test_cc3pniconvection.input</tt>
 */
template <class TypeTag>
class ThreePNIConvectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using IapwsH2O = H2O<Scalar>;
    using VtkOutputModule = typename GET_PROP_TYPE(TypeTag, VtkOutputModule);

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dimWorld : 0 };

    enum {
        // index of the primary variables
        pressureIdx = Indices::pressureIdx,
        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,
        temperatureIdx = Indices::temperatureIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;


public:
    ThreePNIConvectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                int, Problem, OutputInterval);
        darcyVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                Scalar, Problem, DarcyVelocity);

        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
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
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& temperatureExact = outputModule.createScalarField("temperatureExact", dofCodim);

        const auto someElement = *(elements(this->gridView()).begin());
        const auto someElemSol = this->model().elementSolution(someElement, this->model().curSol());
        const auto someInitSol = initialAtPos(someElement.geometry().center());

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density(wPhaseIdx);
        const auto heatCapacityW = IapwsH2O::liquidHeatCapacity(someInitSol[temperatureIdx], someInitSol[pressureIdx]);
        const auto storageW =  densityW*heatCapacityW*porosity;
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);
        std::cout << "storage: " << storageTotal << '\n';

        using std::max;
        const Scalar time = max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);
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
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param scvf The subcontrolvolume face
     *  Negative values mean influx.
     */
    PrimaryVariables neumann(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        const auto globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

        if(globalPos[0] < eps_)
        {
            values[conti0EqIdx] = -darcyVelocity_*volVars.density(wPhaseIdx);
            values[energyEqIdx] = -darcyVelocity_*volVars.density(wPhaseIdx)
                                     *IapwsH2O::liquidEnthalpy(temperatureHigh_, volVars.pressure(wPhaseIdx));
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position for which the source should be evaluated
     *
     * Returns the rate mass of a component is generated or annihilate
     * per volume unit. Positive values mean that mass is created,
     * negative ones mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = pressureLow_; // initial condition for the pressure
        values[swIdx] = 1.0;  // initial condition for the wetting phase saturation
        values[snIdx] = 1e-10;  // initial condition for the non-wetting phase saturation
        values[temperatureIdx] = temperatureLow_;
        return values;
    }

    // \}

private:
    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif
