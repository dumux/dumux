// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 */

#ifndef DUMUX_1PNI_CONDUCTION_PROBLEM_HH
#define DUMUX_1PNI_CONDUCTION_PROBLEM_HH

#include <cmath>
#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidmatrixinteractions/thermalconductivityaverage.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side there is a Dirichlet boundary condition with an increased
 * temperature and on the right hand side a Dirichlet boundary with constant pressure,
 * saturation and temperature is applied.
 *
 * The results are compared to an analytical solution for a diffusion process.
 * This problem uses the \ref OnePModel and \ref NIModel model.
 */
template <class TypeTag>
class OnePNIConductionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    // the analytical solution uses the volume-averaged effective thermal conductivity
    static_assert(std::is_same_v<ThermalConductivityModel, ThermalConductivityAverage<Scalar>>,
                  "The analytical solution assumes ThermalConductivityAverage");

    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    OnePNIConductionProblem(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup)
    : ParentType(gridGeometry, paramGroup)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = getParam<std::string>("Problem.Name");
        temperatureHigh_ = 300.0;
        temperatureExact_.resize(gridGeometry->numDofs());

        // all fluid and solid properties are constant, so the total volumetric heat capacity
        // and the effective thermal conductivity are constant in space and time
        const auto initValues = initial_();
        const Scalar temperature = initValues[temperatureIdx];
        const Scalar pressure = initValues[pressureIdx];
        const Scalar porosity = this->spatialParams().porosityAtPos(gridGeometry->bBoxMin());

        totalHeatCapacity_ = porosity*FluidSystem::density(temperature, pressure)*FluidSystem::heatCapacity(temperature, pressure)
                             + (1.0 - porosity)*SolidSystem::density(temperature)*SolidSystem::heatCapacity(temperature);
        // same as in volvars effectiveThermalConductivity() method, but with const properties, this just needs to be computed once
        effectiveThermalConductivity_ = porosity*FluidSystem::thermalConductivity(temperature, pressure)
                                        + (1.0 - porosity)*SolidSystem::thermalConductivity(temperature);
    }

    //! Get the analytical temperature
    const std::vector<Scalar>& getExactTemperature()
    {
        return temperatureExact_;
    }

    /*!
     * \brief Update the analytical temperature
     * The results are compared to an analytical solution for a diffusion process:
      \f[
         T = T_{high} + (T_{init} - T_{high})\,\mathrm{erf} \left(0.5\sqrt{\frac{x^2 C_{tot}}{t \lambda_{eff}}}\right)
      \f]
     * The solution vector is not needed, as all properties are constant.
     */
    void updateExactTemperature(const SolutionVector&, Scalar time)
    {
        using std::max;
        time = max(time, 1e-10);
        const Scalar temperatureInit = initial_()[temperatureIdx];
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto x = scv.dofPosition()[0];
                using std::erf;
                using std::sqrt;
                temperatureExact_[scv.dofIndex()] = temperatureHigh_ + (temperatureInit - temperatureHigh_)
                                                    *erf(0.5*sqrt(x*x*totalHeatCapacity_/time/effectiveThermalConductivity_));
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
        BoundaryTypes bcTypes;

        if(globalPos[0] < eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
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
        PrimaryVariables priVars(initial_());
        if (globalPos[0] < eps_)
            priVars[temperatureIdx] = temperatureHigh_;
        return priVars;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

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
        return initial_();
    }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_() const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1.0e5;
        priVars[temperatureIdx] = 290.0;
        return priVars;
    }

    Scalar temperatureHigh_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    std::vector<Scalar> temperatureExact_;
    Scalar totalHeatCapacity_;
    Scalar effectiveThermalConductivity_;
};

} // end namespace Dumux

#endif
