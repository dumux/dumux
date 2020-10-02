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
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 */

#ifndef DUMUX_1PNI_CONDUCTION_PROBLEM_HH
#define DUMUX_1PNI_CONDUCTION_PROBLEM_HH

#include <cmath>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class OnePNIConductionProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePNIConduction { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePNIConductionBox { using InheritsFrom = std::tuple<OnePNIConduction, BoxModel>; };
struct OnePNIConductionCCTpfa { using InheritsFrom = std::tuple<OnePNIConduction, CCTpfaModel>; };
struct OnePNIConductionCCMpfa { using InheritsFrom = std::tuple<OnePNIConduction, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConduction> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConduction> { using type = OnePNIConductionProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConduction>
{
    using type = FluidSystems::OnePLiquid<GetPropType<TypeTag, Properties::Scalar>,
                                          Components::H2O<GetPropType<TypeTag, Properties::Scalar>> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConduction>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNISpatialParams<GridGeometry, Scalar>;
};
}

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
 * The results are compared to an analytical solution for a diffusion process:
  \f[
     T =T_{high} + (T_{init} - T_{high})erf \left(0.5\sqrt{\frac{x^2 S_{total}}{t \lambda_{eff}}}\right)
 \f]
 *
 * This problem uses the \ref OnePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1pniconduction -ParameterFile ./test_box1pniconduction.input</tt> or <br>
 * <tt>./test_cc1pniconduction -ParameterFile ./test_cc1pniconduction.input</tt>
 */
template <class TypeTag>
class OnePNIConductionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using IapwsH2O = Components::H2O<Scalar>;

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
    }

    //! Get the analytical temperature
    const std::vector<Scalar>& getExactTemperature()
    {
        return temperatureExact_;
    }

    //! Udpate the analytical temperature
    void updateExactTemperature(const SolutionVector& curSol, Scalar time)
    {
        const auto someElement = *(elements(this->gridGeometry().gridView()).begin());

        auto someElemSol = elementSolution(someElement, curSol, this->gridGeometry());
        const auto someInitSol = initialAtPos(someElement.geometry().center());

        auto someFvGeometry = localView(this->gridGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = IapwsH2O::liquidHeatCapacity(someInitSol[temperatureIdx], someInitSol[pressureIdx]);
        const auto densityS = volVars.solidDensity();
        const auto heatCapacityS = volVars.solidHeatCapacity();
        const auto storage = densityW*heatCapacityW*porosity + densityS*heatCapacityS*(1 - porosity);
        const auto effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars);
        using std::max;
        time = max(time, 1e-10);
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
               auto globalIdx = scv.dofIndex();
               const auto& globalPos = scv.dofPosition();
               using std::erf;
               using std::sqrt;
               temperatureExact_[globalIdx] = temperatureHigh_ + (someInitSol[temperatureIdx] - temperatureHigh_)
                                              *erf(0.5*sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));

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
};

} // end namespace Dumux

#endif
