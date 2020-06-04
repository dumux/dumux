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
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected in a unsaturated porous medium.
 *
 * Tests a 2p2c model assuming chemical nonequiibrium.
 */

#ifndef DUMUX_TWOPTWOC_NONEQUILIBRIUM_PROBLEM_HH
#define DUMUX_TWOPTWOC_NONEQUILIBRIUM_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dune/common/parametertreeparser.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class TwoPTwoCChemicalNonequilibriumProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TwoPTwoCChemicalNonequilibrium { using InheritsFrom = std::tuple<TwoPTwoCNINonEquil>; };
struct TwoPTwoCChemicalNonequilibriumBox { using InheritsFrom = std::tuple<TwoPTwoCChemicalNonequilibrium, BoxModel>; };
struct TwoPTwoCChemicalNonequilibriumCC { using InheritsFrom = std::tuple<TwoPTwoCChemicalNonequilibrium, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = TwoPTwoCChemicalNonequilibriumProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:

    using type = FluidSystems::H2OAir<Scalar,
                                      Components::TabulatedComponent<Components::H2O<Scalar>>,
                                      FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>,
                                      true /*useKelvinEquation*/>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>;
};

// decide which type to use for floating values (double / quad)
template<class TypeTag>
struct Scalar<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium> { using type = double; };
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
public:
    static const TwoPFormulation value = TwoPFormulation::p0s1;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableThermalNonEquilibrium<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ static constexpr bool value = false; };

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = FouriersLaw<TypeTag>; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = Dumux::EnergyLocalResidual<TypeTag>; };

} // end namespace Properties

/*!
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected in a unsaturated porous medium.
 *
 * This test compares a mpnc problem with a 2p2c problem.
 */
template <class TypeTag>
class TwoPTwoCChemicalNonequilibriumProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    enum { dofCodim = isBox ? dim : 0 };
public:
    TwoPTwoCChemicalNonequilibriumProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);
        name_ = getParam<std::string>("Problem.Name");
    }

    void setGridVariables(std::shared_ptr<GridVariables> gridVariables)
    { gridVariables_ = gridVariables; }

    const GridVariables& gridVariables() const
    { return *gridVariables_; }

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

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \name Boundary conditions
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
         if (globalPos[0] <  this->gridGeometry().bBoxMin()[0] + eps_)
            bcTypes.setAllDirichlet();
         else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet oundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = initial_(globalPos);
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * Negative values mean influx.
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        Scalar boundaryLayerThickness = 0.0016;
        //right side
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        {
            Scalar moleFracH2OInside = volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx);
            Scalar moleFracRefH2O = 0.0;
            Scalar evaporationRate = volVars.diffusionCoefficient(FluidSystem::gasPhaseIdx, FluidSystem::AirIdx, FluidSystem::H2OIdx)
                                     * (moleFracH2OInside - moleFracRefH2O)
                                     / boundaryLayerThickness
                                     * volVars.molarDensity(FluidSystem::gasPhaseIdx);
            values[Indices::conti0EqIdx+2] = evaporationRate;

            values[Indices::energyEqIdx] = FluidSystem::enthalpy(volVars.fluidState(), FluidSystem::gasPhaseIdx) * evaporationRate;
            values[Indices::energyEqIdx] += FluidSystem::thermalConductivity(volVars.fluidState(), FluidSystem::gasPhaseIdx)
                                            * (volVars.temperature() - temperature_)/boundaryLayerThickness;
        }
        return values;
    }

    // \}

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

    /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        vtk.addField(xEquilxwn_, "xEquil^Air_liq");
        vtk.addField(xEquilxnw_, "xEquil^H2O_gas");
    }

    void updateVtkFields(const SolutionVector& curSol)
    {
        const auto& gridView = this->gridGeometry().gridView();
        xEquilxwn_.resize(gridView.size(dofCodim));
        xEquilxnw_.resize(gridView.size(dofCodim));
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                xEquilxwn_[dofIdxGlobal] = volVars.xEquil(0,1);
                xEquilxnw_[dofIdxGlobal] = volVars.xEquil(1,0);
            }
        }
    }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5; // water pressure
        values[Indices::switchIdx] = 0.8; // gas saturation
        values[2] = 5e-4; // xwn higher than equil, equil is 3.4e-5
        values[3] = 1e-2; // xnw lower than 1.3e-2
        values[Indices::temperatureIdx] = temperature_;
        values.setState(Indices::bothPhases);

        return values;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;

    std::shared_ptr<GridVariables> gridVariables_;
    std::vector<double> xEquilxnw_;
    std::vector<double> xEquilxwn_;
};
} // end namespace Dumux

#endif
