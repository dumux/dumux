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
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */

#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/staticinteractionvolume.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "spatialparams.hh"

namespace Dumux {

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

template <class TypeTag>
class InjectionProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Injection { using InheritsFrom = std::tuple<TwoPTwoC>; };
struct InjectionBox { using InheritsFrom = std::tuple<Injection, BoxModel>; };
struct InjectionCCTpfa { using InheritsFrom = std::tuple<Injection, CCTpfaModel>; };
struct InjectionCCMpfa { using InheritsFrom = std::tuple<Injection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection> { using type = InjectionProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection>
{
    using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Injection> { static constexpr bool value = true; };

// Enable caching or not (reference solutions created without caching)
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };

// use the static interaction volume around interior vertices in the mpfa test
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::InjectionCCMpfa>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // structured two-d grid
    static constexpr int numIvScvs = 4;
    static constexpr int numIvScvfs = 4;

    // use the default traits
    using Traits = CCMpfaODefaultStaticInteractionVolumeTraits< NodalIndexSet, Scalar, numIvScvs, numIvScvfs >;
public:
    using type = CCMpfaOStaticInteractionVolume< Traits >;
};
} // end namespace Properties

/*!
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one (\f$ K=10e-12\f$) for \f$ y<22m\f$ and one with a lower
 * permeablility (\f$ K=10e-13\f$) in the rest of the domain.
 *
 * A mixture of Nitrogen and Water vapor, which is composed according to the
 * prevailing conditions (temperature, pressure) enters a water-filled aquifer.
 * This is realized with a solution-dependent Neumann boundary condition at the
 * right boundary (\f$ 5m<y<15m\f$). The aquifer is situated 2700m below sea level.
 * The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the lower permeable aquitard.
 *
 * The model is able to use either mole or mass fractions. The property useMoles
 * can be set to either true or false in the problem file. Make sure that the
 * according units are used in the problem set-up.
 * The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2c</tt> or
 * <tt>./test_cc2p2c</tt>
 */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // primary variable indices
    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    // phase presence
    enum { wPhaseOnly = Indices::firstPhaseOnly };

    // equation indices
    enum
    {
        contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + FluidSystem::N2Idx,
    };

    // phase indices
    enum
    {
        gasPhaseIdx = FluidSystem::N2Idx,
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    InjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");
        temperature_        = getParam<Scalar>("Problem.InitialTemperature");
        depthBOR_           = getParam<Scalar>("Problem.DepthBOR");
        name_               = getParam<std::string>("Problem.Name");

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole-fractions"<<std::endl;
        else
            std::cout<<"problem uses mass-fractions"<<std::endl;
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
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K]
     */
    Scalar temperature() const
    { return temperature_; }

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
        if (globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.ipGlobal();

        Scalar injectedPhaseMass = 1e-3;
        Scalar moleFracW = elemVolVars[scvf.insideScvIdx()].moleFraction(gasPhaseIdx, H2OIdx);
        if (globalPos[1] < 14 - eps_ && globalPos[1] > 6.5 - eps_)
        {
            values[contiN2EqIdx] = -(1-moleFracW)*injectedPhaseMass/FluidSystem::molarMass(N2Idx); //mole/(m^2*s) -> kg/(s*m^2)
            values[contiH2OEqIdx] = -moleFracW*injectedPhaseMass/FluidSystem::molarMass(H2OIdx); //mole/(m^2*s) -> kg/(s*m^2)
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * The internal method for the initial condition
     *
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);

        Scalar densityW = FluidSystem::H2O::liquidDensity(temperature_, 1e5);

        Scalar pl = 1e5 - densityW*this->spatialParams().gravity(globalPos)[1]*(depthBOR_ - globalPos[1]);
        Scalar moleFracLiquidN2 = pl*0.95/BinaryCoeff::H2O_N2::henry(temperature_);
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidN2;

        Scalar meanM =
            FluidSystem::molarMass(H2OIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(N2Idx)*moleFracLiquidN2;
        if(useMoles)
        {
            //mole-fraction formulation
            priVars[switchIdx] = moleFracLiquidN2;
        }
        else
        {
            //mass fraction formulation
            Scalar massFracLiquidN2 = moleFracLiquidN2*FluidSystem::molarMass(N2Idx)/meanM;
            priVars[switchIdx] = massFracLiquidN2;
        }
        priVars[pressureIdx] = pl;
        return priVars;
    }

    Scalar temperature_;
    Scalar depthBOR_;
    static constexpr Scalar eps_ = 1e-6;

    int nTemperature_;
    int nPressure_;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    std::string name_;
};

} // end namespace Dumux

#endif
