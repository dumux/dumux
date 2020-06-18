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
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 */

#ifndef DUMUX_HETEROGENEOUS_PROBLEM_HH
#define DUMUX_HETEROGENEOUS_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/co2/model.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/discretization/box/scvftoscvboundarytypes.hh>

#include "spatialparams.hh"
#include "co2tables.hh"

// per default use isothermal model
#ifndef ISOTHERMAL
#define ISOTHERMAL 1
#endif

namespace Dumux {

template <class TypeTag>
class HeterogeneousProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Heterogeneous { using InheritsFrom = std::tuple<TwoPTwoCCO2>; };
struct HeterogeneousBox { using InheritsFrom = std::tuple<Heterogeneous, BoxModel>; };
struct HeterogeneousCCTpfa { using InheritsFrom = std::tuple<Heterogeneous, CCTpfaModel>; };
} // end namespace TTag

//Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Heterogeneous> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Heterogeneous> { using type = HeterogeneousProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Heterogeneous>
{
    using type = HeterogeneousSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Heterogeneous>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        HeterogeneousCO2Tables::CO2Tables,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
};

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Heterogeneous> { static constexpr bool value = false; };

#if !ISOTHERMAL
// Create new type tags
namespace TTag {
struct HeterogeneousNI { using InheritsFrom = std::tuple<TwoPTwoCCO2NI>; };
struct HeterogeneousNIBox { using InheritsFrom = std::tuple<HeterogeneousNI, BoxModel>; };
struct HeterogeneousNICCTpfa { using InheritsFrom = std::tuple<HeterogeneousNI, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::HeterogeneousNI> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::HeterogeneousNI> { using type = HeterogeneousProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HeterogeneousNI>
{
    using type = HeterogeneousSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HeterogeneousNI>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        HeterogeneousCO2Tables::CO2Tables,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
};

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::HeterogeneousNI> { static constexpr bool value = false; };
#endif
} // end namespace Properties

/*!
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 *
 * The domain is sized 200m times 100m and consists of four layers, a
 * permeable reservoir layer at the bottom, a barrier rock layer with reduced
 * permeability, another reservoir layer and at the top a barrier rock layer
 * with a very low permeablility.
 *
 * CO2 is injected at the permeable bottom layer from the left side.
 * The domain is initially filled with brine.
 *
 * The grid is unstructered and permeability and porosity for the elements are
 * read in from the grid file. The grid file also contains so-called boundary
 * IDs which can be used assigned during the grid creation in order to differentiate
 * between different parts of the boundary.
 * These boundary ids can be imported into the problem where the boundary
 * conditions can then be assigned accordingly.
 *
 * The model is able to use either mole or mass fractions. The property useMoles
 * can be set to either true or false in the problem file. Make sure that the
 * according units are used in the problem setup.
 * The default setting for useMoles is false.
 *
 * To run the simulation execute the following line in shell (works with the
 * box and cell centered spatial discretization method):
 * <tt>./test_ccco2 </tt> or <tt>./test_boxco2 </tt>
 */
template <class TypeTag >
class HeterogeneousProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // copy some indices for convenience
    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // phase presence index
        firstPhaseOnly = Indices::firstPhaseOnly,

        // component indices
        BrineIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        contiCO2EqIdx = conti0EqIdx + CO2Idx
    };

#if !ISOTHERMAL
    enum {
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
    };
#endif

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using CO2 = Components::CO2<Scalar, HeterogeneousCO2Tables::CO2Tables>;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

    // the discretization method we are using
    static constexpr auto discMethod = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;

    // world dimension to access gravity vector
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    template<class SpatialParams>
    HeterogeneousProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    , injectionTop_(1)
    , injectionBottom_(2)
    , dirichletBoundary_(3)
    , noFlowBoundary_(4)
    {
        nTemperature_       = getParam<int>("FluidSystem.NTemperature");
        nPressure_          = getParam<int>("FluidSystem.NPressure");
        pressureLow_        = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_       = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("FluidSystem.TemperatureHigh");
        depthBOR_           = getParam<Scalar>("Problem.DepthBOR");
        name_               = getParam<std::string>("Problem.Name");
        injectionRate_      = getParam<Scalar>("Problem.InjectionRate");

#if !ISOTHERMAL
        injectionPressure_ = getParam<Scalar>("Problem.InjectionPressure");
        injectionTemperature_ = getParam<Scalar>("Problem.InjectionTemperature");
#endif

        // set the spatial parameters by reading the DGF grid file
        this->spatialParams().getParamsFromGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        // stating in the console whether mole or mass fractions are used
        if (gridGeometry->gridView().comm().rank() == 0)
        {
            if (useMoles)
                std::cout << "-- Problem uses mole fractions." << std::endl;
            else
                std::cout << "-- Problem uses mass fractions." << std::endl;
        }

        // precompute the boundary types for the box method from the cell-centered boundary types
        scvfToScvBoundaryTypes_.computeBoundaryTypes(*this);
    }

    /*!
     * \brief Appends all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class VTKWriter>
    void addFieldsToWriter(VTKWriter& vtk)
    {
        const auto numElements = this->gridGeometry().gridView().size(0);
        const auto numDofs = this->gridGeometry().numDofs();

        vtkKxx_.resize(numElements);
        vtkPorosity_.resize(numElements);
        vtkBoxVolume_.resize(numDofs, 0.0);

        vtk.addField(vtkKxx_, "Kxx");
        vtk.addField(vtkPorosity_, "cellwisePorosity");
        vtk.addField(vtkBoxVolume_, "boxVolume");

#if !ISOTHERMAL
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(BrineIdx); }, "enthalpyW");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(CO2Idx); }, "enthalpyN");
#else
        vtkTemperature_.resize(numDofs, 0.0);
        vtk.addField(vtkTemperature_, "T");
#endif

        const auto& gridView = this->gridGeometry().gridView();
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                vtkBoxVolume_[dofIdxGlobal] += scv.volume();
#if ISOTHERMAL
                vtkTemperature_[dofIdxGlobal] = initialTemperatureField_(scv.dofPosition());
#endif
            }

            vtkKxx_[eIdx] = this->spatialParams().permeability(eIdx);
            vtkPorosity_[eIdx] = 1- this->spatialParams().inertVolumeFraction(eIdx);
        }

        // communicate box volume at process boundaries for vertex-centered scheme (box)
        if constexpr (isBox)
        {
            if (gridView.comm().size() > 1)
            {
                VectorCommDataHandleSum<typename GridGeometry::VertexMapper, std::vector<Scalar>, GridView::dimension>
                sumVolumeHandle(this->gridGeometry().vertexMapper(), vtkBoxVolume_);
                gridView.communicate(sumVolumeHandle,
                                     Dune::InteriorBorder_InteriorBorder_Interface,
                                     Dune::ForwardCommunication);
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
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The global position
     *
     * This problem assumes a geothermal gradient with
     * a surface temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return initialTemperatureField_(globalPos); }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    { return scvfToScvBoundaryTypes_.boundaryTypes(scv); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub-control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes bcTypes;
        const auto boundaryId = scvf.boundaryFlag();

        if (boundaryId < 1 || boundaryId > 4)
            DUNE_THROW(Dune::InvalidStateException, "Invalid boundary ID: " << boundaryId);

        if (boundaryId == dirichletBoundary_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \return the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

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
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
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
        const auto boundaryId = scvf.boundaryFlag();

        NumEqVector fluxes(0.0);
         // kg/(m^2*s) or mole/(m^2*s) depending on useMoles
        if (boundaryId == injectionBottom_)
        {
            fluxes[contiCO2EqIdx] = useMoles ? -injectionRate_/FluidSystem::molarMass(CO2Idx) : -injectionRate_;
#if !ISOTHERMAL
            // energy fluxes are always mass specific
            fluxes[energyEqIdx] = -injectionRate_/*kg/(m^2 s)*/*CO2::gasEnthalpy(
                                    injectionTemperature_, injectionPressure_)/*J/kg*/; // W/(m^2)
#endif
        }

        return fluxes;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values at a position.
     *
     * \param globalPos The global position
     * \return The initial values for the conservation equations in
     *           \f$ [ \textnormal{unit of primary variables} ] \f$
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

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
        PrimaryVariables values(0.0);
        values.setState(firstPhaseOnly);

        const Scalar temp = initialTemperatureField_(globalPos);
        const Scalar densityW = FluidSystem::Brine::liquidDensity(temp, 1e7);

        const Scalar moleFracLiquidCO2 = 0.00;
        const Scalar moleFracLiquidBrine = 1.0 - moleFracLiquidCO2;

        const Scalar meanM = FluidSystem::molarMass(BrineIdx)*moleFracLiquidBrine
                             + FluidSystem::molarMass(CO2Idx)*moleFracLiquidCO2;

        if(useMoles) // mole-fraction formulation
            values[switchIdx] = moleFracLiquidCO2;
        else // mass-fraction formulation
            values[switchIdx] = moleFracLiquidCO2*FluidSystem::molarMass(CO2Idx)/meanM;

        values[pressureIdx] = 1.0e5 - densityW*this->spatialParams().gravity(globalPos)[dimWorld-1]*(depthBOR_ - globalPos[dimWorld-1]);

#if !ISOTHERMAL
        values[temperatureIdx] = temp;
#endif
        return values;
    }

    Scalar initialTemperatureField_(const GlobalPosition globalPos) const
    {
        return 283.0 + (depthBOR_ - globalPos[dimWorld-1])*0.03;
    }

    Scalar depthBOR_;
    Scalar injectionRate_;

#if !ISOTHERMAL
    Scalar injectionPressure_, injectionTemperature_;
#endif

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    int injectionTop_;
    int injectionBottom_;
    int dirichletBoundary_;
    int noFlowBoundary_;

    // vtk output
    std::vector<Scalar> vtkKxx_, vtkPorosity_, vtkBoxVolume_, vtkTemperature_;
    ScvfToScvBoundaryTypes<BoundaryTypes, discMethod> scvfToScvBoundaryTypes_;
};

} // end namespace Dumux

#endif
