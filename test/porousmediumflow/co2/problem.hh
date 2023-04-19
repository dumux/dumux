// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected into a reservoir.
 */

#ifndef DUMUX_HETEROGENEOUS_PROBLEM_HH
#define DUMUX_HETEROGENEOUS_PROBLEM_HH

#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/parallel/vectorcommdatahandle.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/discretization/box/scvftoscvboundarytypes.hh>
#include <dumux/common/gridcapabilities.hh>


namespace Dumux {

/*!
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected into a reservoir.
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
 * IDs which can be used during the grid creation in order to differentiate
 * between different parts of the boundary.
 * These boundary IDs can be imported into the problem where the boundary
 * conditions can then be assigned accordingly.
 *
 * The model is able to use either mole or mass fractions. The property useMoles
 * can be set to either true or false in the problem file. Make sure that the
 * according units are used in the problem setup.
 * The default setting for useMoles is false.
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
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using CO2 = typename FluidSystem::CO2;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

    // the discretization method we are using
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;

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
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                vtkBoxVolume_[dofIdxGlobal] += scv.volume();
#if ISOTHERMAL
                vtkTemperature_[dofIdxGlobal] = this->spatialParams().temperatureAtPos(scv.dofPosition());
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
                if constexpr (Detail::canCommunicate<typename GridView::Traits::Grid, GridView::dimension>)
                    gridView.communicate(sumVolumeHandle,
                                         Dune::InteriorBorder_InteriorBorder_Interface,
                                         Dune::ForwardCommunication);
                else
                    DUNE_THROW(Dune::InvalidStateException, "Cannot call addFieldsToWriter on multiple processes for a grid that cannot communicate codim-" << GridView::dimension << "-entities.");
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

        const Scalar temp = this->spatialParams().temperatureAtPos(globalPos);
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
    ScvfToScvBoundaryTypes<BoundaryTypes, DiscretizationMethod> scvfToScvBoundaryTypes_;
};

} // end namespace Dumux

#endif
