// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLE_BLOODFLOW_SOLVER_HH
#define DUMUX_EXAMPLE_BLOODFLOW_SOLVER_HH

#include <memory>

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/reorderingdofmapper.hh>

#include <dumux/io/format.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux::BloodFlowSolver {

template<class TypeTag>
class BloodFlowProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    static constexpr int dim = GridView::dimension;

public:
    template<class GridData>
    BloodFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                     std::shared_ptr<GridData> gridData)
    : ParentType(gridGeometry)
    {
        // read spatial parameters
        this->spatialParams().readGridParams(*gridData);

        // mapping from element index to boundaryFlag
        referencePressure_ = getParam<Scalar>("Problem.ReferencePressure", 0.0);
        setBoundaryDomainMarker_(*gridData);
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return { 0.0 }; }

    /////////////////////////////////////////////////////////////////////////////////
    /////// CC Tpfa /////////////////////////////////////////////////////////////////
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = this->bloodPressureBoundaryScvf(scvf.index());
        return values;
    }

    //! blood pressure for boundary sub control volume faces
    Scalar bloodPressureBoundaryScvf(std::size_t scvfIndex) const
    {
        const auto p = bloodPressureBoundary_[scvfIndex];
        if (std::isnan(p))
            DUNE_THROW(Dune::InvalidStateException, "Invalid scvf index for boundary pressure!");
        return p;
    }

    // compute volume fluxes as input for the transport model
    std::vector<Scalar> computeVolumeFluxes(const SolutionVector& sol, const GridVariables& gridVars)
    {
        const auto& gg = this->gridGeometry();
        std::vector<Scalar> volumeFlux(gg.numScvf(), 0.0);

        using FluxVariables =  GetPropType<TypeTag, Properties::FluxVariables>;
        auto upwindTerm = [](const auto& volVars) { return volVars.mobility(Indices::conti0EqIdx); };

        auto fvGeometry = localView(gg);
        auto elemVolVars = localView(gridVars.curGridVolVars());
        auto elemFluxVars = localView(gridVars.gridFluxVarsCache());
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVars.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto scvfIdx = scvf.index();

                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFlux[scvfIdx] = fluxVars.advectiveFlux(Indices::conti0EqIdx, upwindTerm);
                }
                else
                {
                    const auto bcTypes = boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet()) // Dirichlet
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                        volumeFlux[scvfIdx] = fluxVars.advectiveFlux(Indices::conti0EqIdx, upwindTerm);
                    }

                    else // Neumann
                    {
                        volumeFlux[scvfIdx] = this->neumann(element, fvGeometry, elemVolVars, 0.0, scvf)[Indices::pressureIdx]
                                              * scvf.area() * elemVolVars[0].extrusionFactor()
                                              / elemVolVars[0].density(); // volume flux from mass flux
                    }
                }
            }
        }

        return volumeFlux;
    }

private:
    template<class GridData>
    void setBoundaryDomainMarker_(const GridData& gridData)
    {
        const auto& gg = this->gridGeometry();

        // this is a bit difficult as we need to map from vertex to scvf.
        // maybe an scvf specialization for 1d could help containing the vertex index?
        bloodPressureBoundary_.resize(gg.numScvf(), std::numeric_limits<Scalar>::quiet_NaN());
        bloodPressureVertex_.resize(gg.gridView().size(GridView::dimension), std::numeric_limits<Scalar>::quiet_NaN());
        auto fvGeometry = localView(gg);
        std::size_t count = 0;
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& intersection : intersections(gg.gridView(), element))
            {
                if (intersection.boundary())
                {
                    // go up to the 0 level vertex
                    auto level0element = intersection.inside();
                    while (level0element.hasFather())
                        level0element = level0element.father();

                    const auto scvfIdxGlobal = [&]()
                    {
                        for (const auto& scvf : scvfs(fvGeometry))
                        {
                            if (scvf.boundary())
                                if (Dune::FloatCmp::eq(scvf.unitOuterNormal() * intersection.centerUnitOuterNormal(), 1.0, eps_))
                                    return scvf.index();
                        }

                        DUNE_THROW(Dune::InvalidStateException, "Boundary scvf not found!");

                    }();

                    const auto vertex = level0element.template subEntity<dim>(intersection.indexInInside());
                    bloodPressureVertex_[gg.vertexMapper().index(vertex)] = gridData.parameters(vertex)[0] - referencePressure_;
                    bloodPressureBoundary_[scvfIdxGlobal] = gridData.parameters(vertex)[0] - referencePressure_;
                    count++;
                }
            }
        }

        std::cout << Fmt::format("Extracted boundary data at {} boundary facets.\n", count);
    }

    static constexpr Scalar eps_ = 1e-7;
    Scalar referencePressure_;
    std::vector<Scalar> bloodPressureBoundary_; // blood pressure on boundary (Reichold/Schmid/Weber/Jenny)
    std::vector<Scalar> bloodPressureVertex_; // blood pressure on boundary (Reichold/Schmid/Weber/Jenny)
};

} //end namespace Dumux

namespace Dumux::BloodFlowSolver {

template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        referencePressure_ = getParam<Scalar>("Problem.ReferencePressure", 0.0);
    }

    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        const auto radius = this->radius(scv.elementIndex());
        return M_PI*radius*radius;
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_[scv.elementIndex()];
    }

    // Return the radius of the circular pipe
    Scalar radius(std::size_t eIdxGlobal) const
    { return radius_[eIdxGlobal];}

    //! Get the radii for e.g. output
    const std::vector<Scalar>& getRadii() const
    { return radius_; }

    //! Get the viscosity factors for e.g. output
    const std::vector<Scalar>& getViscosityFactors() const
    { return viscosityFactor_; }

    //! Get the boundary pressure at each vertex
    const std::vector<Scalar>& getBoundaryPressure() const
    { return boundaryPressure_; }

    //! Read params from dgf
    template<class GridData>
    void readGridParams(const GridData& gridData)
    {
        const auto& gg = this->gridGeometry();
        auto numElements = gg.gridView().size(0);
        radius_.resize(numElements);
        permeability_.resize(numElements);
        viscosityFactor_.resize(numElements);

        // gridView is a leafGridView. Parameters are only set on level 0.
        // elements have to inherit spatial parameters from their father.
        for (const auto& element : elements(gg.gridView()))
        {
            auto level0element = element;
            while (level0element.hasFather())
                level0element = level0element.father();

            auto eIdx = gg.elementMapper().index(element);
            const auto& params = gridData.parameters(level0element);
            radius_[eIdx] = params[0];

            constexpr Scalar dischargeHematocrit = 0.45;
            viscosityFactor_[eIdx] = etaVivo_(radius_[eIdx], dischargeHematocrit);
            constexpr Scalar gamma = 2; // quadratic velocity profile (Poiseuille flow)
            permeability_[eIdx] = radius_[eIdx]*radius_[eIdx]
                /(2.0*viscosityFactor_[eIdx]*(2.0+gamma));
        }

        boundaryPressure_.resize(gg.gridView().size(dim));
        for (const auto& vertex : vertices(gg.gridView()))
        {
            const auto vIdx = gg.vertexMapper().index(vertex);
            boundaryPressure_[vIdx] = gridData.parameters(vertex)[0] - referencePressure_;
        }

        for (const auto& element : elements(gg.gridView()))
            for (const auto& intersection : intersections(gg.gridView(), element))
                if (!intersection.boundary())
                    boundaryPressure_[
                        gg.vertexMapper().subIndex(element, intersection.indexInInside(), dim)
                    ] = 1e6;

    }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 37.0; } // Body temperature

private:
    Scalar etaVivo_(const Scalar radius, const Scalar dischargeHematocrit)
    {
        const auto D = 2*radius*1e6*1.187; // times rat blood factor (55fl RBC), convert to µm
        const auto scaling = (D/(D - 1.1))*(D/(D - 1.1));
        const auto C = corrFactor_(D);
        using std::exp; using std::pow;
        return scaling*(1.0 + scaling*(etaVivoH45_(D) - 1.0)*(pow(1.0 - dischargeHematocrit, C) - 1.0)/(pow(1.0 - 0.45, C) - 1.0));
    }

    Scalar etaVivoH45_(const Scalar D) // in µm
    {
        using std::exp; using std::pow;
        return 6.0*exp(-0.085*D) - 2.44*exp(-0.06*pow(D, 0.645)) + 3.2;
    }

    Scalar corrFactor_(const Scalar D) // in µm
    {
        using std::exp; using std::pow;
        const auto scaling = 1.0/(1.0 + 1e-11*pow(D, 12));
        return (0.8 + exp(-0.075*D))*(scaling - 1.0) + scaling;
    }

    std::vector<Scalar> radius_;
    std::vector<Scalar> permeability_;
    std::vector<Scalar> viscosityFactor_;
    std::vector<Scalar> boundaryPressure_;

    Scalar referencePressure_;
};

} // end namespace Dumux

namespace Dumux::Properties {
namespace TTag { struct BloodFlowCCTpfa { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; }; }

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = true; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::BloodFlowCCTpfa> { static constexpr bool value = false; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BloodFlowCCTpfa>
{ using type = Dune::FoamGrid<1, 3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::BloodFlowCCTpfa>
{ using type = BloodFlowSolver::BloodFlowProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BloodFlowCCTpfa>
{ using type = BloodFlowSolver::BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BloodFlowCCTpfa>
{ using type = FluidSystems::OnePLiquid<double, Components::Constant<0, double> >; };

// Set the local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BloodFlowCCTpfa>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

// if we have pt scotch use the reordering dof mapper to optimally sort the dofs (cc)
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::BloodFlowCCTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using ElementMapper = ReorderingDofMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using MapperTraits = DefaultMapperTraits<GridView, ElementMapper, VertexMapper>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>>;
};

} // end namespace Dumux::Properties

namespace Dumux {

template<class GridGeometry, class GridData>
std::vector<double> computeBloodVolumeFluxes(const GridGeometry& gridGeometry, const GridData& gridData)
{
    // Define the sub problem type tags
    using TypeTag = Properties::TTag::BloodFlowCCTpfa;

    // start the timer
    Dune::Timer timer;

    ////////////////////////////////////////////////////////////
    // Flow: prepare variables and problem
    ////////////////////////////////////////////////////////////
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, gridData);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector sol;
    problem->applyInitialSolution(sol);

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // debug output
    using VTKOutputModule = VtkOutputModule<GridVariables, SolutionVector>;
    std::unique_ptr<VTKOutputModule> vtkWriter;

    static const bool writeVTK = getParam<bool>("BloodFlow.WriteVTK", false);
    if (writeVTK)
    {
        vtkWriter = std::make_unique<VTKOutputModule>(*gridVariables, sol, problem->name());
        GetPropType<TypeTag, Properties::IOFields>::initOutputModule(*vtkWriter);
        vtkWriter->addVelocityOutput(
            std::make_shared<PorousMediumFlowVelocityOutput<
                GridVariables, GetPropType<TypeTag, Properties::FluxVariables>
            >>(*gridVariables)
        );
        vtkWriter->addField(problem->spatialParams().getRadii(), "radius");
        vtkWriter->addField(problem->spatialParams().getViscosityFactors(), "vessel tortuosity factor");
        vtkWriter->addField(problem->spatialParams().getBoundaryPressure(), "pBC");
        vtkWriter->write(0.0, Dune::VTK::OutputType::appendedraw);
    }

    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(sol);

    if (writeVTK)
        vtkWriter->write(1.0, Dune::VTK::OutputType::appendedraw);
    std::cout << "Flow computation took " << timer.elapsed() << " seconds." << std::endl;

    return problem->computeVolumeFluxes(sol, *gridVariables);
}

} // end namespace Dumux

#endif
