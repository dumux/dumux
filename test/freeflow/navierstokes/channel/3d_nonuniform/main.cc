// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief 3D Channel flow test for the staggered grid (Navier-)Stokes model
 */

#include <config.h>

#include <iostream>
#include <random>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/math.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include "properties.hh"

#ifndef USE_STOKES_SOLVER
#define USE_STOKES_SOLVER 0
#endif

template<class Vector, class MomGG, class MassGG, class MomP, class MomIdx, class MassIdx>
auto dirichletDofs(std::shared_ptr<MomGG> momentumGridGeometry,
                   std::shared_ptr<MassGG> massGridGeometry,
                   std::shared_ptr<MomP> momentumProblem,
                   MomIdx momentumIdx, MassIdx massIdx)
{
    Vector dirichletDofs;
    dirichletDofs[momentumIdx].resize(momentumGridGeometry->numDofs());
    dirichletDofs[massIdx].resize(massGridGeometry->numDofs());
    dirichletDofs = 0.0;

    auto fvGeometry = localView(*momentumGridGeometry);
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        fvGeometry.bind(element);

        for (const auto& intersection : intersections(momentumGridGeometry->gridView(), element))
        {
            if(intersection.boundary() == false)
                continue;

            const auto bcTypes = momentumProblem->boundaryTypes(fvGeometry, intersection);
            if (bcTypes.hasDirichlet())
            {
                for (auto&& localDof : localDofs(fvGeometry, intersection))
                {
                        for (int i = 0; i < bcTypes.size(); ++i)
                            if (bcTypes.isDirichlet(i))
                                dirichletDofs[momentumIdx][localDof.index()][i] = 1.0;
                }
            }
        }
    }

    return dirichletDofs;
}

#if HAVE_DUNE_FOAMGRID
// Compute the wall shear stress
template<class VelSolutionVector, class GridVariables, class CouplingManager, class Assembler>
void computeWallShearStress(
    const VelSolutionVector& curSol,
    const GridVariables& gridVariables,
    CouplingManager& couplingManager,
    const Assembler& assembler
){
    const auto& problem = gridVariables.curGridVolVars().problem();
    const auto& gg = problem.gridGeometry();
    const auto& gv = gg.gridView();

    using GridView = typename GridVariables::GridGeometry::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename GridVariables::GridGeometry::GlobalCoordinate;

    using Grid = Dune::FoamGrid<dim-1, dimWorld>;
    Dune::GridFactory<Grid> factory;
    {
        std::vector<std::vector<std::vector<unsigned int>>> elems;
        elems.reserve(gg.numBoundaryScvf());
        std::vector<bool> insertedVertex(gv.size(dim), false);
        std::vector<int> vertexIndexMap(gv.size(dim), -1);
        std::vector<Dune::FieldVector<double, 3>> points;
        points.reserve(gg.numBoundaryScvf()*4);
        std::size_t boundaryElementIndex = 0;
        for (const auto& element : elements(gv))
        {
            for (const auto& intersection : intersections(gv, element))
            {
                if (intersection.boundary())
                {
                    const auto insideIdx = intersection.indexInInside();
                    const auto refElement = referenceElement(element);
                    const auto numVertices = refElement.size(insideIdx, 1, dim);
                    if (numVertices == 3)
                        elems.push_back({std::vector<unsigned int>{}});
                    else if (numVertices == 4)
                        elems.push_back({std::vector<unsigned int>{},
                                            std::vector<unsigned int>{}}); // add two triangles
                    else
                        DUNE_THROW(Dune::NotImplemented,
                            "Wall shear stress for boundary type with " << numVertices << " corners"
                        );


                    for (int i = 0; i < numVertices; ++i)
                    {
                        const auto localVIdx = refElement.subEntity(insideIdx, 1, i, dim);
                        const auto& vertex = element.template subEntity<dim>(localVIdx);
                        const auto vIdx = gg.vertexMapper().index(vertex);
                        if (!insertedVertex[vIdx])
                        {
                            vertexIndexMap[vIdx] = points.size();
                            points.push_back(vertex.geometry().corner(0));
                            insertedVertex[vIdx] = true;
                        }

                        if (numVertices == 3)
                            elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);

                        else if (numVertices == 4)
                        {
                            // add two triangles
                            if (i == 0)
                                elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);
                            else if (i == 1 || i == 2)
                            {
                                elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);
                                elems[boundaryElementIndex][1].push_back(vertexIndexMap[vIdx]);
                            }
                            else
                                elems[boundaryElementIndex][1].push_back(vertexIndexMap[vIdx]);
                        }
                        else
                            DUNE_THROW(Dune::NotImplemented,
                                "Wall shear stress for boundary type with " << numVertices << " corners"
                            );
                    }

                    ++boundaryElementIndex;
                }
            }
        }

        for (const auto& p : points)
            factory.insertVertex(p);
        for (const auto& e : elems)
            for (const auto& ee : e)
                factory.insertElement(Dune::GeometryTypes::simplex(dim-1), ee);
    }

    auto bGrid = factory.createGrid();

    using ElementSet = Dumux::GridViewGeometricEntitySet<typename Grid::LeafGridView, 0>;
    using Tree = Dumux::BoundingBoxTree<ElementSet>;
    Tree tree(std::make_shared<ElementSet>(bGrid->leafGridView()));

    std::vector<GlobalPosition> wallShearStress(bGrid->leafGridView().size(0));

    auto fvGeometry = localView(gg);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    for (const auto& element : elements(gv))
    {
        const auto h = Dumux::diameter(element.geometry());
        couplingManager.bindCouplingContext(
            CouplingManager::freeFlowMomentumIndex, element, assembler
        );
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, curSol);

        for (const auto& intersection : intersections(gv, element))
        {
            if (intersection.boundary())
            {

                for (const auto& qpData : Dumux::CVFE::quadratureRule(fvGeometry, intersection, Dumux::QuadratureRules::DuneQuadrature<2>()))
                {
                    const auto& ipData = qpData.ipData();

                    // assemble wss = sigma*n - (sigma*n*n)*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto& elementGeometry = fvGeometry.elementGeometry();
                    const auto gradV = evalGradientsAtLocalPos(element, elementGeometry, gg, elemSol, ipData.local());

                    auto D = gradV;
                    for(int i=0; i < gradV.size(); i++)
                        for(int j=0; j<gradV[i].size(); j++)
                            D[i][j] += gradV[j][i];

                    GlobalPosition sigman(0.0);
                    const auto normal = ipData.unitOuterNormal();
                    for(int i=0; i < normal.size(); i++)
                        sigman[i] = D[i]*normal;

                    sigman *= -problem.effectiveViscosity(element, fvGeometry, ipData);
                    sigman += problem.pressure(element, fvGeometry, ipData)*normal;

                    const auto normalStress = (sigman * normal)*normal;
                    const auto wss = sigman - normalStress;

                    // make sure we hit all triangles of this face
                    const auto points = Dumux::EmbeddedCoupling::circlePoints(ipData.global(), normal, 0.05*h, 4);
                    for (const auto& p : points)
                        for (auto b : intersectingEntities(p, tree))
                            wallShearStress[b] = wss;
                }
            }
        }
    }

    Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView> bMapper(
        bGrid->leafGridView(), Dune::mcmgLayout(Dune::Codim<0>())
    );
    using Field = Dumux::Vtk::Field<typename Grid::LeafGridView>;
    Dune::VTKWriter<typename Grid::LeafGridView> writer(bGrid->leafGridView());
    writer.addCellData(Field(
        bGrid->leafGridView(), bMapper, wallShearStress, "wallShearStress", dimWorld
    ).get());
    writer.write(problem.name() + "_wall_shear_stress");
}
#endif

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

    // initialize MPI+X
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // time this test
    Dune::Timer timer;

    // create a grid
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;
    gridManager.init();
    auto gridData = gridManager.getGridData();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager, gridData);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager, gridData);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    std::cout << "Total number of dofs: "
        << massGridGeometry->numDofs() + momentumGridGeometry->numDofs()*MomentumGridGeometry::GridView::dimension
        << std::endl;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // the linearize and solve
#if USE_STOKES_SOLVER
    using Matrix = typename Assembler::JacobianMatrix;
    using Vector = typename Assembler::ResidualType;
    using LinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
    auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
    auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
#else
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
#endif

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
    nonLinearSolver.solve(x);

    // write vtk output
    vtkWriter.write(1.0);

    // compute influx and outflux
    couplingManager->updateSolution(x);
    momentumProblem->computeFluxes(x[momentumIdx]);

#if HAVE_DUNE_FOAMGRID
    computeWallShearStress(x[momentumIdx], *momentumGridVariables, *couplingManager, *assembler);
#endif

    timer.stop();
    const auto& comm = leafGridView.comm();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
