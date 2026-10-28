// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>

#include <dune/common/indices.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#if DUMUX_HAVE_MUMPS
#include <dumux/linear/mumpssolver.hh>
#endif

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/solidmechanics/plate/foeppl_von_karman/localresidual.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include <dumux/solidmechanics/plate/foeppl_von_karman/holegauge.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief Deflection sampled along a line through the sheet
 *
 * \param scanAxis coordinate along which the line runs
 * \param station position of the line along the other coordinate
 */
template<class GridGeometry, class SolutionVector>
std::vector<std::pair<double, double>>
lineProfile(const GridGeometry& gridGeometry, const SolutionVector& x,
            int scanAxis, double station, double tolerance, int deformationIdx)
{
    const auto fixedAxis = 1 - scanAxis;
    std::vector<std::pair<double, double>> profile;
    for (const auto& vertex : vertices(gridGeometry.gridView()))
    {
        const auto pos = vertex.geometry().center();
        if (std::abs(pos[fixedAxis] - station) < tolerance)
            profile.emplace_back(pos[scanAxis], x[gridGeometry.vertexMapper().index(vertex)][deformationIdx]);
    }
    std::sort(profile.begin(), profile.end());
    return profile;
}

/*!
 * \brief Positions of the wrinkle crests along a scan line
 *
 * A crest is a local maximum of the deflection, which is how the wrinkles are counted
 * in the experiment. Maxima below the threshold are discarded so that numerical ripple
 * is not counted.
 */
std::vector<double> crestPositions(const std::vector<std::pair<double, double>>& profile,
                                   double threshold)
{
    std::vector<double> crests;
    for (std::size_t i = 1; i + 1 < profile.size(); ++i)
        if (profile[i].second > profile[i-1].second && profile[i].second >= profile[i+1].second
            && profile[i].second > threshold)
            crests.push_back(profile[i].first);
    return crests;
}

//! Half the peak-to-peak deflection of a scan line between two positions
double profileAmplitude(const std::vector<std::pair<double, double>>& profile, double lower, double upper)
{
    double minW = 0.0, maxW = 0.0;
    for (const auto& [s, w] : profile)
        if (s > lower && s < upper)
        {
            minW = std::min(minW, w);
            maxW = std::max(maxW, w);
        }
    return 0.5*(maxW - minW);
}

/*!
 * \brief Bending and membrane strain energy of the sheet
 *
 * One-point quadrature per element, exact for the piecewise linear parts of the fields;
 * the bubble of the rotations has no gradient at the element centre.
 */
template<class RotationGridGeometry, class DeformationGridGeometry, class InPlaneGridGeometry,
         class SolutionVector, class Sheet, std::size_t r, std::size_t d, std::size_t i>
std::array<double, 2> strainEnergy(const RotationGridGeometry& rotationGridGeometry,
                                   const DeformationGridGeometry& deformationGridGeometry,
                                   const InPlaneGridGeometry& inPlaneGridGeometry,
                                   const SolutionVector& x, const Sheet& sheet,
                                   Dune::index_constant<r> rotationIdx,
                                   Dune::index_constant<d> deformationIdx,
                                   Dune::index_constant<i> inPlaneIdx)
{
    const auto D = sheet.bendingStiffness();
    const auto nu = sheet.poissonRatio();
    double bending = 0.0, membrane = 0.0;
    for (const auto& element : elements(rotationGridGeometry.gridView()))
    {
        const auto geometry = element.geometry();
        const auto center = geometry.center();
        const auto gradRotation = evalGradients(
            element, geometry, rotationGridGeometry,
            elementSolution(element, x[rotationIdx], rotationGridGeometry), center
        );
        const auto gradU = evalGradients(
            element, geometry, inPlaneGridGeometry,
            elementSolution(element, x[inPlaneIdx], inPlaneGridGeometry), center
        );
        const auto gradW = evalGradients(
            element, geometry, deformationGridGeometry,
            elementSolution(element, x[deformationIdx], deformationGridGeometry), center
        )[FoepplVonKarmanPlateIndices::verticalDeformationIdx];

        const auto N = vonKarmanStressResultant(gradU, gradW, sheet.membraneStiffness(), nu);
        double curvatureNorm = 0.0, curvatureTrace = 0.0, membraneDensity = 0.0;
        for (int a = 0; a < 2; ++a)
        {
            curvatureTrace += gradRotation[a][a];
            for (int b = 0; b < 2; ++b)
            {
                const auto curvature = 0.5*(gradRotation[a][b] + gradRotation[b][a]);
                curvatureNorm += curvature*curvature;
                const auto strain = 0.5*(gradU[a][b] + gradU[b][a] + gradW[a]*gradW[b]);
                membraneDensity += 0.5*N[a][b]*strain;
            }
        }
        const auto volume = geometry.volume();
        bending += 0.5*D*((1.0 - nu)*curvatureNorm + nu*curvatureTrace*curvatureTrace)*volume;
        membrane += membraneDensity*volume;
    }
    return {bending, membrane};
}

//! Mean spacing of the crests in the central half of the scan line
double meanCrestSpacing(const std::vector<double>& crests, double lower, double upper)
{
    std::vector<double> central;
    std::copy_if(crests.begin(), crests.end(), std::back_inserter(central),
                 [&](double s){ return s > lower && s < upper; });
    if (central.size() < 2)
        return 0.0;
    return (central.back() - central.front())/(central.size() - 1);
}

//! In-plane stress resultants evaluated at the element centres
template<class InPlaneGridGeometry, class DeformationGridGeometry, class SolutionVector, class Sheet,
         std::size_t i, std::size_t j>
std::array<std::vector<double>, 3>
stressResultants(const InPlaneGridGeometry& inPlaneGridGeometry,
                 const DeformationGridGeometry& deformationGridGeometry,
                 const SolutionVector& x, const Sheet& sheet,
                 Dune::index_constant<i> inPlaneIdx, Dune::index_constant<j> deformationIdx)
{
    const auto numElements = inPlaneGridGeometry.gridView().size(0);
    std::array<std::vector<double>, 3> resultants{
        std::vector<double>(numElements), std::vector<double>(numElements), std::vector<double>(numElements)
    };

    for (const auto& element : elements(inPlaneGridGeometry.gridView()))
    {
        const auto geometry = element.geometry();
        const auto center = geometry.center();
        const auto gradU = evalGradients(
            element, geometry, inPlaneGridGeometry,
            elementSolution(element, x[inPlaneIdx], inPlaneGridGeometry), center
        );
        const auto gradW = evalGradients(
            element, geometry, deformationGridGeometry,
            elementSolution(element, x[deformationIdx], deformationGridGeometry), center
        )[FoepplVonKarmanPlateIndices::verticalDeformationIdx];

        const auto N = vonKarmanStressResultant(
            gradU, gradW, sheet.membraneStiffness(), sheet.poissonRatio()
        );

        const auto eIdx = inPlaneGridGeometry.elementMapper().index(element);
        resultants[0][eIdx] = N[0][0];
        resultants[1][eIdx] = N[1][1];
        resultants[2][eIdx] = N[0][1];
    }
    return resultants;
}

/*!
 * \brief Orientation of the wrinkles from the structure tensor of the slope
 *
 * Over the central half of the sheet, \f$ \int\nabla w\otimes\nabla w\,dA \f$ has its
 * largest eigenvector along the wrinkle normal. Returns the angle of that normal to the
 * \f$ x \f$ axis in degrees, in \f$ (-90, 90] \f$.
 */
template<class DeformationGridGeometry, class SolutionVector, std::size_t j>
double wrinkleNormalAngle(const DeformationGridGeometry& gridGeometry, const SolutionVector& x,
                          double lower, double upper, Dune::index_constant<j> deformationIdx)
{
    Dune::FieldMatrix<double, 2, 2> tensor(0.0);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto geometry = element.geometry();
        const auto center = geometry.center();
        if (center[0] < lower || center[0] > upper)
            continue;
        const auto gradW = evalGradients(
            element, geometry, gridGeometry, elementSolution(element, x[deformationIdx], gridGeometry), center
        )[FoepplVonKarmanPlateIndices::verticalDeformationIdx];
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                tensor[a][b] += gradW[a]*gradW[b]*geometry.volume();
    }
    // eigenvector of the larger eigenvalue of the symmetric 2x2 tensor
    const auto angle = 0.5*std::atan2(2.0*tensor[0][1], tensor[0][0] - tensor[1][1]);
    return angle*180.0/M_PI;
}

//! Area-weighted mean of an element field over the central half of the sheet
template<class GridGeometry>
double centralMean(const GridGeometry& gridGeometry, const std::vector<double>& field,
                   double lower, double upper)
{
    double sum = 0.0, area = 0.0;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto geometry = element.geometry();
        if (geometry.center()[0] < lower || geometry.center()[0] > upper)
            continue;
        sum += field[gridGeometry.elementMapper().index(element)]*geometry.volume();
        area += geometry.volume();
    }
    return sum/area;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using RotationTypeTag = Properties::TTag::WrinklingTestRotation;
    using DeformationTypeTag = Properties::TTag::WrinklingTestDeformation;
    using InPlaneTypeTag = Properties::TTag::WrinklingTestInPlane;
    using CommonTypeTag = Properties::TTag::WrinklingTestCommon;

    using Grid = GetPropType<CommonTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using RotationGridGeometry = GetPropType<RotationTypeTag, Properties::GridGeometry>;
    using DeformationGridGeometry = GetPropType<DeformationTypeTag, Properties::GridGeometry>;
    using InPlaneGridGeometry = GetPropType<InPlaneTypeTag, Properties::GridGeometry>;
    auto rotationGridGeometry = std::make_shared<RotationGridGeometry>(leafGridView);
    auto deformationGridGeometry = std::make_shared<DeformationGridGeometry>(leafGridView);
    auto inPlaneGridGeometry = std::make_shared<InPlaneGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<RotationTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(
        rotationGridGeometry, deformationGridGeometry, inPlaneGridGeometry
    );

    using Scalar = GetPropType<DeformationTypeTag, Properties::Scalar>;
    using GlobalPosition = typename Grid::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    auto sheet = std::make_shared<WrinklingSheet<Scalar, GlobalPosition>>();

    const auto loads = getParam<std::vector<Scalar>>("Problem.LoadSteps");
    sheet->setLoad(loads.front());

    using RotationProblem = GetPropType<RotationTypeTag, Properties::Problem>;
    auto rotationProblem = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManager, sheet);

    using DeformationProblem = GetPropType<DeformationTypeTag, Properties::Problem>;
    auto deformationProblem = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManager, sheet);

    using InPlaneProblem = GetPropType<InPlaneTypeTag, Properties::Problem>;
    auto inPlaneProblem = std::make_shared<InPlaneProblem>(inPlaneGridGeometry, couplingManager, sheet);

    constexpr auto rotationIdx = CouplingManager::rotationIdx;
    constexpr auto deformationIdx = CouplingManager::deformationIdx;
    constexpr auto inPlaneIdx = CouplingManager::inPlaneIdx;
    using Traits = MultiDomainTraits<RotationTypeTag, DeformationTypeTag, InPlaneTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;

    SolutionVector x;
    rotationProblem->applyInitialSolution(x[rotationIdx]);
    deformationProblem->applyInitialSolution(x[deformationIdx]);
    inPlaneProblem->applyInitialSolution(x[inPlaneIdx]);

    // the rotations start from the slope of the seed at the vertices and the bubbles from
    // zero; the dofs are addressed through the sub-control volumes because the bubble
    // scheme numbers its element dofs before its vertex dofs
    {
        x[rotationIdx] = 0.0;
        auto fvGeometry = localView(*rotationGridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (scv.localDofIndex() >= element.geometry().corners())
                    continue;
                const auto rotation = sheet->seedRotation(scv.dofPosition());
                for (int dir = 0; dir < GlobalPosition::dimension; ++dir)
                    x[rotationIdx][scv.dofIndex()][dir] = rotation[dir];
            }
        }
    }

    using RotationGridVariables = GetPropType<RotationTypeTag, Properties::GridVariables>;
    auto rotationGridVariables = std::make_shared<RotationGridVariables>(rotationProblem, rotationGridGeometry);

    using DeformationGridVariables = GetPropType<DeformationTypeTag, Properties::GridVariables>;
    auto deformationGridVariables = std::make_shared<DeformationGridVariables>(deformationProblem, deformationGridGeometry);

    using InPlaneGridVariables = GetPropType<InPlaneTypeTag, Properties::GridVariables>;
    auto inPlaneGridVariables = std::make_shared<InPlaneGridVariables>(inPlaneProblem, inPlaneGridGeometry);

    couplingManager->init(rotationProblem, deformationProblem, inPlaneProblem, x);
    std::vector<std::size_t> secondClampDofs;
    for (const auto& vertex : vertices(leafGridView))
        if (sheet->isSecondClamp(vertex.geometry().center()))
            secondClampDofs.push_back(deformationGridGeometry->vertexMapper().index(vertex));
    if (sheet->floatingClampGauge())
        couplingManager->setHole(*deformationGridGeometry, secondClampDofs);
    rotationGridVariables->init(x[rotationIdx]);
    deformationGridVariables->init(x[deformationIdx]);
    inPlaneGridVariables->init(x[inPlaneIdx]);

    VtkOutputModule<DeformationGridVariables, std::tuple_element_t<1, SolutionVector>>
        vtkWriter(*deformationGridVariables, x[deformationIdx], deformationProblem->name());
    vtkWriter.addVolumeVariable([](const auto& v){ return v.verticalDeformation(); }, "w");

    VtkOutputModule<InPlaneGridVariables, std::tuple_element_t<2, SolutionVector>>
        inPlaneVtkWriter(*inPlaneGridVariables, x[inPlaneIdx], deformationProblem->name() + "_inplane");
    inPlaneVtkWriter.addVolumeVariable([](const auto& v){ return v.displacement(0); }, "ux");
    inPlaneVtkWriter.addVolumeVariable([](const auto& v){ return v.displacement(1); }, "uy");
    std::array<std::vector<double>, 3> resultants{
        std::vector<double>(leafGridView.size(0)),
        std::vector<double>(leafGridView.size(0)),
        std::vector<double>(leafGridView.size(0))
    };
    inPlaneVtkWriter.addField(resultants[0], "Nxx", Vtk::FieldType::element);
    inPlaneVtkWriter.addField(resultants[1], "Nyy", Vtk::FieldType::element);
    inPlaneVtkWriter.addField(resultants[2], "Nxy", Vtk::FieldType::element);

    // The relaxation integrates c dw/dtau = -dE/dw implicitly in pseudo-time, which descends
    // the energy and settles on a stable equilibrium; a static Newton solve at these loads
    // seeks a stationary point of very high index and has no basin of attraction.
    auto xOld = x;
    const auto initialStep = getParam<Scalar>("Relaxation.InitialTimeStepSize", 1e-3);
    // a later load starts from a relaxed state, which tolerates a larger first step than the seed
    const auto continuationStep = getParam<Scalar>("Relaxation.ContinuationTimeStepSize", initialStep);
    const auto maxTimeStepSize = getParam<Scalar>("Relaxation.MaxTimeStepSize", 1e100);
    const auto maxSteps = getParam<int>("Relaxation.MaxSteps", 0);
    // with a duration every load step is relaxed over the same pseudo-time, so the damped sheet
    // is loaded at a constant rate and a step size pinned small by a transient does not end
    // the relaxation early; the step count then only caps it
    const auto loadStepDuration = getParam<Scalar>("Relaxation.LoadStepDuration", 0.0);
    // steady when the deflection moves by less than this fraction of the thickness in a
    // step taken at the largest pseudo-time step allowed
    const auto steadyStateTolerance = getParam<Scalar>("Relaxation.SteadyStateTolerance", 1e-2);
    const auto finalStaticSolve = getParam<bool>("Relaxation.FinalStaticSolve", false);
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, initialStep, 1e100);

    using Assembler = HoleGaugeAssembler<MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(rotationProblem, deformationProblem, inPlaneProblem),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry, inPlaneGridGeometry),
        std::make_tuple(rotationGridVariables, deformationGridVariables, inPlaneGridVariables),
        couplingManager, timeLoop, xOld
    );
    if (sheet->floatingClampGauge())
        assembler->setHole(secondClampDofs, FoepplVonKarmanPlateIndices::shearGradPotentialIdx,
                           FoepplVonKarmanPlateIndices::deformationEqIdx);

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;

    // the rest of the run, for either direct solver
    const auto run = [&](auto linearSolver) -> int
    {
    using LinearSolver = typename decltype(linearSolver)::element_type;
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    const auto length = sheet->length();
    const auto width = sheet->width();
    const auto thickness = sheet->thickness();
    const bool shear = sheet->loading() == SheetLoading::shear;

    // scan across the wrinkles: along the sheet for shear, across it for stretch
    const auto scanAxis = shear ? 0 : 1;
    const auto scanLength = shear ? length : width;
    const auto station = shear ? 0.5*width : 0.5*length;
    const auto cells = getParam<std::vector<int>>("Grid.Cells");
    const auto tolerance = 0.25*(shear ? width/cells[1] : length/cells[0]);
    // crests lower than this fraction of the central amplitude are not counted
    const auto crestThreshold = getParam<Scalar>("Problem.CrestThreshold", 0.25);

    using Indices = Dumux::FoepplVonKarmanPlateIndices;
    const auto scanProfile = [&]()
    {
        return lineProfile(
            *deformationGridGeometry, x[deformationIdx], scanAxis, station, tolerance,
            Indices::verticalDeformationIdx
        );
    };
    // the wrinkle amplitude over the central half of the scan line; the free edges curl
    // and would dominate an extreme taken over the whole sheet
    const auto centralAmplitude = [&](const auto& profile)
    { return profileAmplitude(profile, 0.25*scanLength, 0.75*scanLength); };
    const auto energy = [&]()
    {
        return strainEnergy(
            *rotationGridGeometry, *deformationGridGeometry, *inPlaneGridGeometry,
            x, *sheet, rotationIdx, deformationIdx, inPlaneIdx
        );
    };

    std::ofstream summary(deformationProblem->name() + "_wrinkles.dat");
    summary << "# load amplitude amplitudeOverThickness numCrests crestSpacing normalAngle wavelength"
               " predictedWavelength centralNxy flatNxy tensionFieldNxy relaxationSteps steadyState"
               " bendingEnergy membraneEnergy centralNxx centralNyy\n";
    summary << std::scientific << std::setprecision(8);

    const auto restore = [&]()
    {
        x = xOld;
        couplingManager->updateSolution(x);
        rotationGridVariables->update(x[rotationIdx]);
        deformationGridVariables->update(x[deformationIdx]);
        inPlaneGridVariables->update(x[inPlaneIdx]);
    };
    const auto accept = [&]()
    {
        xOld = x;
        rotationGridVariables->advanceTimeStep();
        deformationGridVariables->advanceTimeStep();
        inPlaneGridVariables->advanceTimeStep();
    };

    // A static solve with a pseudo-time step this large has no damping left in it. The
    // imperfection ladder, if given, coarsens the pattern by forcing a large deflection and
    // then lowering the forcing again while staying on the branch found.
    const auto imperfections = getParam<std::vector<Scalar>>(
        "Problem.ImperfectionSteps", std::vector<Scalar>{sheet->imperfection()}
    );
    timeLoop->setMaxTimeStepSize(1e100);
    timeLoop->setTimeStepSize(1e100);
    for (const auto imperfection : imperfections)
    {
        if (imperfection == 0.0 && imperfections.size() == 1)
            break;
        sheet->setImperfection(imperfection);
        std::cout << "\n### Branch search: imperfection " << imperfection << " ###\n";
        nonLinearSolver->solve(x);
        accept();
    }

    bool failed = false;
    Scalar previousLoad = loads.front();
    for (std::size_t step = 0; step < loads.size(); ++step)
    {
        const auto load = loads[step];
        sheet->setLoad(load);

        // the in-plane response beyond the locked-in strain is close to linear in the applied
        // load, so rescaling it gives a good initial guess and the right driven edge values;
        // from zero load there is nothing to rescale and the homogeneous state is added instead
        if (step > 0)
        {
            for (const auto& vertex : vertices(leafGridView))
            {
                const auto pos = vertex.geometry().center();
                auto& u = x[inPlaneIdx][inPlaneGridGeometry->vertexMapper().index(vertex)];
                const auto locked = sheet->preStrainDisplacement(pos);
                const auto affine = sheet->affineDisplacement(pos, load);
                for (int dir = 0; dir < GlobalPosition::dimension; ++dir)
                {
                    u[dir] -= locked[dir];
                    u[dir] = previousLoad != 0.0 ? u[dir]*(load/previousLoad) : u[dir] + affine[dir];
                    u[dir] += locked[dir];
                }
            }
            couplingManager->updateSolution(x);
            inPlaneGridVariables->update(x[inPlaneIdx]);
            accept();
        }
        previousLoad = load;

        std::cout << "\n### Applied load " << load << " ###\n";
        int relaxationSteps = 0;
        Scalar elapsed = 0.0;
        bool steady = maxSteps == 0;
        if (maxSteps > 0)
        {
            timeLoop->setMaxTimeStepSize(maxTimeStepSize);
            timeLoop->setTimeStepSize(step == 0 ? initialStep : continuationStep);
            {
                const auto [bending, membrane] = energy();
                std::cout << "Initial state: central A/t " << centralAmplitude(scanProfile())/thickness
                          << ", energy " << bending + membrane
                          << " (bending " << bending << ", membrane " << membrane << ")" << std::endl;
            }
            std::cout << std::setw(6) << "step" << std::setw(14) << "dtau"
                      << std::setw(14) << "max |w|/t" << std::setw(14) << "step change/t"
                      << std::setw(14) << "central A/t" << std::setw(16) << "energy" << std::endl;
            for (; relaxationSteps < maxSteps; ++relaxationSteps)
            {
                try {
                    nonLinearSolver->solve(x, *timeLoop);
                }
                catch (const Dune::Exception& e) {
                    timeLoop->setTimeStepSize(0.25*timeLoop->timeStepSize());
                    restore();
                    std::cout << "step " << relaxationSteps << " failed, retrying with dtau = "
                              << timeLoop->timeStepSize() << std::endl;
                    continue;
                }

                // the rate of the deflection alone: the potentials are of the order of the
                // bending moment and change by far more than the deflection does
                Scalar change = 0.0;
                for (std::size_t i = 0; i < x[deformationIdx].size(); ++i)
                    change = std::max(change, std::abs(
                        x[deformationIdx][i][Dumux::FoepplVonKarmanPlateIndices::verticalDeformationIdx]
                        - xOld[deformationIdx][i][Dumux::FoepplVonKarmanPlateIndices::verticalDeformationIdx]));
                accept();

                Scalar maxW = 0.0;
                for (const auto& dof : x[deformationIdx])
                    maxW = std::max(maxW, std::abs(dof[Dumux::FoepplVonKarmanPlateIndices::verticalDeformationIdx]));
                const auto [bending, membrane] = energy();
                std::cout << std::setw(6) << relaxationSteps << std::setw(14) << timeLoop->timeStepSize()
                          << std::setw(14) << maxW/thickness << std::setw(14) << change/thickness
                          << std::setw(14) << centralAmplitude(scanProfile())/thickness
                          << std::setw(16) << std::setprecision(10) << bending + membrane
                          << std::setprecision(6) << std::endl;

                elapsed += timeLoop->timeStepSize();
                if (change < steadyStateTolerance*thickness && timeLoop->timeStepSize() >= 0.5*maxTimeStepSize)
                {
                    steady = true;
                    ++relaxationSteps;
                    std::cout << "Steady state after " << relaxationSteps << " steps" << std::endl;
                    break;
                }
                if (loadStepDuration > 0.0 && elapsed >= (1.0 - 1e-10)*loadStepDuration)
                {
                    ++relaxationSteps;
                    break;
                }
                timeLoop->advanceTimeStep();
                auto next = nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize());
                if (loadStepDuration > 0.0)
                    next = std::min(next, loadStepDuration - elapsed);
                timeLoop->setTimeStepSize(next);
            }
            if (!steady)
                std::cout << "*** no steady state after " << relaxationSteps << " relaxation steps over pseudo-time "
                          << elapsed << " ***" << std::endl;
        }

        // the relaxed state is a stable equilibrium, where the static solve converges and
        // removes the last of the pseudo-time damping
        if (finalStaticSolve)
        {
            timeLoop->setMaxTimeStepSize(1e100);
            timeLoop->setTimeStepSize(1e100);
            try {
                nonLinearSolver->solve(x);
                accept();
            }
            catch (const Dune::Exception& e) {
                restore();
                std::cout << "static solve from the relaxed state failed, keeping the relaxed state" << std::endl;
                steady = false;
            }
        }

        rotationGridVariables->update(x[rotationIdx]);
        deformationGridVariables->update(x[deformationIdx]);
        inPlaneGridVariables->update(x[inPlaneIdx]);
        resultants = stressResultants(
            *inPlaneGridGeometry, *deformationGridGeometry, x, *sheet, inPlaneIdx, deformationIdx
        );
        vtkWriter.write(load);
        inPlaneVtkWriter.write(load);

        const auto profile = scanProfile();
        const auto amplitude = centralAmplitude(profile);
        const auto crests = crestPositions(profile, crestThreshold*amplitude);
        const auto spacing = meanCrestSpacing(crests, 0.25*scanLength, 0.75*scanLength);

        // the crest spacing along the scan line is the wavelength divided by the cosine of
        // the angle between the wrinkle normal and the scan direction
        const auto normalAngle = wrinkleNormalAngle(
            *deformationGridGeometry, x, 0.25*length, 0.75*length, deformationIdx
        );
        const auto normalToScan = shear ? normalAngle : normalAngle - 90.0;
        const auto wavelength = spacing*std::abs(std::cos(normalToScan*M_PI/180.0));

        const auto centralNxx = centralMean(*inPlaneGridGeometry, resultants[0], 0.25*length, 0.75*length);
        const auto centralNyy = centralMean(*inPlaneGridGeometry, resultants[1], 0.25*length, 0.75*length);
        const auto centralNxy = centralMean(*inPlaneGridGeometry, resultants[2], 0.25*length, 0.75*length);
        const auto flatNxy = sheet->flatShearResultant();
        const auto tensionFieldNxy = sheet->tensionFieldShearResultant();
        const auto predicted = sheet->predictedWavelength(load);
        const auto [bendingEnergy, membraneEnergy] = energy();
        std::cout << "Applied load: " << load << std::endl;
        std::cout << "Strain energy: " << bendingEnergy + membraneEnergy
                  << " (bending " << bendingEnergy << ", membrane " << membraneEnergy << ")" << std::endl;
        std::cout << "Shear resultant over the central half: " << centralNxy
                  << " (flat sheet " << flatNxy << ", tension field " << tensionFieldNxy << ")" << std::endl;
        std::cout << "Normal resultants over the central half: " << centralNxx << ", " << centralNyy
                  << " (locked-in strain times the membrane stiffness "
                  << sheet->membraneStiffness()*sheet->preStrain() << ")" << std::endl;
        std::cout << "Wrinkle amplitude: " << amplitude << " (A/t = " << amplitude/thickness << ")" << std::endl;
        std::cout << "Number of crests along the scan line: " << crests.size() << std::endl;
        std::cout << "Wrinkle normal at " << normalAngle << " degrees to the x axis" << std::endl;
        std::cout << "Mean crest spacing: " << spacing
                  << ", wavelength across the wrinkles: " << wavelength
                  << " (bending-tension balance " << predicted << ")" << std::endl;

        summary << load << " " << amplitude << " " << amplitude/thickness << " " << crests.size()
                << " " << spacing << " " << normalAngle << " " << wavelength << " " << predicted
                << " " << centralNxy << " " << flatNxy << " " << tensionFieldNxy << " "
                << relaxationSteps << " " << steady << " " << bendingEnergy << " " << membraneEnergy
                << " " << centralNxx << " " << centralNyy << "\n";
        summary.flush();

        std::ofstream profileFile(
            deformationProblem->name() + "_profile_" + std::to_string(step) + ".dat"
        );
        profileFile << "# applied load " << load << "\n# position w\n";
        profileFile << std::scientific << std::setprecision(8);
        for (const auto& [s, w] : profile)
            profileFile << s << " " << w << "\n";

        // the gates, each active only if its parameter is set
        const auto check = [&](const std::string& what, Scalar deviation, Scalar allowed)
        {
            if (allowed > 0.0 && deviation > allowed)
            {
                std::cerr << what << " deviates by " << deviation << ", more than the "
                          << allowed << " allowed" << std::endl;
                failed = true;
            }
        };
        if (!steady && getParam<bool>("Problem.RequireSteadyState", false))
        {
            std::cerr << "The relaxation did not reach a steady state" << std::endl;
            failed = true;
        }
        if (shear)
            check("The central shear resultant from the tension-field value",
                  std::abs(centralNxy/tensionFieldNxy - 1.0),
                  getParam<Scalar>("Problem.ShearResultantTolerance", -1.0));
        const auto expectedNormalAngle = getParam<Scalar>("Problem.NormalAngle", sheet->tensionFieldNormalAngle());
        check("The wrinkle normal angle from " + std::to_string(expectedNormalAngle) + " degrees",
              std::abs(normalAngle - expectedNormalAngle),
              getParam<Scalar>("Problem.NormalAngleTolerance", -1.0));
        // the sheet tensioned across its width carries N_yy = C eps0 away from the free ends,
        // which the wrinkles of a sheared sheet relieve: this gates the shear-free state
        if (const auto preStrainResultant = sheet->membraneStiffness()*sheet->preStrain(); preStrainResultant > 0.0)
            check("The resultant across the sheet from the locked-in strain",
                  std::abs(centralNyy/preStrainResultant - 1.0),
                  getParam<Scalar>("Problem.PreStrainResultantTolerance", -1.0));
        check("The wavelength from the bending-tension balance", std::abs(wavelength/predicted - 1.0),
              getParam<Scalar>("Problem.WavelengthTolerance", -1.0));
        const auto minAmplitude = getParam<Scalar>("Problem.MinAmplitudeOverThickness", -1.0);
        if (minAmplitude > 0.0 && amplitude/thickness < minAmplitude)
        {
            std::cerr << "The wrinkle amplitude " << amplitude/thickness
                      << " thicknesses is below the " << minAmplitude << " required" << std::endl;
            failed = true;
        }
    }

    return failed ? 1 : 0;
    };

    const auto solverType = getParam<std::string>("LinearSolver.Type", "umfpack");
    if (solverType == "umfpack")
        return run(std::make_shared<UMFPackIstlSolver<SeqLinearSolverTraits, LAT>>());
#if DUMUX_HAVE_MUMPS
    if (solverType == "mumps")
    {
        // built from the grid geometries so that a distributed grid gives the parallel solver
        using LinearSolver = DirectSolverMumps<LinearSolverTraits<DeformationGridGeometry>, LAT>;
        return run(std::make_shared<LinearSolver>(
            std::make_tuple(rotationGridGeometry, deformationGridGeometry, inPlaneGridGeometry)
        ));
    }
#endif
    DUNE_THROW(Dune::IOError, "Unknown LinearSolver.Type '" << solverType << "', expected umfpack or mumps");
}
