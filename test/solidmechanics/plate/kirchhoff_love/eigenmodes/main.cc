// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/discretization/localview.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/io/format.hh>

#include "eigenmodes.hh"
#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using RotationTypeTag = Properties::TTag::KLPlateEigenmodesRotation;
    using DeformationTypeTag = Properties::TTag::KLPlateEigenmodesDeformation;
    using CommonTypeTag = Properties::TTag::KLPlateEigenmodesCommon;

    using Grid = GetPropType<CommonTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using RotationGridGeometry = GetPropType<RotationTypeTag, Properties::GridGeometry>;
    using DeformationGridGeometry = GetPropType<DeformationTypeTag, Properties::GridGeometry>;
    auto rotationGridGeometry = std::make_shared<RotationGridGeometry>(leafGridView);
    auto deformationGridGeometry = std::make_shared<DeformationGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<RotationTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(rotationGridGeometry, deformationGridGeometry);

    using RotationProblem = GetPropType<RotationTypeTag, Properties::Problem>;
    auto rotationProblem = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManager);

    using DeformationProblem = GetPropType<DeformationTypeTag, Properties::Problem>;
    auto deformationProblem = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManager);

    constexpr auto rotationIdx = CouplingManager::rotationIdx;
    constexpr auto deformationIdx = CouplingManager::deformationIdx;
    using Traits = MultiDomainTraits<RotationTypeTag, DeformationTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    rotationProblem->applyInitialSolution(x[rotationIdx]);
    deformationProblem->applyInitialSolution(x[deformationIdx]);

    using RotationGridVariables = GetPropType<RotationTypeTag, Properties::GridVariables>;
    auto rotationGridVariables = std::make_shared<RotationGridVariables>(rotationProblem, rotationGridGeometry);

    using DeformationGridVariables = GetPropType<DeformationTypeTag, Properties::GridVariables>;
    auto deformationGridVariables = std::make_shared<DeformationGridVariables>(deformationProblem, deformationGridGeometry);

    couplingManager->init(rotationProblem, deformationProblem, x);
    rotationGridVariables->init(x[rotationIdx]);
    deformationGridVariables->init(x[deformationIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(rotationProblem, deformationProblem),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry),
        std::make_tuple(rotationGridVariables, deformationGridVariables),
        couplingManager
    );

    // build lumped mass vector: only the vertical deformation (w) DOF carries inertia
    // M_j = linearMassDensity * controlVolumeArea_j for the w-component, zero elsewhere
    using DeformationIndices = typename GetPropType<DeformationTypeTag, Properties::ModelTraits>::Indices;
    const auto linearMassDensity = getParam<double>("Problem.LinearMassDensity");
    SolutionVector massVec;
    rotationProblem->applyInitialSolution(massVec[rotationIdx]);
    deformationProblem->applyInitialSolution(massVec[deformationIdx]);
    massVec = 0.0;
    {
        auto fvGeometry = localView(*deformationGridGeometry);
        for (const auto& element : elements(deformationGridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                massVec[deformationIdx][scv.dofIndex()][DeformationIndices::verticalDeformationIdx]
                    += linearMassDensity * scv.volume();
        }
    }

    // assemble the Jacobian (stiffness matrix) at zero displacement
    assembler->assembleJacobianAndResidual(x);
    Dune::MatrixAdapter<
        std::decay_t<decltype(assembler->jacobian())>,
        SolutionVector, SolutionVector
    > op(assembler->jacobian());

    const int numEigenmodes = getParam<int>("EigenSolver.NumberOfEigenmodes", 3);
    const auto [eigenvalues, eigenvectors] = Dumux::smallestNEigenmodes(op, assembler->residual(), massVec, numEigenmodes);

    // non-dimensional frequency: Omega = omega * a^2 * sqrt(rho*h/D)
    // Reference values for CCCC square plate from Leissa (1973) "Vibration of Plates", Table 4.3
    // Ω = ω·a²·√(ρh/D), modes ordered by increasing frequency
    static constexpr std::array<double, 12> leissa1973CCCC = {
        35.985,  // (1,1)
        73.413,  // (1,2)
        73.413,  // (2,1)
       108.27,   // (2,2)
       131.64,   // (1,3)
       131.64,   // (3,1)
       165.15,   // (2,3)
       165.15,   // (3,2)
       210.53,   // (1,4)
       210.53,   // (4,1)
       220.06,   // (3,3)
       242.66,   // (2,4)
    };

    std::cout << "Non-dimensional frequencies Omega = omega*a^2*sqrt(rho*h/D):\n";
    std::cout << Fmt::format("  {:>6}  {:>12}  {:>12}  {:>10}\n", "mode", "computed", "Leissa(1973)", "rel.error");
    const auto a = getParam<std::vector<double>>("Grid.UpperRight")[0];
    const auto D = rotationProblem->D({0.0, 0.0});
    const auto nonDimFactor = a*a*std::sqrt(linearMassDensity/D);
    for (int i = 0; i < numEigenmodes; ++i)
    {
        const double omega = std::sqrt(std::real(eigenvalues[i]))*nonDimFactor;
        if (i < static_cast<int>(leissa1973CCCC.size()))
        {
            const double ref = leissa1973CCCC[i];
            std::cout << Fmt::format("  {:>6}  {:>12.4f}  {:>12.3f}  {:>10.4f}\n",
                i, omega, ref, (omega - ref)/ref);
        }
        else
            std::cout << Fmt::format("  {:>6}  {:>12.4f}\n", i, omega);
    }

    // write eigenmodes to VTK
    VtkOutputModule<DeformationGridVariables, std::tuple_element_t<1, SolutionVector>>
        vtkWriter(*deformationGridVariables, x[deformationIdx], deformationProblem->name());

    std::vector<std::vector<double>> hmodes(numEigenmodes);
    for (int i = 0; i < numEigenmodes; ++i)
    {
        hmodes[i].resize(deformationGridGeometry->numDofs());
        for (std::size_t j = 0; j < deformationGridGeometry->numDofs(); ++j)
            hmodes[i][j] = eigenvectors[i][deformationIdx][j][DeformationIndices::verticalDeformationIdx];
        vtkWriter.addField(hmodes[i], Fmt::format("h_mode_{:0>3d}", i));
    }
    vtkWriter.write(0.0);

    return 0;
}
