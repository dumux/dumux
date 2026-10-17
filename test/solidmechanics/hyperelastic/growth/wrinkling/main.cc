// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Onset of growth-induced wrinkling: the growth stretch of the film is increased along the
 *        flat solution until the smallest real eigenvalue of the Jacobian changes sign, and the
 *        critical stretch is then located by the Illinois method
 *
 * If `Problem.BranchAmplitudes` is given, the wrinkled branch is traced beyond the onset: for each
 * amplitude a, the top corner at x = 0 is displaced by a from the flat state, and the growth stretch
 * at which its reaction vanishes is found by the secant method, starting from the previous state
 * moved along the eigenvector of the onset. The stretches are fitted by
 * \f$ \vartheta = \vartheta_0 + c_2 a^2 + c_4 a^4 \f$, whose intercept has to agree with the
 * critical stretch from the eigenvalues.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/eigensolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonlinearSolver(assembler, linearSolver);

    using EigenSolver = ShiftInvertEigenSolver<typename Assembler::JacobianMatrix, SolutionVector>;
    const EigenSolver eigenSolver;
    const int numEigenvalues = getParam<int>("EigenSolver.NumEigenvalues");

    // the smallest real part of the eigenvalues closest to zero of the Jacobian at the flat solution
    // for the given growth stretch, and its eigenvector
    const auto smallestEigenvalue = [&](Scalar stretch)
    {
        problem->setGrowthStretch(stretch);
        nonlinearSolver.solve(x);
        assembler->assembleJacobianAndResidual(x);
        const auto result = eigenSolver.solve(assembler->jacobian(), x, numEigenvalues);
        std::size_t k = 0;
        for (std::size_t i = 1; i < result.eigenvalues.size(); ++i)
            if (result.eigenvalues[i].real() < result.eigenvalues[k].real())
                k = i;
        return std::make_pair(result.eigenvalues[k].real(), result.eigenvectorsReal[k]);
    };

    // bracket the sign change
    const auto step = getParam<Scalar>("Problem.GrowthStep");
    Scalar stretchLow = 1.0, stretchHigh = 1.0;
    auto [eigLow, modeLow] = smallestEigenvalue(stretchLow);
    auto [eigHigh, mode] = std::make_pair(eigLow, modeLow);
    while (eigHigh > 0.0)
    {
        stretchLow = stretchHigh;
        eigLow = eigHigh;
        stretchHigh += step;
        std::tie(eigHigh, mode) = smallestEigenvalue(stretchHigh);
        std::cout << "growth stretch " << stretchHigh << ": smallest eigenvalue " << eigHigh << std::endl;
    }

    // Illinois method
    const auto tolerance = getParam<Scalar>("Problem.StretchTolerance");
    int side = 0;
    Scalar stretch = stretchHigh;
    while (stretchHigh - stretchLow > tolerance)
    {
        stretch = (stretchLow*eigHigh - stretchHigh*eigLow)/(eigHigh - eigLow);
        Scalar eig;
        std::tie(eig, mode) = smallestEigenvalue(stretch);
        if (eig > 0.0)
        {
            stretchLow = stretch; eigLow = eig;
            if (side == -1) eigHigh *= 0.5;
            side = -1;
        }
        else
        {
            stretchHigh = stretch; eigHigh = eig;
            if (side == 1) eigLow *= 0.5;
            side = 1;
        }
        if (eig == 0.0)
            break;
    }

    // number of half wavelengths of the mode along the top surface
    std::vector<std::pair<Scalar, Scalar>> top;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto pos = ipData(elemDisc, localDof).global();
            if (pos[1] > gridGeometry->bBoxMax()[1] - 1e-8)
                top.emplace_back(pos[0], mode[localDof.dofIndex()][1]);
        }
    }
    std::sort(top.begin(), top.end());
    top.erase(std::unique(top.begin(), top.end(), [](auto a, auto b) { return std::abs(a.first - b.first) < 1e-12; }), top.end());
    Scalar maxAmplitude = 0.0;
    for (const auto& p : top)
        maxAmplitude = std::max(maxAmplitude, std::abs(p.second));
    int signChanges = 0;
    Scalar lastSign = 0.0;
    for (const auto& p : top)
    {
        if (std::abs(p.second) < 1e-6*maxAmplitude)
            continue;
        const Scalar sign = p.second > 0.0 ? 1.0 : -1.0;
        if (lastSign != 0.0 && sign != lastSign)
            ++signChanges;
        lastSign = sign;
    }

    Scalar hMax = 0.0;
    for (const auto& element : elements(leafGridView))
    {
        const auto& geometry = element.geometry();
        hMax = std::max(hMax, (geometry.corner(geometry.corners() - 1) - geometry.corner(0)).two_norm());
    }

    std::cout << std::setprecision(12)
              << "Max element diameter: " << hMax << "\n"
              << "Critical growth stretch: " << stretch << "\n"
              << "Wrinkle half-wavelengths: " << signChanges << "\n";

    const auto writeVtk = [&]
    {
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
        {
            if (getParam<bool>("Problem.EnableVtkOutput", false))
            {
                VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
                vtkWriter.addField(x, "d");
                vtkWriter.addField(mode, "mode");
                vtkWriter.write(1.0);
            }
        }
    };

    if (!hasParam("Problem.BranchAmplitudes"))
    {
        writeVtk();
        return 0;
    }

    const auto criticalStretch = stretch;
    const auto amplitudes = getParam<std::vector<Scalar>>("Problem.BranchAmplitudes");
    const auto pointDof = problem->pointDofIndex();
    const auto reaction = [&](Scalar stretch, Scalar amplitude)
    {
        problem->setGrowthStretch(stretch);
        problem->setPointDisplacement((problem->filmThicknessStretch(stretch) - 1.0)*problem->filmThickness() + amplitude);
        problem->enablePointConstraint(true);
        nonlinearSolver.solve(x);
        problem->enablePointConstraint(false);
        assembler->assembleResidual(x);
        return assembler->residual()[pointDof][1];
    };

    std::vector<Scalar> branchStretches;
    for (std::size_t i = 0; i < amplitudes.size(); ++i)
    {
        // the initial guess moves along the eigenvector of the onset by the increment of the amplitude
        x.axpy((amplitudes[i] - (i > 0 ? amplitudes[i-1] : 0.0))/mode[pointDof][1], mode);

        // first guesses from the parabolic shape of the branch
        Scalar s0 = criticalStretch + 1e-4;
        if (i > 0)
        {
            const auto ratio = amplitudes[i]/amplitudes[i-1];
            s0 = criticalStretch + (branchStretches[i-1] - criticalStretch)*ratio*ratio;
        }
        Scalar s1 = s0 + 1e-4*std::max(std::abs(s0 - criticalStretch), Scalar(1e-3));
        Scalar r0 = reaction(s0, amplitudes[i]);
        Scalar r1 = reaction(s1, amplitudes[i]);
        for (int iter = 0; iter < 50 && std::abs(s1 - s0) > 1e-13; ++iter)
        {
            const Scalar s2 = s1 - r1*(s1 - s0)/(r1 - r0);
            s0 = s1; r0 = r1;
            s1 = s2; r1 = reaction(s1, amplitudes[i]);
        }
        branchStretches.push_back(s1);
        std::cout << "Branch amplitude " << amplitudes[i] << ": growth stretch " << s1 << "\n";
    }

    // least-squares fit of theta_0 + c_2 a^2 + c_4 a^4
    Dune::FieldMatrix<Scalar, 3, 3> normal(0.0);
    Dune::FieldVector<Scalar, 3> rhs(0.0), coefficients;
    for (std::size_t i = 0; i < amplitudes.size(); ++i)
    {
        const Scalar a2 = amplitudes[i]*amplitudes[i];
        const Dune::FieldVector<Scalar, 3> basis{1.0, a2, a2*a2};
        for (int r = 0; r < 3; ++r)
        {
            rhs[r] += basis[r]*branchStretches[i];
            for (int c = 0; c < 3; ++c)
                normal[r][c] += basis[r]*basis[c];
        }
    }
    normal.solve(coefficients, rhs);
    std::cout << "Branch intercept: " << coefficients[0] << "\n"
              << "Branch intercept deviation: " << std::abs(coefficients[0] - criticalStretch) << "\n"
              << "Branch coefficient c2: " << coefficients[1] << "\n"
              << "Branch coefficient c4: " << coefficients[2] << "\n";

    writeVtk();
    return 0;
}
