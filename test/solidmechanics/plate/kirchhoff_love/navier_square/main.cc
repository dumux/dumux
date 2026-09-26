// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief Navier's double sine series of the simply supported square plate under a uniform load
 *
 * The load \f$ F \f$ acts against \f$ w \f$, so
 * \f$ w = -\frac{16FL^4}{\pi^6D}\sum_{m,n\,\mathrm{odd}}\frac{\sin(m\pi x/L)\sin(n\pi y/L)}{mn(m^2+n^2)^2} \f$.
 */
class NavierSquare
{
public:
    NavierSquare(double F, double D, double nu, double L, int maxIndex)
    : F_(F), D_(D), nu_(nu), L_(L), maxIndex_(maxIndex) {}

    double deflection(double x, double y, int maxIndex = 201) const
    {
        std::vector<double> sx, sy;
        for (int m = 1; m <= maxIndex; m += 2)
        {
            sx.push_back(std::sin(m*M_PI*x/L_));
            sy.push_back(std::sin(m*M_PI*y/L_));
        }
        double sum = 0.0;
        for (int i = 0, m = 1; m <= maxIndex; ++i, m += 2)
            for (int j = 0, n = 1; n <= maxIndex; ++j, n += 2)
                sum += sx[i]*sy[j]/(m*n*square(m*m + n*n));
        return -16.0*F_*L_*L_*L_*L_/(std::pow(M_PI, 6)*D_)*sum;
    }

    //! The moment tensor \f$ \mathbf{M} = -D\{(1-\nu)\nabla\nabla w + \nu\Delta w\mathbf{I}\} \f$
    Dune::FieldMatrix<double, 2, 2> moment(double x, double y) const
    {
        double wxx = 0.0, wyy = 0.0, wxy = 0.0;
        for (int m = 1; m <= maxIndex_; m += 2)
            for (int n = 1; n <= maxIndex_; n += 2)
            {
                const auto sx = std::sin(m*M_PI*x/L_), sy = std::sin(n*M_PI*y/L_);
                const auto cx = std::cos(m*M_PI*x/L_), cy = std::cos(n*M_PI*y/L_);
                const auto denominator = square(m*m + n*n);
                wxx += double(m)/n*sx*sy/denominator;
                wyy += double(n)/m*sx*sy/denominator;
                wxy -= cx*cy/denominator;
            }
        const auto scale = 16.0*F_*L_*L_/(std::pow(M_PI, 4)*D_);
        wxx *= scale; wyy *= scale; wxy *= scale;
        Dune::FieldMatrix<double, 2, 2> M(0.0);
        M[0][0] = -D_*(wxx + nu_*wyy);
        M[1][1] = -D_*(wyy + nu_*wxx);
        M[0][1] = M[1][0] = -D_*(1.0 - nu_)*wxy;
        return M;
    }

private:
    static double square(double v) { return v*v; }
    double F_, D_, nu_, L_;
    int maxIndex_;
};

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using RotationTypeTag = Properties::TTag::KLNavierSquareTestRotation;
    using DeformationTypeTag = Properties::TTag::KLNavierSquareTestDeformation;
    using CommonTypeTag = Properties::TTag::KLNavierSquareTestCommon;

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

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);
    nonLinearSolver->solve(x);

    if (getParam<bool>("Vtk.Write", false))
    {
        VtkOutputModule<DeformationGridVariables, std::tuple_element_t<1, SolutionVector>>
            vtkWriter(*deformationGridVariables, x[deformationIdx], deformationProblem->name());
        vtkWriter.addVolumeVariable([](const auto& v){ return v.verticalDeformation(); }, "w");
        vtkWriter.addVolumeVariable([](const auto& v){ return v.shearCurlPotential(); }, "psi");
        vtkWriter.addVolumeVariable([](const auto& v){ return v.shearGradPotential(); }, "phi");
        vtkWriter.write(1.0);
    }

    const auto nu = getParam<double>("Problem.PoissonRatio");
    const auto E = getParam<double>("Problem.E");
    const auto t = getParam<double>("Problem.Thickness");
    const auto D = E*t*t*t/(12.0*(1.0 - nu*nu));
    const auto F = getParam<double>("Problem.Force");
    const auto L = getParam<double>("Problem.Length", 1.0);
    const NavierSquare navier(F, D, nu, L, getParam<int>("Problem.SeriesTerms", 2001));
    using Position = Dune::FieldVector<double, 2>;

    // The gradient of the rotations jumps across elements, so the moment at a vertex is the
    // average over the elements that touch it.
    const auto momentAt = [&](const Position& point)
    {
        Dune::FieldMatrix<double, 2, 2> M(0.0);
        int count = 0;
        for (const auto& element : elements(rotationGridGeometry->gridView()))
        {
            const auto geo = element.geometry();
            bool touches = false;
            for (int c = 0; c < geo.corners(); ++c)
                if ((geo.corner(c) - point).two_norm() < 1e-8*L)
                    touches = true;
            if (!touches)
                continue;
            const auto elemSol = elementSolution(element, x[rotationIdx], *rotationGridGeometry);
            const auto g = evalGradients(element, geo, *rotationGridGeometry, elemSol, point);
            M[0][0] += -D*(g[0][0] + nu*g[1][1]);
            M[1][1] += -D*(g[1][1] + nu*g[0][0]);
            M[0][1] += -D*(1.0 - nu)*0.5*(g[0][1] + g[1][0]);
            ++count;
        }
        if (count == 0)
            DUNE_THROW(Dune::InvalidStateException, "Moment evaluation requires a grid vertex at " << point);
        M /= count;
        M[1][0] = M[0][1];
        return M;
    };

    const Position centre(0.5*L);
    const auto wCentreRef = navier.deflection(centre[0], centre[1]);
    double wCentre = 0.0, maxError = 0.0;
    for (const auto& vertex : vertices(leafGridView))
    {
        const auto pos = vertex.geometry().center();
        const auto w = x[deformationIdx][deformationGridGeometry->vertexMapper().index(vertex)][1];
        if ((pos - centre).two_norm() < 1e-8*L)
            wCentre = w;
        maxError = std::max(maxError, std::abs(w - navier.deflection(pos[0], pos[1])));
    }

    const auto MCentre = momentAt(centre);
    const auto MCentreRef = navier.moment(centre[0], centre[1]);
    const auto MCornerRef = navier.moment(0.0, 0.0);

    std::cout << std::setprecision(10)
              << "Mapping: " << getParam<std::string>("Problem.Mapping") << "\n"
              << "Cells: " << getParam<std::vector<int>>("Grid.Cells")[0] << "\n"
              << "Centre deflection: " << wCentre << " reference " << wCentreRef << "\n"
              << "Maximum deflection error: " << maxError/std::abs(wCentreRef) << "\n"
              << "Centre moment M11: " << MCentre[0][0] << " reference " << MCentreRef[0][0] << "\n";

    // At a right-angled corner the Kirchhoff corner force is twice the twisting moment.
    const std::vector<Position> corners{{0.0, 0.0}, {L, 0.0}, {L, L}, {0.0, L}};
    for (const auto& corner : corners)
    {
        const auto M = momentAt(corner);
        const auto MRef = navier.moment(corner[0], corner[1]);
        std::cout << "Corner (" << corner[0] << "," << corner[1] << ") M12: " << M[0][1]
                  << " reference " << MRef[0][1] << "\n";
    }
    std::cout << "Corner force 2|M12|/(F L^2): " << 2.0*std::abs(MCornerRef[0][1])/(F*L*L) << std::endl;

    return 0;
}
