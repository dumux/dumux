// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Thick cylindrical shell under pressure — Reese, Wriggers, Reddy (2000),
// Comput. Struct. 75, 291-304, §4.1 ("cylindric shell"); also reproduced as the
// F-bar locking benchmark in Elguedj et al. (2008), §5.4. Solved here with PQ2 (T2)
// single-field hyperelasticity and a compressible neo-Hookean material model.
//
// The load is applied incrementally (load stepping) to aid Newton convergence in
// this large-deformation regime, mirroring the OGS time ramp.
//
#include <config.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/assembler.hh>

#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/gridwriter.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::ThickCylindricalShell;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    Dumux::initialize(argc, argv);
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
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    static constexpr int dim = GridGeometry::GridView::dimension;
    using DisplacementVector = Dune::FieldVector<Scalar, dim>;

    // IO::GridWriter writes vertex DOFs as P1 point field (edge-midpoint DOFs are
    // skipped in the visual output but are still used in the computation).
    IO::GridWriter writer{
        IO::Format::pvd_with(IO::Format::vtu.with({
            .encoder = IO::Encoding::ascii,
            .compressor = IO::Compression::none,
            .data_format = IO::VTK::DataFormat::inlined
        })),
        leafGridView,
        getParam<std::string>("Problem.Name"),
        IO::order<1>
    };
    writer.setPointField("d", [&](const auto& vertex) -> DisplacementVector {
        return x[gridGeometry->dofMapper().index(vertex)];
    }, IO::Precision::float64);
    writer.write(0.0);

    // ---------------------------------------------------------------------
    // Crisfield cylindrical arc-length (Riks) path-following.
    //
    // The shell snaps through a limit point at ~37.5% load, so load control
    // cannot reach full load. Here the load factor lambda is an unknown,
    // constrained by a fixed arc length |Delta x| = dl, so the solver can trace
    // *around* the limit point. We continue until the path returns to lambda=1.
    //
    // Residual convention (DuMux): J*delta = r, x <- x - delta. With a dead load
    //   r(x,lambda) = r_int(x) - lambda*Q,   Q = external load vector (unit factor)
    // the correction is  delta_x = -J^{-1} r + dLam * J^{-1} Q.
    // ---------------------------------------------------------------------
    using Vec = SolutionVector;
    assembler->setLinearSystem();

    // External load vector at unit factor: Q = r(lambda=0) - r(lambda=1) at x=0.
    problem->setLoadFactor(0.0); assembler->assembleResidual(x); Vec Q = assembler->residual();
    problem->setLoadFactor(1.0); assembler->assembleResidual(x); { Vec r1 = assembler->residual(); Q -= r1; }
    const Scalar normQ = Q.two_norm();
    const Scalar resTol = 1e-7*normQ;
    std::cout << "[arc] ||Q|| = " << normQ << ", resTol = " << resTol << std::endl;

    const auto numLoadSteps = getParam<int>("Problem.LoadSteps", 40);
    const int maxArcSteps   = getParam<int>("Problem.MaxArcSteps", 2000);
    const int maxIter       = getParam<int>("Newton.MaxSteps", 30);
    const int targetIter    = 4;   // desired correctors per step (arc-length adaption)

    Scalar lambda = 0.0;
    Vec x0 = x;                    // last converged solution (zero)
    Vec dXprev = x; dXprev = 0.0;  // previous converged increment (predictor sign)

    // initial arc length ~ first nominal load step's tangential displacement
    problem->setLoadFactor(0.0);
    assembler->assembleJacobianAndResidual(x0);
    Vec dut = x; dut = 0.0;
    const auto solRes = linearSolver->solve(assembler->jacobian(), dut, Q);
    { Scalar m = 0.0; std::size_t mi = 0;
      for (std::size_t i = 0; i < dut.size(); ++i) if (dut[i].two_norm() > m) { m = dut[i].two_norm(); mi = i; }
      std::cout << "[arc] initial tangent: solve converged=" << solRes.converged
                << ", ||dut||=" << dut.two_norm() << ", max|dut|=" << m*1e3 << " mm at dof " << mi
                << " u=[" << dut[mi][0]*1e3 << "," << dut[mi][1]*1e3 << "," << dut[mi][2]*1e3 << "]" << std::endl; }
    Scalar dl = (1.0/numLoadSteps)*dut.two_norm();
    const Scalar dlMin = dl*1e-4;

    int outStep = 0;
    bool reachedFull = false;
    for (int step = 0; step < maxArcSteps && !reachedFull; ++step)
    {
        // ---- predictor: tangential step of length dl ----
        problem->setLoadFactor(lambda);
        assembler->assembleJacobianAndResidual(x0);
        dut = 0.0;
        linearSolver->solve(assembler->jacobian(), dut, Q);
        Scalar sign = (step > 0 && (dut*dXprev) < 0.0) ? -1.0 : 1.0;
        Scalar dLamTot = sign*dl/dut.two_norm();
        Vec dX = dut; dX *= dLamTot;                 // total increment this step
        Scalar lam = lambda + dLamTot;
        x = x0; x += dX;
        if (step == 0) std::cout << "[pred] dLamTot=" << dLamTot << " lam=" << lam << " |dX|=" << dX.two_norm() << std::endl;

        // ---- corrector: Newton with the cylindrical constraint |dX| = dl ----
        bool converged = false;
        int it = 0;
        for (; it < maxIter; ++it)
        {
            problem->setLoadFactor(lam);
            assembler->assembleJacobianAndResidual(x);
            Vec r = assembler->residual();
            const Scalar rnorm = r.two_norm();
            using std::isnan;
            if (isnan(rnorm)) break;                 // element inverted -> cut dl, retry
            if (rnorm < resTol) { converged = true; break; }

            Vec dur = x; dur = 0.0; linearSolver->solve(assembler->jacobian(), dur, r);
            dut = 0.0;              linearSolver->solve(assembler->jacobian(), dut, Q);

            Vec v = dX; v -= dur;                    // v = dX - J^{-1} r
            const Scalar a = dut*dut;
            const Scalar b = 2.0*(v*dut);
            const Scalar c = (v*v) - dl*dl;
            const Scalar disc = b*b - 4.0*a*c;
            if (a <= 0.0 || disc < 0.0) break;       // arc misses Newton line -> cut dl, retry

            const Scalar sq = std::sqrt(disc);
            const Scalar root1 = (-b + sq)/(2.0*a);
            const Scalar root2 = (-b - sq)/(2.0*a);
            Vec t1 = v; t1.axpy(root1, dut);         // candidate new total increments
            Vec t2 = v; t2.axpy(root2, dut);
            const Scalar dLam = ((t1*dX) >= (t2*dX)) ? root1 : root2;  // keep aligned with dX
            if (step == 0) std::cout << "  [corr it" << it << "] rnorm=" << rnorm
                << " |dur|=" << dur.two_norm() << " roots=" << root1 << "," << root2
                << " t1.dX=" << (t1*dX) << " t2.dX=" << (t2*dX) << " -> dLam=" << dLam
                << " newlam=" << (lam+dLam) << std::endl;

            dX = v; dX.axpy(dLam, dut);              // dX = v + dLam*dut
            lam += dLam;
            x = x0; x += dX;
        }

        if (!converged)
        {
            x = x0;                                  // roll back to last converged state
            dl *= 0.5;
            if (dl < dlMin) { std::cout << "  --> dl below threshold; STOPPING at lambda=" << lambda << "\n"; break; }
            continue;
        }

        // accept
        x0 = x; dXprev = dX; lambda = lam;
        writer.write(double(++outStep));

        Scalar stepMax = 0.0; std::size_t stepMaxDof = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i].two_norm() > stepMax) { stepMax = x[i].two_norm(); stepMaxDof = i; }
        std::cout << "  --> step " << outStep << ": lambda = " << lambda
                  << ", max |u| = " << stepMax*1e3 << " mm (corr its = " << it << ", dl = " << dl << ")"
                  << std::endl;

        if (lambda >= 1.0) reachedFull = true;
        else
        {
            // adapt the arc length toward the target corrector count
            dl *= std::sqrt(Scalar(targetIter)/Scalar(std::max(1, it)));
            dl = std::min(dl, (1.0/numLoadSteps)*dut.two_norm()*4.0);
        }
    }

    std::cout << "\nArc-length finished: lambda = " << lambda
              << (reachedFull ? " (full load reached)\n" : " (did NOT reach full load)\n");
    problem->setLoadFactor(lambda);

    // Report the maximum displacement magnitude and its location.
    Scalar maxNorm = 0.0;
    std::size_t maxDof = 0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const Scalar norm = x[i].two_norm();
        if (norm > maxNorm) { maxNorm = norm; maxDof = i; }
    }

    // Point A (Reese Fig. 1): inner-radius node on the loaded top crown at the
    // z = 0 symmetry plane, i.e. closest to (0, +r_i, 0).
    const Dune::FieldVector<Scalar, dim> pointA{0.0, 0.008, 0.0};
    Scalar bestDist = std::numeric_limits<Scalar>::max();
    std::size_t dofA = 0;
    for (const auto& vertex : vertices(leafGridView))
    {
        const auto idx = gridGeometry->dofMapper().index(vertex);
        const auto d = (vertex.geometry().corner(0) - pointA).two_norm();
        if (d < bestDist) { bestDist = d; dofA = idx; }
    }

    std::cout << "\n"
              << "Thick cylindrical shell — PQ2 (T2), compressible neo-Hookean\n"
              << std::string(70, '-') << "\n"
              << std::scientific << std::setprecision(8)
              << "max |u|        = " << maxNorm << " m\n"
              << "u at max-dof   = [" << x[maxDof][0] << ", " << x[maxDof][1] << ", " << x[maxDof][2] << "] m\n"
              << "u at point A   = [" << x[dofA][0] << ", " << x[dofA][1] << ", " << x[dofA][2] << "] m\n"
              << "  |vertical (u_y) at A| = " << std::abs(x[dofA][1])*1e3 << " mm\n"
              << std::string(70, '-') << "\n"
              << "Reference: Reese et al. (2000) Fig. 2(a) / Elguedj et al. (2008) Fig. 23,\n"
              << "  vertical displacement at point A converges to ~16.5 mm.\n"
              << "PQ2 is locking-free; result should be in that range.\n";

    return 0;
}
