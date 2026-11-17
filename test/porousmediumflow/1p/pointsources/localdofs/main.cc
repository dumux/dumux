// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Point sources for control-volume, hybrid and finite element schemes.
 *
 * Three checks with the single-phase model, Darcy flow with unit coefficients and zero pressure on the boundary:
 * - Residual entries: At a uniform pressure the residual of an element consists of the point sources only.
 *   A point source \f$ q \f$ at \f$ \mathbf{x} \f$ contributes \f$ -q \f$ to the control volume that contains \f$ \mathbf{x} \f$
 *   and \f$ -q \varphi_i(\mathbf{x}) \f$ to every local dof \f$ i \f$ without control volume. The Lagrange basis functions
 *   \f$ \varphi_i \f$ are evaluated independently of the discretization. A point on a face between two elements is
 *   split equally between them.
 * - Global balance: The rows of the control volumes conserve, so the reaction at the boundary
 *   adds up to the strength of the point source. The same sum is reported for finite element schemes.
 * - Green's function: In one dimension the Galerkin solution of \f$ -p'' = q \delta(x - x_0) \f$ is exact at the vertices.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/concepts.hh>

#include <dumux/io/grid/gridmanager_oned.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#if HAVE_DUNE_UGGRID
#include <dumux/io/grid/gridmanager_ug.hh>
#endif

#include "properties.hh"

#ifndef ORDER
#error "ORDER, the polynomial order of the discretization, has to be defined"
#endif

namespace Dumux::Test {

struct Checker
{
    int failures = 0;

    void check(bool ok, const std::string& what)
    {
        if (!ok)
        {
            ++failures;
            std::cerr << "FAILED: " << what << std::endl;
        }
    }
};

//! the 1D Lagrange polynomial of the node with index j on the equidistant nodes 0, 1/order, ..., 1
inline double lagrange1D(int order, int j, double xi)
{
    double value = 1.0;
    for (int m = 0; m <= order; ++m)
        if (m != j)
            value *= (xi - double(m)/order)/(double(j - m)/order);
    return value;
}

/*!
 * \brief The Lagrange basis function of order 1 or 2 with the given node, at a point in the reference element
 * \note tensor product polynomials on cubes, barycentric polynomials on simplices
 */
template<class LocalPosition>
double lagrangeBasisFunction(bool simplex, const LocalPosition& node, const LocalPosition& xi)
{
    constexpr int dim = LocalPosition::dimension;
    if (!simplex || dim == 1)
    {
        double value = 1.0;
        for (int k = 0; k < dim; ++k)
            value *= lagrange1D(ORDER, std::lround(node[k]*ORDER), xi[k]);
        return value;
    }

    const auto barycentric = [] (const LocalPosition& x)
    {
        std::array<double, dim + 1> lambda;
        lambda[0] = 1.0;
        for (int k = 0; k < dim; ++k)
        {
            lambda[k + 1] = x[k];
            lambda[0] -= x[k];
        }
        return lambda;
    };

    const auto lambda = barycentric(xi);
    const auto nodeLambda = barycentric(node);
    if constexpr (ORDER == 1)
    {
        for (int a = 0; a <= dim; ++a)
            if (std::abs(nodeLambda[a] - 1.0) < 1e-12)
                return lambda[a];
    }
    else
    {
        std::array<int, 2> vertices{-1, -1};
        int numVertices = 0;
        for (int a = 0; a <= dim; ++a)
            if (nodeLambda[a] > 1e-12)
                vertices[numVertices++] = a;
        if (numVertices == 1)
            return lambda[vertices[0]]*(2.0*lambda[vertices[0]] - 1.0);
        return 4.0*lambda[vertices[0]]*lambda[vertices[1]];
    }
    DUNE_THROW(Dune::InvalidStateException, "Node is not a Lagrange node");
}

//! a random position in the interior of the reference element that is not close to the boundaries of the vertex control volumes
template<class LocalPosition, class Rng>
LocalPosition randomPosition(bool simplex, Rng& rng)
{
    constexpr int dim = LocalPosition::dimension;
    LocalPosition xi(0.0);
    if constexpr (dim > 0)
    {
        std::uniform_real_distribution<double> uniform(0.05, 0.95);
        while (true)
        {
            if (!simplex || dim == 1)
            {
                for (auto& x : xi)
                    x = uniform(rng);
                if (std::all_of(xi.begin(), xi.end(), [] (double x) { return std::abs(x - 0.5) > 0.03; }))
                    return xi;
            }
            else
            {
                std::array<double, dim + 1> lambda;
                double sum = 0.0;
                for (auto& l : lambda)
                {
                    l = -std::log(uniform(rng));
                    sum += l;
                }
                for (auto& l : lambda)
                    l /= sum;
                for (int k = 0; k < dim; ++k)
                    xi[k] = lambda[k + 1];
                auto sorted = lambda;
                std::sort(sorted.begin(), sorted.end());
                if (sorted[dim] - sorted[dim - 1] > 0.03)
                    return xi;
            }
        }
    }
    return xi;
}

/*!
 * \brief The weight that decides which vertex control volume contains a point in the reference element
 * \note The control volume of a vertex is the part of a simplex where its barycentric coordinate is the largest, and the
 *       part of a cube that is closer to it than to the other vertices in every coordinate direction.
 */
template<class LocalPosition>
double controlVolumeWeight(bool simplex, const LocalPosition& vertex, const LocalPosition& xi)
{
    constexpr int dim = LocalPosition::dimension;
    if (simplex && dim > 1)
    {
        double lambda = 1.0;
        for (int k = 0; k < dim; ++k)
            lambda -= xi[k];
        for (int k = 0; k < dim; ++k)
            if (vertex[k] > 0.5)
                lambda = xi[k];
        return lambda;
    }

    double weight = 1.0;
    for (int k = 0; k < dim; ++k)
        weight *= vertex[k] > 0.5 ? double(xi[k] > 0.5) : double(xi[k] < 0.5);
    return weight;
}

//! the local dofs of the bound element that own a control volume
template<class ElementDiscretization>
std::vector<bool> controlVolumeDofs(const ElementDiscretization& elemDisc)
{
    std::vector<bool> isControlVolume(Detail::LocalDofs::numLocalDofs(elemDisc), false);
    if constexpr (Experimental::Concepts::FVElementDiscretization<ElementDiscretization>)
        for (const auto& scv : scvs(elemDisc))
            isControlVolume[scv.localDofIndex()] = true;
    return isControlVolume;
}

} // end namespace Dumux::Test

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

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

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    constexpr int dim = Grid::dimension;
    constexpr bool hasScvs = Experimental::Concepts::FVElementDiscretization<ElementDiscretization>;

    Test::Checker checker;
    std::mt19937 rng(20261113);
    const LocalResidual localResidual(problem.get());

    // the residual of every element at the current solution
    const auto elementResiduals = [&] ()
    {
        std::vector<typename LocalResidual::ElementResidualVector> residuals;
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            elemDisc.bind(element);
            const auto elemVars = localView(gridVariables->curGridVars()).bind(element, elemDisc, x);
            residuals.push_back(localResidual.evalFluxAndSource(element, elemDisc, elemVars));
        }
        return residuals;
    };

    // the residual of an element with the point sources of strength q at the positions xi (in local coordinates) in it:
    // the control volume that contains the point or all control volumes together (weights sum to one) and the other local dofs
    struct Source { LocalPosition xi; double strength; };
    const auto checkResidual = [&] (const std::string& what, const auto& residuals,
                                    const std::vector<std::vector<Source>>& sources, bool singleOwner)
    {
        double maxDifference = 0.0;
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            elemDisc.bind(element);
            const bool simplex = element.type().isSimplex();
            const auto isControlVolume = Test::controlVolumeDofs(elemDisc);

            double sumControlVolumes = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto i = localDof.index();
                const auto node = ipData(elemDisc, localDof).local();

                double expected = 0.0;
                for (const auto& source : sources[eIdx])
                {
                    if (!isControlVolume[i])
                        expected -= source.strength*Test::lagrangeBasisFunction(simplex, node, source.xi);
                    else if (singleOwner)
                    {
                        // the control volume with the largest weight among all vertices of the element
                        bool owns = true;
                        const auto weight = Test::controlVolumeWeight(simplex, node, source.xi);
                        for (const auto& other : localDofs(elemDisc))
                            if (isControlVolume[other.index()]
                                && Test::controlVolumeWeight(simplex, ipData(elemDisc, other).local(), source.xi) > weight)
                                owns = false;
                        if (owns)
                            expected -= source.strength;
                    }
                }

                if (isControlVolume[i])
                    sumControlVolumes += residuals[eIdx][i][0];
                if (!isControlVolume[i] || singleOwner)
                    maxDifference = std::max(maxDifference, std::abs(residuals[eIdx][i][0] - expected));
            }

            if (!singleOwner)
            {
                double total = 0.0;
                for (const auto& source : sources[eIdx])
                    total -= source.strength;
                maxDifference = std::max(maxDifference, hasScvs ? std::abs(sumControlVolumes - total) : 0.0);
            }
        }
        std::cout << what << ": largest difference to the expected residual entries " << maxDifference << std::endl;
        checker.check(maxDifference < 1e-12, "residual entries of " + what);
    };

    // a point source in the interior of each element
    {
        std::vector<GlobalPosition> positions;
        std::vector<double> strengths;
        std::vector<std::vector<Source>> sources(leafGridView.size(0));
        for (const auto& element : elements(leafGridView))
        {
            const auto xi = Test::randomPosition<LocalPosition>(element.type().isSimplex(), rng);
            positions.push_back(element.geometry().global(xi));
            strengths.push_back(1.0 + 0.37*positions.size());
            sources[gridGeometry->elementMapper().index(element)].push_back({xi, strengths.back()});
        }
        problem->setPointSources(positions, strengths);
        checkResidual("points in the interior of elements", elementResiduals(), sources, /*singleOwner=*/true);
    }

    // a point source on every face between two elements, which is split equally between them
    {
        std::vector<GlobalPosition> positions;
        std::vector<double> strengths;
        std::vector<std::vector<Source>> sources(leafGridView.size(0));
        for (const auto& element : elements(leafGridView))
            for (const auto& intersection : intersections(leafGridView, element))
            {
                const auto eIdx = gridGeometry->elementMapper().index(element);
                if (!intersection.neighbor() || gridGeometry->elementMapper().index(intersection.outside()) < eIdx)
                    continue;

                using FacePosition = Dune::FieldVector<double, dim - 1>;
                const auto xiFace = Test::randomPosition<FacePosition>(intersection.type().isSimplex(), rng);
                positions.push_back(intersection.geometry().global(xiFace));
                strengths.push_back(1.0 + 0.37*positions.size());
                const auto nIdx = gridGeometry->elementMapper().index(intersection.outside());
                sources[eIdx].push_back({intersection.geometryInInside().global(xiFace), 0.5*strengths.back()});
                sources[nIdx].push_back({intersection.geometryInOutside().global(xiFace), 0.5*strengths.back()});
            }
        problem->setPointSources(positions, strengths);
        checkResidual("points on faces between elements", elementResiduals(), sources, /*singleOwner=*/false);
    }

    // a single point source in the interior of the domain
    const auto position = getParam<GlobalPosition>("PointSource.Position");
    const auto strength = getParam<double>("PointSource.Strength");
    problem->setPointSources({position}, {strength});

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    // the rows of the control volumes (reactionControlVolumes) and of the other local dofs (reactionOthers)
    SolutionVector reactionControlVolumes(gridGeometry->numDofs()), reactionOthers(gridGeometry->numDofs());
    reactionControlVolumes = 0.0;
    reactionOthers = 0.0;
    {
        const auto residuals = elementResiduals();
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            elemDisc.bind(element);
            const auto isControlVolume = Test::controlVolumeDofs(elemDisc);
            for (const auto& localDof : localDofs(elemDisc))
                (isControlVolume[localDof.index()] ? reactionControlVolumes : reactionOthers)[localDof.dofIndex()]
                    += residuals[eIdx][localDof.index()];
        }
    }

    double boundaryControlVolumes = 0.0, boundaryOthers = 0.0, maxInterior = 0.0;
    for (std::size_t dofIdx = 0; dofIdx < gridGeometry->numDofs(); ++dofIdx)
    {
        if (gridGeometry->dofOnBoundary(dofIdx))
        {
            boundaryControlVolumes += reactionControlVolumes[dofIdx][0];
            boundaryOthers += reactionOthers[dofIdx][0];
        }
        else
            maxInterior = std::max({maxInterior, std::abs(reactionControlVolumes[dofIdx][0]), std::abs(reactionOthers[dofIdx][0])});
    }

    // the outflow is the negative sum of the reactions of all rows that carry a conservation law
    const auto outflow = hasScvs ? -boundaryControlVolumes : -boundaryOthers;
    std::cout << std::setprecision(12)
              << "point source strength " << strength << ", outflow through the boundary " << outflow
              << ", relative difference " << (outflow - strength)/strength << "\n"
              << "largest residual of an interior row " << maxInterior << std::endl;
    checker.check(maxInterior < 1e-9*std::abs(strength), "the interior rows are solved");
    if constexpr (hasScvs)
        checker.check(std::abs(outflow - strength) < 1e-9*std::abs(strength), "the outflow equals the point source strength");

    // in one dimension the Galerkin solution is exact at the vertices
    if constexpr (dim == 1 && !hasScvs)
    {
        const auto x0 = position[0];
        const auto green = [&] (double xPos) { return strength*(xPos <= x0 ? xPos*(1.0 - x0) : x0*(1.0 - xPos)); };

        double maxError = 0.0, maxErrorOfNodes = 0.0;
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            elemDisc.bind(element);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto ip = ipData(elemDisc, localDof);
                const auto error = std::abs(x[localDof.dofIndex()][0] - green(ip.global()[0]));
                const bool isVertex = ip.local()[0] < 1e-12 || ip.local()[0] > 1.0 - 1e-12;
                (isVertex ? maxError : maxErrorOfNodes) = std::max(isVertex ? maxError : maxErrorOfNodes, error);
            }
        }
        std::cout << "largest error at the vertices " << maxError << ", at other nodes " << maxErrorOfNodes << std::endl;
        checker.check(maxError < 1e-12*strength, "the solution is exact at the vertices");
    }

    if (checker.failures > 0)
    {
        std::cerr << checker.failures << " checks failed" << std::endl;
        return 1;
    }

    if (Dune::MPIHelper::instance().rank() == 0)
        Parameters::print();

    return 0;
}
