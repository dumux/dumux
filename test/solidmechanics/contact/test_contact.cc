// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Kinematics, law and search of the penalty contact
 *
 * - closest points on segments and triangles satisfy the projection inequality,
 * - the penalty pressure is the negative derivative of its energy density and is C1,
 * - a straight strip has no contact candidates,
 * - a strip rolled into a spiral with the turns a gap apart has that gap, up to the sagitta of
 *   the facets and the blending of nearby facets, with the right sign for separated and for
 *   overlapping turns,
 * - on a randomly perturbed spiral: the residual is the gradient of the contact energy, the
 *   forces have zero sum and zero moment, the Jacobian is the derivative of the residual and is
 *   symmetric, energy, residual and candidates are invariant under a superposed rigid motion, the
 *   tree search finds the same candidates as testing all pairs, a step is covered by the
 *   candidates exactly if no node moves by more than half the search margin, and the contact force
 *   on a node sliding along the opposite wall changes continuously, also across concave corners.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/solidmechanics/contact/closestpoint.hh>
#include <dumux/solidmechanics/contact/penalty.hh>
#include <dumux/solidmechanics/contact/penaltycontact.hh>
#include <dumux/solidmechanics/contact/search.hh>
#include <dumux/solidmechanics/contact/surface.hh>

namespace Dumux {

void check(bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, message);
}

template<int dimWorld, std::size_t numCorners>
void testClosestPoint()
{
    using Point = Dune::FieldVector<double, dimWorld>;
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const auto randomPoint = [&](double scale){ Point p; for (auto& v : p) v = scale*dist(gen); return p; };

    std::array<int, numCorners> numHits{};
    for (int sample = 0; sample < 4000; ++sample)
    {
        std::array<Point, numCorners> c;
        for (auto& corner : c)
            corner = randomPoint(1.0);
        const auto p = randomPoint(2.0);
        const auto cp = Contact::closestPointOnSimplex(p, c);

        Point y(0.0);
        double sum = 0.0;
        for (std::size_t j = 0; j < numCorners; ++j)
        {
            check(cp.weights[j] >= 0.0, "negative weight");
            y.axpy(cp.weights[j], c[j]);
            sum += cp.weights[j];
        }
        check(std::abs(sum - 1.0) < 1e-12, "weights do not sum to one");
        check((y - cp.point).two_norm() < 1e-12, "closest point differs from the weighted corners");
        check(std::abs((p - y).two_norm2() - cp.squaredDistance) < 1e-12, "wrong squared distance");
        // projection onto a convex set: (p - y).(z - y) <= 0 for all z in the simplex
        for (const auto& corner : c)
            check((p - y)*(corner - y) < 1e-10, "not the closest point");
        ++numHits[std::count(cp.weights.begin(), cp.weights.end(), 0.0)];
    }
    for (std::size_t k = 0; k < numCorners; ++k)
        check(numHits[k] > 0, "a feature type was never the closest");
    std::cout << "closest point on " << numCorners << " corners in " << dimWorld << "D checked" << std::endl;
}

void testPenaltyLaw()
{
    const Contact::PenaltyLaw<double> law(3.0, 0.2);
    const double h = 1e-6;
    for (int i = 0; i <= 100; ++i)
    {
        const double g = -0.3 + 0.6*i/100.0 + 1.234e-3;
        const double dPsi = (law.energyDensity(g + h) - law.energyDensity(g - h))/(2.0*h);
        check(std::abs(dPsi + law.pressure(g)) < 1e-8, "pressure is not the negative energy derivative");
        const double dP = (law.pressure(g + h) - law.pressure(g - h))/(2.0*h);
        check(std::abs(dP - law.pressureDerivative(g)) < 1e-8, "wrong pressure derivative");
    }
    for (const double g : {0.0, 0.2})
    {
        check(std::abs(law.pressure(g + 1e-12) - law.pressure(g - 1e-12)) < 1e-10, "pressure jumps");
        check(std::abs(law.pressureDerivative(g + 1e-12) - law.pressureDerivative(g - 1e-12)) < 1e-10, "stiffness jumps");
        check(std::abs(law.energyDensity(g + 1e-12) - law.energyDensity(g - 1e-12)) < 1e-10, "energy jumps");
    }
    check(law.pressure(0.25) == 0.0 && law.energyDensity(0.25) == 0.0 && law.pressure(-0.1) > 0.0, "wrong support");
    std::cout << "penalty law checked" << std::endl;
}

constexpr double length = 1.0, thickness = 0.05, outerRadius = 0.2;
constexpr double turnAngle = 2.6*std::numbers::pi;
constexpr int cellsAlong = 80;

template<int dim>
auto makeStrip()
{
    Dune::FieldVector<double, dim> upperRight(thickness);
    upperRight[0] = length;
    std::array<int, dim> cells;
    cells.fill(2);
    cells[0] = cellsAlong;
    return std::make_unique<Dune::YaspGrid<dim>>(upperRight, cells);
}

// the strip angle and the radius of the reference point X after rolling up
double spiralAngle(double X0)
{ return turnAngle*X0/length; }

double spiralRadius(double X0, double X1, double gap)
{ return outerRadius - (thickness + gap)*spiralAngle(X0)/(2.0*std::numbers::pi) - X1; }

// rolls the strip into a spiral whose consecutive turns are separated by gap (negative: overlap)
template<class GridView, class Mapper>
auto spiralDisplacement(const GridView& gridView, const Mapper& mapper, double gap)
{
    constexpr int dim = GridView::dimension;
    Dune::BlockVector<Dune::FieldVector<double, dim>> u(mapper.size());
    for (const auto& vertex : vertices(gridView))
    {
        const auto X = vertex.geometry().corner(0);
        const double theta = spiralAngle(X[0]);
        const double r = spiralRadius(X[0], X[1], gap);
        auto x = X;
        x[0] = r*std::cos(theta);
        x[1] = r*std::sin(theta);
        u[mapper.index(vertex)] = x - X;
    }
    return u;
}

template<int dim>
void testSpiral()
{
    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LeafGridView;
    using Surface = Contact::ContactSurface<GridView>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, dim>>;
    using Point = Dune::FieldVector<double, dim>;

    const auto grid = makeStrip<dim>();
    const auto gridView = grid->leafGridView();
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(gridView, Dune::mcmgVertexLayout());
    const auto surface = std::make_shared<Surface>(gridView, mapper);
    Contact::PenaltyContact<Surface> contact(surface);
    const double offset = contact.law().offset();
    const double largeRadius = 1.5*thickness;

    // straight strip: no candidates, even with a radius larger than the thickness
    {
        Vector u(mapper.size());
        u = 0.0;
        const auto x = surface->positions(u);
        const auto n = surface->normals(x);
        for (const double radius : {contact.searchRadius(), largeRadius})
        {
            const auto candidates = Contact::findCandidateFacets(*surface, x, radius, contact.candidateFilter(x, n));
            check(candidates == Contact::findCandidateFacetsBruteForce(*surface, x, radius, contact.candidateFilter(x, n)),
                  "tree and brute-force search differ on the straight strip");
            for (const auto& c : candidates)
                check(c.empty(), "candidates on the straight strip");
        }
    }

    // unperturbed spiral: the gap between the turns up to the sagitta of the facets
    const double facetLength = turnAngle/cellsAlong*outerRadius;
    const double minRadius = outerRadius - 1.3*(thickness + offset) - thickness;
    const double sagitta = facetLength*facetLength/(8.0*minRadius);
    check(sagitta < 0.25*offset, "test setup: sagitta too large to distinguish the gaps");
    const double smoothingBias = 2.0*getParam<double>("Contact.SmoothingLength");
    for (const double gap : {0.5*offset, -0.5*offset})
    {
        const auto u = spiralDisplacement(gridView, mapper, gap);
        contact.updateCandidates(u);
        const auto pairs = contact.activePairs(u);
        std::vector<bool> found(surface->numNodes(), false);
        for (const auto& p : pairs)
            found[p.node] = true;

        int numInterior = 0;
        for (std::size_t s = 0; s < surface->numNodes(); ++s)
        {
            const auto& X = surface->referencePosition(s);
            const double theta = spiralAngle(X[0]);
            const bool outerTurnInside = X[1] > thickness - 1e-10 && theta > 0.1*std::numbers::pi && theta < 0.5*std::numbers::pi;
            const bool innerTurnOutside = X[1] < 1e-10 && theta > 2.1*std::numbers::pi && theta < 2.5*std::numbers::pi;
            if (!outerTurnInside && !innerTurnOutside)
                continue;
            ++numInterior;
            check(found[s], "a node facing the next turn has no contact pair");
        }
        check(numInterior > 20, "test setup: too few facing nodes");
        for (const auto& p : pairs)
        {
            const auto& X = surface->referencePosition(p.node);
            const double theta = spiralAngle(X[0]);
            const bool interior = (theta > 0.1*std::numbers::pi && theta < 0.5*std::numbers::pi)
                                  || (theta > 2.1*std::numbers::pi && theta < 2.5*std::numbers::pi);
            if (interior && std::abs(p.gap - gap) >= 1.5*sagitta + smoothingBias)
                DUNE_THROW(Dune::Exception, "gap " << p.gap << " of node " << p.node << " at " << X
                           << " differs from the distance of the turns " << gap);
        }
        std::cout << dim << "D spiral with gap " << gap << ": " << numInterior << " facing nodes, "
                  << pairs.size() << " active pairs" << std::endl;
    }

    // perturbed spiral
    auto u = spiralDisplacement(gridView, mapper, 0.2*offset);
    std::mt19937 gen(11);
    std::uniform_real_distribution<double> dist(-0.3*offset, 0.3*offset);
    for (auto& block : u)
        for (auto& v : block)
            v += dist(gen);
    contact.updateCandidates(u);

    {
        const auto x = surface->positions(u);
        const auto n = surface->normals(x);
        for (const double radius : {contact.searchRadius(), largeRadius})
            check(Contact::findCandidateFacets(*surface, x, radius, contact.candidateFilter(x, n))
                  == Contact::findCandidateFacetsBruteForce(*surface, x, radius, contact.candidateFilter(x, n)),
                  "tree and brute-force search differ");
    }

    const auto pairs = contact.activePairs(u);
    const auto numPenetrated = std::count_if(pairs.begin(), pairs.end(), [](const auto& p){ return p.gap < 0.0; });
    check(numPenetrated > 0 && numPenetrated < static_cast<long>(pairs.size()), "test setup: need separated and penetrated pairs");

    Vector r(mapper.size());
    r = 0.0;
    contact.addResidual(r, u);

    // gradient of the energy
    const double h = 1e-7;
    double maxResidual = 0.0;
    for (const auto& block : r)
        maxResidual = std::max(maxResidual, block.infinity_norm());
    double maxGradientError = 0.0;
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
    {
        const auto dof = surface->dofIndex(s);
        for (int k = 0; k < dim; ++k)
        {
            auto uPlus = u, uMinus = u;
            uPlus[dof][k] += h;
            uMinus[dof][k] -= h;
            const double gradient = (contact.energy(uPlus) - contact.energy(uMinus))/(2.0*h);
            maxGradientError = std::max(maxGradientError, std::abs(gradient - r[dof][k]));
        }
    }
    check(maxGradientError < 1e-6*maxResidual, "residual is not the energy gradient");

    // zero net force and moment
    const auto x = surface->positions(u);
    Point netForce(0.0);
    Dune::FieldVector<double, dim == 2 ? 1 : 3> netMoment(0.0);
    double scale = 0.0;
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
    {
        const auto& f = r[surface->dofIndex(s)];
        netForce += f;
        if constexpr (dim == 2)
            netMoment[0] += x[s][0]*f[1] - x[s][1]*f[0];
        else
            netMoment += crossProduct(x[s], f);
        scale += f.two_norm()*(1.0 + x[s].two_norm());
    }
    check(netForce.two_norm() < 1e-13*scale && netMoment.two_norm() < 1e-13*scale, "net force or moment of the contact forces");

    // Jacobian against differences of the residual, and its symmetry
    std::vector<int> local(mapper.size(), -1);
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
        local[surface->dofIndex(s)] = s;
    const std::size_t size = surface->numNodes()*dim;
    Dune::DynamicMatrix<double> jacobian(size, size, 0.0), difference(size, size, 0.0);
    contact.addJacobian(u, [&](std::size_t row, std::size_t col, const auto& block)
    {
        check(local[row] >= 0 && local[col] >= 0, "Jacobian entry outside the surface");
        for (int i = 0; i < dim; ++i)
            for (int k = 0; k < dim; ++k)
                jacobian[local[row]*dim + i][local[col]*dim + k] += block[i][k];
    });
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
    {
        const auto dof = surface->dofIndex(s);
        for (int k = 0; k < dim; ++k)
        {
            auto uPlus = u, uMinus = u;
            uPlus[dof][k] += h;
            uMinus[dof][k] -= h;
            Vector rPlus(mapper.size()), rMinus(mapper.size());
            rPlus = 0.0;
            rMinus = 0.0;
            contact.addResidual(rPlus, uPlus);
            contact.addResidual(rMinus, uMinus);
            for (std::size_t t = 0; t < surface->numNodes(); ++t)
                for (int i = 0; i < dim; ++i)
                {
                    const auto tDof = surface->dofIndex(t);
                    difference[t*dim + i][s*dim + k] = (rPlus[tDof][i] - rMinus[tDof][i])/(2.0*h);
                }
        }
    }
    double maxJacobian = 0.0, maxJacobianError = 0.0, maxAsymmetry = 0.0;
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = 0; j < size; ++j)
        {
            maxJacobian = std::max(maxJacobian, std::abs(jacobian[i][j]));
            maxJacobianError = std::max(maxJacobianError, std::abs(jacobian[i][j] - difference[i][j]));
            maxAsymmetry = std::max(maxAsymmetry, std::abs(jacobian[i][j] - jacobian[j][i]));
        }
    check(maxJacobianError < 1e-4*maxJacobian, "Jacobian is not the residual derivative");
    check(maxAsymmetry < 1e-4*maxJacobian, "Jacobian is not symmetric");

    // superposed rigid motion
    Dune::FieldMatrix<double, dim, dim> rotation(0.0);
    const double angle = 0.7;
    rotation[0][0] = std::cos(angle); rotation[0][1] = -std::sin(angle);
    rotation[1][0] = std::sin(angle); rotation[1][1] = std::cos(angle);
    if constexpr (dim == 3)
        rotation[2][2] = 1.0;
    if constexpr (dim == 3)
    {
        Dune::FieldMatrix<double, 3, 3> tilt(0.0);
        tilt[0][0] = 1.0;
        tilt[1][1] = std::cos(0.4); tilt[1][2] = -std::sin(0.4);
        tilt[2][1] = std::sin(0.4); tilt[2][2] = std::cos(0.4);
        rotation.leftmultiply(tilt);
    }
    Point translation(0.3);
    Vector uMoved(mapper.size());
    for (const auto& vertex : vertices(gridView))
    {
        const auto dof = mapper.index(vertex);
        const auto X = vertex.geometry().corner(0);
        Point moved;
        rotation.mv(X + u[dof], moved);
        uMoved[dof] = moved + translation - X;
    }
    const double energy = contact.energy(u);
    check(std::abs(contact.energy(uMoved) - energy) < 1e-10*energy, "energy changes under a rigid motion");
    Vector rMoved(mapper.size());
    rMoved = 0.0;
    contact.addResidual(rMoved, uMoved);
    double maxRotationError = 0.0;
    for (std::size_t dof = 0; dof < r.size(); ++dof)
    {
        Point rotated;
        rotation.mv(r[dof], rotated);
        maxRotationError = std::max(maxRotationError, (rotated - rMoved[dof]).two_norm());
    }
    check(maxRotationError < 1e-9*maxResidual, "residual does not rotate with the body");
    std::vector<std::vector<std::size_t>> candidates(surface->numNodes());
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
        candidates[s] = contact.candidates(s);
    contact.updateCandidates(uMoved);
    for (std::size_t s = 0; s < surface->numNodes(); ++s)
        check(candidates[s] == contact.candidates(s), "candidates change under a rigid motion");

    // a step is covered by the candidates if no surface node moves by more than half the search margin
    {
        const double margin = contact.searchRadius() - offset;
        auto uStep = u;
        uStep[surface->dofIndex(0)][0] += 0.49*margin;
        check(contact.candidatesCover(u, uStep), "a step below half the margin is not covered");
        uStep[surface->dofIndex(0)][0] += 0.02*margin;
        check(!contact.candidatesCover(u, uStep), "a step beyond half the margin is covered");
    }

    // continuity: every active node slides along its own surface over one facet length
    double maxJump = 0.0;
    std::size_t worstNode = 0;
    {
        contact.updateCandidates(u);
        const auto n = surface->normals(surface->positions(u));
        const int numSamples = 800;
        const double slideLength = length/cellsAlong;
        for (const auto& p : contact.activePairs(u))
        {
            const auto dof = surface->dofIndex(p.node);
            Point tangent(0.0);
            tangent[0] = -n.node[p.node][1];
            tangent[1] = n.node[p.node][0];
            if (tangent.two_norm() < 1e-8)
                continue;
            tangent /= tangent.two_norm();
            Point previous(0.0);
            double maxForce = 0.0, maxStep = 0.0;
            for (int i = 0; i <= numSamples; ++i)
            {
                auto uSlid = u;
                uSlid[dof].axpy(slideLength*(double(i)/numSamples - 0.5), tangent);
                Vector rSlid(mapper.size());
                rSlid = 0.0;
                contact.addResidual(rSlid, uSlid);
                maxForce = std::max(maxForce, rSlid[dof].two_norm());
                if (i > 0)
                    maxStep = std::max(maxStep, (rSlid[dof] - previous).two_norm());
                previous = rSlid[dof];
            }
            if (maxForce > 0.0 && maxStep/maxForce > maxJump)
            {
                maxJump = maxStep/maxForce;
                worstNode = p.node;
            }
        }
        if (maxJump >= 0.1)
            DUNE_THROW(Dune::Exception, "contact force jumps by " << maxJump << " of its maximum while node " << worstNode
                       << " at " << surface->referencePosition(worstNode) << " slides along the opposite wall");
    }

    std::cout << dim << "D perturbed spiral: " << pairs.size() << " active pairs (" << numPenetrated << " penetrated), "
              << "gradient error " << maxGradientError/maxResidual << ", Jacobian error " << maxJacobianError/maxJacobian
              << ", asymmetry " << maxAsymmetry/maxJacobian << ", rigid motion error " << maxRotationError/maxResidual << ", largest force change while sliding " << maxJump << std::endl;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    Dumux::Parameters::init(argc, argv, [](auto& p){
        p["Contact.PenaltyStiffness"] = "100";
        p["Contact.Offset"] = "0.008";
        p["Contact.SmoothingLength"] = "4e-4";
    });

    Dumux::testClosestPoint<2, 2>();
    Dumux::testClosestPoint<3, 2>();
    Dumux::testClosestPoint<2, 3>();
    Dumux::testClosestPoint<3, 3>();
    Dumux::testPenaltyLaw();
    Dumux::testSpiral<2>();
    Dumux::testSpiral<3>();
    return 0;
}
