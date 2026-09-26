// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test L2-projections between broken spaces on non-matching, topology-less entity sets.
 */
#include <config.h>

#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectionentityset.hh>
#include <dumux/discretization/projection/brokenbasis.hh>
#include <dumux/discretization/projection/projector.hh>

namespace {

using Segment = Dune::MultiLinearGeometry<double, 1, 2>;
using SegmentSet = Dumux::GeometriesEntitySet<Segment>;
using Coefficients = Dune::BlockVector<Dune::FieldVector<double, 1>>;

//! Subdivide the unit segment along the x-axis into n pieces
auto makeInterface(std::size_t n)
{
    using P = Dune::FieldVector<double, 2>;
    std::vector<Segment> segments;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double x0 = double(i)/double(n);
        const double x1 = double(i + 1)/double(n);
        segments.emplace_back(Dune::GeometryTypes::line, std::vector<P>{{x0, 0.0}, {x1, 0.0}});
    }
    return std::make_shared<SegmentSet>(std::move(segments));
}

auto makeGlue(std::shared_ptr<const SegmentSet> domain, std::shared_ptr<const SegmentSet> target)
{
    Dumux::IntersectionEntitySet<SegmentSet, SegmentSet> glue;
    glue.build(domain, target);
    return glue;
}

/*!
 * \brief Interpolate an analytic function into a broken basis.
 * \note Lagrange bases are nodal, so interpolation is evaluation at the local nodes.
 */
template<class Basis, class F>
Coefficients interpolate(const Basis& basis, const F& f)
{
    Coefficients u(basis.size());
    u = 0.0;
    auto localView = basis.localView();
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto& fe = localView.tree().finiteElement();
        const auto geometry = entity.geometry();
        std::vector<double> nodalValues;
        fe.localInterpolation().interpolate(
            [&] (const auto& local) { return f(geometry.global(local)); }, nodalValues
        );
        for (std::size_t i = 0; i < nodalValues.size(); ++i)
            u[localView.index(i)] = nodalValues[i];
    }
    return u;
}

/*!
 * \brief Interpolate an analytic function into a broken basis, shifted by the entity index.
 *
 * The result lies in the basis by construction but is discontinuous between entities,
 * so it lies in no conforming space and in no broken space over another subdivision.
 */
template<class Basis, class F>
Coefficients interpolateBroken(const Basis& basis, const F& f)
{
    Coefficients u(basis.size());
    u = 0.0;
    auto localView = basis.localView();
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto jump = double(basis.entitySet().index(entity));
        const auto geometry = entity.geometry();
        std::vector<double> nodalValues;
        localView.tree().finiteElement().localInterpolation().interpolate(
            [&] (const auto& local) { return f(geometry.global(local)) + jump; }, nodalValues
        );
        for (std::size_t i = 0; i < nodalValues.size(); ++i)
            u[localView.index(i)] = nodalValues[i];
    }
    return u;
}

/*!
 * \brief The L2 projection of a discrete function of the domain basis onto the target
 *        basis, assembled independently of the projector over the same intersections with
 *        a quadrature rule far beyond what any of the integrands requires.
 */
template<class DomainBasis, class TargetBasis, class Glue>
Dune::DynamicVector<double> referenceProjection(const DomainBasis& domainBasis,
                                                const TargetBasis& targetBasis,
                                                const Glue& glue,
                                                const Coefficients& u)
{
    const auto n = targetBasis.size();
    Dune::DynamicMatrix<double> massMatrix(n, n, 0.0);
    Dune::DynamicVector<double> rhs(n, 0.0);
    auto targetLocalView = targetBasis.localView();
    auto domainLocalView = domainBasis.localView();
    std::vector<Dune::FieldVector<double, 1>> targetValues, domainValues;
    for (const auto& is : intersections(glue))
    {
        if (is.numDomainNeighbors() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Expected one domain neighbour per intersection");

        const auto& targetEntity = is.targetEntity(0);
        const auto& domainEntity = is.domainEntity(0);
        targetLocalView.bind(targetEntity);
        domainLocalView.bind(domainEntity);
        const auto& targetLocalBasis = targetLocalView.tree().finiteElement().localBasis();
        const auto& domainLocalBasis = domainLocalView.tree().finiteElement().localBasis();

        const auto isGeometry = is.geometry();
        const auto& quad = Dune::QuadratureRules<double, 1>::rule(isGeometry.type(), 20);
        for (const auto& qp : quad)
        {
            const auto weight = qp.weight()*isGeometry.integrationElement(qp.position());
            const auto global = isGeometry.global(qp.position());
            targetLocalBasis.evaluateFunction(targetEntity.geometry().local(global), targetValues);
            domainLocalBasis.evaluateFunction(domainEntity.geometry().local(global), domainValues);

            double sourceValue = 0.0;
            for (std::size_t j = 0; j < domainValues.size(); ++j)
                sourceValue += u[domainLocalView.index(j)]*domainValues[j][0];

            for (std::size_t i = 0; i < targetValues.size(); ++i)
            {
                rhs[targetLocalView.index(i)] += weight*targetValues[i][0]*sourceValue;
                for (std::size_t j = 0; j < targetValues.size(); ++j)
                    massMatrix[targetLocalView.index(i)][targetLocalView.index(j)]
                        += weight*targetValues[i][0]*targetValues[j][0];
            }
        }
    }

    Dune::DynamicVector<double> reference(n, 0.0);
    massMatrix.solve(reference, rhs);
    return reference;
}

//! Largest coefficient difference between a projection and its reference
double maxDifference(const Coefficients& projected, const Dune::DynamicVector<double>& reference)
{
    using std::abs;
    using std::max;
    double diff = 0.0;
    for (std::size_t i = 0; i < projected.size(); ++i)
        diff = max(diff, abs(projected[i][0] - reference[i]));
    return diff;
}

//! Evaluate a discrete function of a broken basis at a global position inside a given entity
template<class Basis>
double evaluate(const Basis& basis, const Coefficients& u, const auto& entity, const auto& globalPos)
{
    auto localView = basis.localView();
    localView.bind(entity);
    const auto& localBasis = localView.tree().finiteElement().localBasis();
    std::vector<Dune::FieldVector<double, 1>> shapeValues;
    localBasis.evaluateFunction(entity.geometry().local(globalPos), shapeValues);
    double result = 0.0;
    for (std::size_t i = 0; i < shapeValues.size(); ++i)
        result += u[localView.index(i)]*shapeValues[i][0];
    return result;
}

//! Integral of a discrete function over the whole entity set
template<class Basis>
double integrate(const Basis& basis, const Coefficients& u)
{
    double integral = 0.0;
    auto localView = basis.localView();
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto& localBasis = localView.tree().finiteElement().localBasis();
        const auto geometry = entity.geometry();
        const auto& quad = Dune::QuadratureRules<double, 1>::rule(geometry.type(), 2*localBasis.order() + 2);
        std::vector<Dune::FieldVector<double, 1>> shapeValues;
        for (const auto& qp : quad)
        {
            localBasis.evaluateFunction(qp.position(), shapeValues);
            double value = 0.0;
            for (std::size_t i = 0; i < shapeValues.size(); ++i)
                value += u[localView.index(i)]*shapeValues[i][0];
            integral += value*qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    return integral;
}

//! Largest deviation of a discrete function from an analytic one, sampled per entity
template<class Basis, class F>
double maxDeviation(const Basis& basis, const Coefficients& u, const F& f)
{
    using std::abs;
    using std::max;
    double error = 0.0;
    for (const auto& entity : basis.entitySet())
    {
        const auto geometry = entity.geometry();
        for (const double x : {0.0, 0.23, 0.5, 0.77, 1.0})
        {
            const auto global = geometry.global({x});
            error = max(error, abs(evaluate(basis, u, entity, global) - f(global)));
        }
    }
    return error;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    // deliberately non-matching subdivisions of the same interface
    const auto domainSet = makeInterface(3);
    const auto targetSet = makeInterface(2);
    const auto glue = makeGlue(domainSet, targetSet);

    if (glue.size() == 0)
        DUNE_THROW(Dune::InvalidStateException, "Glue contains no intersections");

    const auto affine = [] (const auto& p) { return 0.3 + 1.7*p[0]; };
    const auto quadratic = [] (const auto& p) { return 0.2 - 0.9*p[0] + 2.1*p[0]*p[0]; };

    // ------------------------------------------------------------------
    // order 1: an affine function lies in both spaces and must be reproduced
    // ------------------------------------------------------------------
    {
        const BrokenLagrangeBasis<SegmentSet, 1> domainBasis{domainSet};
        const BrokenLagrangeBasis<SegmentSet, 1> targetBasis{targetSet};

        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;

        const auto u = interpolate(domainBasis, affine);
        const auto projected = projector.project(u, params);

        if (projected.size() != targetBasis.size())
            DUNE_THROW(Dune::InvalidStateException, "Projected vector has wrong size");

        const auto error = maxDeviation(targetBasis, projected, affine);
        if (error > 1e-12)
            DUNE_THROW(Dune::MathError, "P1 affine reproduction failed, error = " << error);
        std::cout << "P1 -> P1 affine reproduction: error = " << error << std::endl;
    }

    // ------------------------------------------------------------------
    // order 2: reproduction of a quadratic, which lies in both spaces. This checks
    // the index maps and the intersection partition; it is insensitive to the
    // quadrature order, which the dedicated check below covers.
    // ------------------------------------------------------------------
    {
        const BrokenLagrangeBasis<SegmentSet, 2> domainBasis{domainSet};
        const BrokenLagrangeBasis<SegmentSet, 2> targetBasis{targetSet};

        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;

        const auto u = interpolate(domainBasis, quadratic);
        const auto projected = projector.project(u, params);

        const auto error = maxDeviation(targetBasis, projected, quadratic);
        if (error > 1e-12)
            DUNE_THROW(Dune::MathError, "P2 quadratic reproduction failed, error = " << error);
        std::cout << "P2 -> P2 quadratic reproduction: error = " << error << std::endl;
    }

    // ------------------------------------------------------------------
    // integral preservation for a function NOT in the target space:
    // constants lie in the target space, so the projection is integral preserving
    // ------------------------------------------------------------------
    {
        const BrokenLagrangeBasis<SegmentSet, 2> domainBasis{domainSet};
        const BrokenLagrangeBasis<SegmentSet, 1> targetBasis{targetSet};

        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;

        const auto u = interpolate(domainBasis, quadratic);
        const auto projected = projector.project(u, params);

        const auto before = integrate(domainBasis, u);
        const auto after = integrate(targetBasis, projected);
        using std::abs;
        if (abs(before - after) > 1e-12)
            DUNE_THROW(Dune::MathError,
                       "Integral not preserved: " << before << " vs " << after);
        std::cout << "P2 -> P1 integral preservation: " << before << " vs " << after << std::endl;

        // the quadratic is not in the P1 target space, so this must NOT be a reproduction
        const auto error = maxDeviation(targetBasis, projected, quadratic);
        if (error < 1e-10)
            DUNE_THROW(Dune::InvalidStateException,
                       "Quadratic unexpectedly reproduced in a P1 space; the test is vacuous");
    }

    // ------------------------------------------------------------------
    // Quadrature exactness. A reproduction test cannot detect under-integration:
    // when the source lies in the target space, the same rule appears in the mass
    // and the projection matrix and cancels exactly. So the source here is
    // discontinuous across the domain mesh, putting it in the domain space but not
    // the target's, and the result is compared against the same assembly carried
    // out with a far higher quadrature order.
    // ------------------------------------------------------------------
    const auto checkExactness = [&] (const auto& domainBasis, const auto& targetBasis,
                                     const auto& source, const std::string& name)
    {
        const auto u = interpolateBroken(domainBasis, source);
        const auto reference = referenceProjection(domainBasis, targetBasis, glue, u);

        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;
        const auto diff = maxDifference(projector.project(u, params), reference);
        if (diff > 1e-12)
            DUNE_THROW(Dune::MathError,
                       name << " projection of a broken source differs from the high-order "
                       "reference by " << diff << ", indicating under-integration");
        std::cout << name << " broken source against high-order reference: difference = " << diff << std::endl;
        return reference;
    };

    // with both spaces at order 2 every intersection integrand is of degree 4
    checkExactness(BrokenLagrangeBasis<SegmentSet, 2>{domainSet},
                   BrokenLagrangeBasis<SegmentSet, 2>{targetSet}, quadratic, "P2 -> P2");

    // with the target at order 2 and the domain at order 1, the mixed integrand is of
    // degree 3 but the target mass matrix of degree 4: a rule chosen for the mixed
    // term alone is not enough
    {
        const BrokenLagrangeBasis<SegmentSet, 1> domainBasis{domainSet};
        const BrokenLagrangeBasis<SegmentSet, 2> targetBasis{targetSet};
        const auto reference = checkExactness(domainBasis, targetBasis, affine, "P1 -> P2");

        // the check has teeth only if a rule of the mixed term's degree is caught
        const auto u = interpolateBroken(domainBasis, affine);
        const auto underIntegrating = makeProjector(domainBasis, targetBasis, glue, 3);
        auto params = underIntegrating.defaultParams();
        params.residualReduction = 1e-16;
        const auto diff = maxDifference(underIntegrating.project(u, params), reference);
        if (diff < 1e-8)
            DUNE_THROW(Dune::MathError,
                       "A rule of order 3 was expected to under-integrate the P2 mass matrix, "
                       "but the projection differs from the reference only by " << diff);
        std::cout << "P1 -> P2 with a rule of order 3 differs from the reference by " << diff
                  << ", as it must" << std::endl;
    }

    // ------------------------------------------------------------------
    // idempotence: projecting an already-projected function changes nothing
    // ------------------------------------------------------------------
    {
        const BrokenLagrangeBasis<SegmentSet, 1> targetBasis{targetSet};
        const auto selfGlue = makeGlue(targetSet, targetSet);
        const auto projector = makeProjector(targetBasis, targetBasis, selfGlue);
        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;

        const auto u = interpolate(targetBasis, quadratic);
        const auto once = projector.project(u, params);
        const auto twice = projector.project(once, params);

        double diff = 0.0;
        using std::abs;
        using std::max;
        for (std::size_t i = 0; i < once.size(); ++i)
            diff = max(diff, abs(once[i][0] - twice[i][0]));
        if (diff > 1e-12)
            DUNE_THROW(Dune::MathError, "Projection is not idempotent, difference = " << diff);
        std::cout << "Idempotence: difference = " << diff << std::endl;
    }

    std::cout << "All broken-space projection checks passed" << std::endl;
    return 0;
}
