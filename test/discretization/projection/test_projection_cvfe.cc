// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief L2 projection between P1 and P2 CVFE spaces on non-matching grids against a high-order reference
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/projection/l2_projection.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/multidomain/glue.hh>

namespace {

using Coefficients = Dune::BlockVector<Dune::FieldVector<double, 1>>;

//! Nodal interpolation of an analytic function into a Lagrange basis
template<class Basis, class F>
Coefficients interpolate(const Basis& basis, const F& f)
{
    Coefficients u(basis.size());
    u = 0.0;
    auto localView = basis.localView();
    for (const auto& element : elements(basis.gridView()))
    {
        localView.bind(element);
        const auto geometry = element.geometry();
        std::vector<double> nodalValues;
        localView.tree().finiteElement().localInterpolation().interpolate(
            [&] (const auto& local) { return f(geometry.global(local)); }, nodalValues
        );
        for (std::size_t i = 0; i < nodalValues.size(); ++i)
            u[localView.index(i)] = nodalValues[i];
    }
    return u;
}

//! L2 projection of a domain-basis function onto the target basis, assembled with an order-20 rule
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
        const auto& quad = Dune::QuadratureRules<double, 2>::rule(isGeometry.type(), 20);
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

//! Largest coefficient difference between two discrete functions of the same basis
double maxDifference(const Coefficients& a, const Coefficients& b)
{
    using std::abs;
    using std::max;
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        diff = max(diff, abs(a[i][0] - b[i][0]));
    return diff;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;
    using P1GridGeometry = BoxFVGridGeometry<double, GridView>;
    using P2GridGeometry = PQ2FVGridGeometry<double, GridView>;
    using P1Basis = FEBasisFromCVFEGridDiscretization<P1GridGeometry>;
    using P2Basis = FEBasisFromCVFEGridDiscretization<P2GridGeometry>;

    // different subdivisions so that the P1 space is not a subspace of the P2 space
    Grid domainGrid{{1.0, 1.0}, {3, 3}};
    Grid targetGrid{{1.0, 1.0}, {2, 2}};
    P1GridGeometry domainGridGeometry{domainGrid.leafGridView()};
    P2GridGeometry targetGridGeometry{targetGrid.leafGridView()};
    const P1Basis domainBasis{domainGridGeometry};
    const P2Basis targetBasis{targetGridGeometry};

    const auto glue = makeGlue(domainGridGeometry, targetGridGeometry);
    if (glue.size() == 0)
        DUNE_THROW(Dune::InvalidStateException, "Glue contains no intersections");

    const auto affine = [] (const auto& p) { return 0.3 + 1.7*p[0] - 0.4*p[1]; };
    const auto smooth = [] (const auto& p) { return std::sin(3.0*p[0])*std::cos(2.0*p[1]) + p[0]*p[1]; };

    using L2Projector = decltype(makeProjector(domainBasis, targetBasis, glue));
    auto params = L2Projector::defaultParams();
    params.residualReduction = 1e-16;

    // an affine function lies in both spaces and must be reproduced whatever the quadrature
    {
        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        const auto projected = projector.project(interpolate(domainBasis, affine), params);
        if (projected.size() != targetBasis.size())
            DUNE_THROW(Dune::InvalidStateException, "Projected vector has wrong size");
        const auto error = maxDifference(projected, interpolate(targetBasis, affine));
        if (error > 1e-12)
            DUNE_THROW(Dune::MathError, "P1 -> P2 affine reproduction failed, error = " << error);
        std::cout << "P1 -> P2 affine reproduction: error = " << error << std::endl;
    }

    // the P1 interpolant of a smooth function is not in the P2 space, so the default rule must match the reference
    {
        const auto u = interpolate(domainBasis, smooth);
        const auto reference = referenceProjection(domainBasis, targetBasis, glue, u);

        const auto projector = makeProjector(domainBasis, targetBasis, glue);
        const auto diff = maxDifference(projector.project(u, params), reference);
        if (diff > 1e-12)
            DUNE_THROW(Dune::MathError,
                       "P1 -> P2 projection differs from the high-order reference by " << diff
                       << ", indicating under-integration");
        std::cout << "P1 -> P2 against high-order reference: difference = " << diff << std::endl;

        // order 4, twice the basis order, under-integrates the Q2 mass matrix on the triangulated overlaps
        const auto underIntegrating = makeProjector(domainBasis, targetBasis, glue, 4);
        const auto diffOrder4 = maxDifference(underIntegrating.project(u, params), reference);
        if (diffOrder4 < 1e-8)
            DUNE_THROW(Dune::MathError,
                       "A rule of order 4 was expected to under-integrate the P2 mass matrix, "
                       "but the projection differs from the reference only by " << diffOrder4);
        std::cout << "P1 -> P2 with a rule of order 4 differs from the reference by " << diffOrder4
                  << ", as it must" << std::endl;
    }

    // an explicit order 20 must reproduce the default rule
    {
        const auto u = interpolate(domainBasis, smooth);
        const auto byDefault = makeProjector(domainBasis, targetBasis, glue).project(u, params);
        const auto byOverride = makeProjector(domainBasis, targetBasis, glue, 20).project(u, params);
        const auto diff = maxDifference(byDefault, byOverride);
        if (diff > 1e-12)
            DUNE_THROW(Dune::MathError, "Default and order-20 projections differ by " << diff);
        std::cout << "Default order against order 20: difference = " << diff << std::endl;
    }

    std::cout << "All CVFE projection checks passed" << std::endl;
    return 0;
}
