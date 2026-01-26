// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief P2/Q2 hierarchical local finite element
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_HIERARCHICAL_LOCAL_FINITE_ELEMENT_HH
#define DUMUX_DISCRETIZATION_PQ2_HIERARCHICAL_LOCAL_FINITE_ELEMENT_HH

#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

namespace Dumux::Detail {

template<class D, class R, unsigned int dim, Dune::GeometryType::Id typeId>
class PQ2HierarchicalLocalBasis
{
    using PQ1FE = std::conditional_t<
        Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim),
        Dune::LagrangeCubeLocalFiniteElement<D, R, dim, 1>,
        Dune::LagrangeSimplexLocalFiniteElement<D, R, dim, 1>
    >;
    using PQ2FE = std::conditional_t<
        Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim),
        Dune::LagrangeCubeLocalFiniteElement<D, R, dim, 2>,
        Dune::LagrangeSimplexLocalFiniteElement<D, R, dim, 2>
    >;

public:
    using Traits = Dune::LocalBasisTraits<
        D, dim, Dune::FieldVector<D, dim>,
        R, 1, Dune::FieldVector<R, 1>,
        Dune::FieldMatrix<R, 1, dim>
    >;

    PQ2HierarchicalLocalBasis()
    {
        const auto refElement = Dune::referenceElement<D, dim>(type());
        numVertices_ = refElement.size(dim);
    }

    static constexpr unsigned int size()
    { return PQ2FE{}.size(); }

    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
        out.resize(size());

        // Evaluate both P1 and P2 bases
        std::vector<typename Traits::RangeType> pq1Values, pq2Values;
        pq1FE_.localBasis().evaluateFunction(x, pq1Values);
        pq2FE_.localBasis().evaluateFunction(x, pq2Values);

        // Hierarchical splitting
        const auto& pq2Coeff = pq2FE_.localCoefficients();
        for (std::size_t i = 0; i < size(); ++i)
        {
            const auto& key = pq2Coeff.localKey(i);
            if (key.codim() == dim) // vertex DOF
                out[i] = pq1Values[key.subEntity()];
            else // remaining DOF
                out[i] = pq2Values[i];
        }
    }

    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
        out.resize(size());

        std::vector<typename Traits::JacobianType> pq1Jac, pq2Jac;
        pq1FE_.localBasis().evaluateJacobian(x, pq1Jac);
        pq2FE_.localBasis().evaluateJacobian(x, pq2Jac);

        // For each P2 DOF, use P1 if it's a vertex, otherwise use P2
        const auto& pq2Coeff = pq2FE_.localCoefficients();
        for (std::size_t i = 0; i < size(); ++i)
        {
            const auto& key = pq2Coeff.localKey(i);
            if (key.codim() == dim) // vertex DOF
                out[i] = pq1Jac[key.subEntity()];
            else // remaining DOF
                out[i] = pq2Jac[i];
        }
    }

    void partial(const std::array<unsigned int, dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
        out.resize(size());
        std::vector<typename Traits::RangeType> pq1Partial, pq2Partial;
        pq1FE_.localBasis().partial(order, in, pq1Partial);
        pq2FE_.localBasis().partial(order, in, pq2Partial);

        const auto& pq2Coeff = pq2FE_.localCoefficients();
        for (std::size_t i = 0; i < size(); ++i)
        {
            const auto& key = pq2Coeff.localKey(i);
            if (key.codim() == dim) // vertex DOF
                out[i] = pq1Partial[key.subEntity()];
            else // remaining DOF
                out[i] = pq2Partial[i];
        }
    }

    static constexpr unsigned int order()
    { return 2; }

    static constexpr Dune::GeometryType type()
    { return { typeId }; }

private:
    PQ1FE pq1FE_;
    PQ2FE pq2FE_;
    std::size_t numVertices_;
};

template<int dim, Dune::GeometryType::Id typeId>
class PQ2HierarchicalLocalCoefficients
{
    using PQ2FE = std::conditional_t<
        Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim),
        Dune::LagrangeCubeLocalFiniteElement<double, double, dim, 2>,
        Dune::LagrangeSimplexLocalFiniteElement<double, double, dim, 2>
    >;

public:
    PQ2HierarchicalLocalCoefficients()
    : p2Coeff_(PQ2FE{}.localCoefficients())
    {}

    static constexpr std::size_t size()
    { return PQ2FE{}.size(); }

    const Dune::LocalKey& localKey(std::size_t i) const
    { return p2Coeff_.localKey(i); }

    static constexpr Dune::GeometryType type()
    { return { typeId }; }

private:
    typename PQ2FE::Traits::LocalCoefficientsType p2Coeff_;
};

template<class LocalBasis>
class PQ2HierarchicalLocalInterpolation
{
public:
    template<typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
        DUNE_THROW(Dune::NotImplemented, "Local interpolation not implemented for PQ2 basis.");
    }
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \brief Hierarchical P2/Q2 finite element
 *
 * Decomposition:
 * - Vertices: P1/Q1 nodal basis functions
 * - Remaining dofs: P2/Q2 nodal basis functions
 *
 * \note This represents a hierarchical splitting but does not give the classical one
 *       with edge bubble functions. Instead, the edge functions are the standard P2/Q2
 *       basis functions..
 */
template<class D, class R, int dim, Dune::GeometryType::Id typeId>
class PQ2HierarchicalLocalFiniteElement
{
    using Basis = Detail::PQ2HierarchicalLocalBasis<D, R, dim, typeId>;
    using Coefficients = Detail::PQ2HierarchicalLocalCoefficients<dim, typeId>;
    using Interpolation = Detail::PQ2HierarchicalLocalInterpolation<Basis>;

    static constexpr Dune::GeometryType gt = Dune::GeometryType{ typeId };
    static_assert(
        gt == Dune::GeometryTypes::cube(dim) || gt == Dune::GeometryTypes::simplex(dim),
        "Only implemented for cubes and simplices"
    );

public:
    using Traits = Dune::LocalFiniteElementTraits<Basis, Coefficients, Interpolation>;

    const typename Traits::LocalBasisType& localBasis() const
    { return basis_; }

    const typename Traits::LocalCoefficientsType& localCoefficients() const
    { return coefficients_; }

    const typename Traits::LocalInterpolationType& localInterpolation() const
    { return interpolation_; }

    static constexpr std::size_t size()
    { return Basis::size(); }

    static constexpr Dune::GeometryType type()
    { return Traits::LocalBasisType::type(); }

private:
    Basis basis_;
    Coefficients coefficients_;
    Interpolation interpolation_;
};

} // end namespace Dumux

#endif
