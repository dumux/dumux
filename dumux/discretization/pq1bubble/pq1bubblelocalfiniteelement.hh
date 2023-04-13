// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Evaluate P1/Q1 basis with bubble function
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_LOCAL_FINITE_ELEMENT_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_LOCAL_FINITE_ELEMENT_HH

#include <array>
#include <vector>
#include <numeric>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

namespace Dumux::Detail {

/*!
 * \brief P1/Q1 + Bubble on the reference element
 *
 * \tparam D Type to represent the field in the domain
 * \tparam R Type to represent the field in the range
 * \tparam dim Dimension of the domain element
 * \tparam typeId The geometry type
 */
template<class D, class R, unsigned int dim, Dune::GeometryType::Id typeId>
class PQ1BubbleLocalBasis
{
    using PQ1FiniteElement = std::conditional_t<
        Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim),
        Dune::LagrangeCubeLocalFiniteElement<D, R, dim, 1>,
        Dune::LagrangeSimplexLocalFiniteElement<D, R, dim, 1>
    >;
    static constexpr std::size_t numDofs
        = Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim) ? (1<<dim)+1 : (dim+1)+1;
public:
    using Traits = Dune::LocalBasisTraits<
        D, dim, Dune::FieldVector<D, dim>,
        R, 1, Dune::FieldVector<R, 1>,
        Dune::FieldMatrix<R, 1, dim>
    >;

    PQ1BubbleLocalBasis()
    : pq1FiniteElement_{}
    {
        // precompute center values for normalization
        const auto& p1Basis = pq1FiniteElement_.localBasis();
        const auto refElement = Dune::referenceElement<typename Traits::DomainFieldType, dim>(type());
        const auto& center = refElement.position(0, 0);
        p1Basis.evaluateFunction(center, pq1AtCenter_);
    }

    /*!
     * \brief Number of shape functions (one for each vertex and one in the element)
     */
    static constexpr unsigned int size()
    { return numDofs; }

    /*!
     * \brief Evaluate all shape functions
     */
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
        out.reserve(size());
        const auto& p1Basis = pq1FiniteElement_.localBasis();
        p1Basis.evaluateFunction(x, out);
        const auto bubble = evaluateBubble_(x);
        out.resize(size());
        out.back() = bubble;
        for (int i = 0; i < numDofs-1; ++i)
            out[i] -= pq1AtCenter_[i]*out.back();
    }

    /*!
     * \brief Evaluate the Jacobians of all shape functions
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
        out.reserve(size());
        const auto& p1Basis = pq1FiniteElement_.localBasis();
        p1Basis.evaluateJacobian(x, out);

        std::vector<typename Traits::RangeType> shapeValues;
        p1Basis.evaluateFunction(x, shapeValues);

        const auto bubbleJacobian = evaluateBubbleJacobian_(x);

        for (int i = 0; i < numDofs-1; ++i)
            for (int k = 0; k < dim; ++k)
                out[i][0][k] -= pq1AtCenter_[i]*bubbleJacobian[0][k];

        out.resize(size());
        out.back() = bubbleJacobian;
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int, dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
        DUNE_THROW(Dune::NotImplemented, "Partial derivatives");
    }

    /*!
     * \brief Evaluate the Jacobians of all shape functions
     * we are actually cubic/quartic but cannot represent all cubic/quartic polynomials
     */
    static constexpr unsigned int order()
    {
        return 1;
    }

    /*!
     * \brief The reference element type
     */
    static constexpr Dune::GeometryType type()
    {
        return { typeId };
    }
private:
    // evaluate bubble function at x
    typename Traits::RangeType evaluateBubble_(const typename Traits::DomainType& x) const
    {
        if constexpr (type() == Dune::GeometryTypes::simplex(dim))
        {
            if constexpr (dim == 2)
                return 27*x[0]*x[1]*(1-x[0]-x[1]);
            else if constexpr (dim == 3)
                return 256*x[0]*x[1]*x[2]*(1-x[0]-x[1]-x[2]);
        }
        else if constexpr (type() == Dune::GeometryTypes::cube(dim))
        {
            if constexpr (dim == 2)
                return 16*x[0]*x[1]*(1-x[0])*(1-x[1]);
            else if constexpr (dim == 3)
                return 64*x[0]*x[1]*x[2]*(1-x[0])*(1-x[1])*(1-x[2]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Bubble function for " << type());
    }

    // evaluate bubble function at x
    typename Traits::JacobianType evaluateBubbleJacobian_(const typename Traits::DomainType& x) const
    {
        if constexpr (type() == Dune::GeometryTypes::simplex(dim))
        {
            if constexpr (dim == 2)
                return {{27*(x[1]*(1-x[0]-x[1]) - x[0]*x[1]),
                         27*(x[0]*(1-x[0]-x[1]) - x[0]*x[1])}};
            else if constexpr (dim == 3)
                return {{256*(x[1]*x[2]*(1-x[0]-x[1]-x[2]) - x[0]*x[1]*x[2]),
                         256*(x[0]*x[2]*(1-x[0]-x[1]-x[2]) - x[0]*x[1]*x[2]),
                         256*(x[0]*x[1]*(1-x[0]-x[1]-x[2]) - x[0]*x[1]*x[2])}};
        }
        else if constexpr (type() == Dune::GeometryTypes::cube(dim))
        {
            if constexpr (dim == 2)
                return {{16*(x[1]*(1-x[0])*(1-x[1]) - x[0]*x[1]*(1-x[1])),
                         16*(x[0]*(1-x[0])*(1-x[1]) - x[0]*x[1]*(1-x[0]))}};
            else if constexpr (dim == 3)
                return {{64*(x[1]*x[2]*(1-x[0])*(1-x[1])*(1-x[2]) - x[0]*x[1]*x[2]*(1-x[1]))*(1-x[2]),
                         64*(x[0]*x[2]*(1-x[0])*(1-x[1])*(1-x[2]) - x[0]*x[1]*x[2]*(1-x[0]))*(1-x[2]),
                         64*(x[0]*x[1]*(1-x[0])*(1-x[1])*(1-x[2]) - x[0]*x[1]*x[2]*(1-x[0]))*(1-x[1])}};
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Bubble function for " << type() << " dim = " << dim);
    }

    PQ1FiniteElement pq1FiniteElement_;
    std::vector<typename Traits::RangeType> pq1AtCenter_;
};

/*!
 * \brief Associations of the P1/Q1 + Bubble degrees of freedom to the facets of the reference element
 * \tparam dim Dimension of the reference element
 * \tparam typeId The geometry type
 */
template<int dim, Dune::GeometryType::Id typeId>
class PQ1BubbleLocalCoefficients
{
    static constexpr std::size_t numDofs
        = Dune::GeometryType{ typeId } == Dune::GeometryTypes::cube(dim) ? (1<<dim)+1 : (dim+1)+1;
public:

    PQ1BubbleLocalCoefficients()
    {
        // Dune::LocalKey is constructed with
        // - the local index with respect to the subentity type
        // - the codim of the subentity the dof is attached to
        // - the local number of the function to evaluate

        // vertex dofs
        for (std::size_t i=0; i<size()-1; i++)
            localKeys_[i] = Dune::LocalKey(i, dim, 0);

        // cell dof
        localKeys_.back() = Dune::LocalKey(0, 0, 0);
    }

    //! Number of coefficients
    static constexpr std::size_t size ()
    { return numDofs; }

    //! Get i-th local key
    const Dune::LocalKey& localKey (std::size_t i) const
    { return localKeys_[i]; }

private:
    std::array<Dune::LocalKey, numDofs> localKeys_;
};

/*!
 * \brief Evaluate the degrees of freedom of a P1 + Bubble basis
 *
 * \tparam LocalBasis The corresponding set of shape functions
 */
template<class LocalBasis>
class PQ1BubbleLocalInterpolation
{
public:
    /*!
     * \brief Evaluate a given function at the vertices and the cell midpoint
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] f Function to evaluate (call operator that gets a local position)
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
        constexpr auto dim = LocalBasis::Traits::dimDomain;

        out.resize(LocalBasis::size());

        const auto refElement = Dune::referenceElement<typename LocalBasis::Traits::DomainFieldType,dim>(LocalBasis::type());

        // Evaluate at the vertices and at the center
        for (int i = 0; i < refElement.size(dim); ++i)
            out[i] = f(refElement.position(i, dim));
        out.back() = f(refElement.position(0, 0));
    }
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \brief P1/Q1 + Bubble finite element
 *
 * \tparam D type used for domain coordinates
 * \tparam R type used for function values
 * \tparam dim dimension of the reference element
 * \tparam typeId The geometry type
 */
template<class D, class R, int dim, Dune::GeometryType::Id typeId>
class PQ1BubbleLocalFiniteElement
{
    using Basis = Detail::PQ1BubbleLocalBasis<D, R, dim, typeId>;
    using Coefficients = Detail::PQ1BubbleLocalCoefficients<dim, typeId>;
    using Interpolation = Detail::PQ1BubbleLocalInterpolation<Basis>;

    static constexpr Dune::GeometryType gt = Dune::GeometryType{ typeId };
    static_assert(
        gt == Dune::GeometryTypes::cube(dim) || gt == Dune::GeometryTypes::simplex(dim),
        "Only implemented for cubes and simplices"
    );

public:
    using Traits = Dune::LocalFiniteElementTraits<Basis, Coefficients, Interpolation>;

    /*!
     * \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
        return basis_;
    }

    /*!
     * \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
        return coefficients_;
    }

    /*!
     * \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
        return interpolation_;
    }

    /*!
     * \brief The number of coefficients in the basis
     */
    static constexpr std::size_t size()
    {
        return Basis::size();
    }

    /*!
     * \brief The reference element type that the local finite element is defined on
     */
    static constexpr Dune::GeometryType type()
    {
        return Traits::LocalBasisType::type();
    }

private:
    Basis basis_;
    Coefficients coefficients_;
    Interpolation interpolation_;
};

} // end namespace Dumux

#endif
