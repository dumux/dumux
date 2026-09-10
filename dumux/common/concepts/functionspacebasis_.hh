// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts for function space bases usable in L2-projections
 * \note The nested-type vocabulary required here, a local view with a tree carrying a finite
 *       element, follows the global basis interface of dune-functions.
 */
#ifndef DUMUX_CONCEPTS_FUNCTIONSPACEBASIS__HH
#define DUMUX_CONCEPTS_FUNCTIONSPACEBASIS__HH

#include <cstddef>
#include <concepts>

namespace Dumux::Concept {

/*!
 * \ingroup Concepts
 * \brief A local finite element with scalar-valued shape functions.
 *
 * The projection assembly stores shape values in scalar-valued buffers, so a
 * vector-valued local basis cannot be used without evaluating componentwise.
 */
template<class FE>
concept ScalarLocalFiniteElement = requires
{
    typename FE::Traits::LocalBasisType;
    typename FE::Traits::LocalBasisType::Traits::RangeFieldType;
    typename FE::Traits::LocalBasisType::Traits::RangeType;
    requires (FE::Traits::LocalBasisType::Traits::RangeType::dimension == 1);
};

/*!
 * \ingroup Concepts
 * \brief A function space basis that an L2-projection can assemble against.
 *
 * The dimension is read either from a `dimension` member or, for bases defined
 * over a grid, from the grid view. Binding is not part of this concept: the
 * entity type a local view accepts is determined by the glue the projection runs
 * over, not by the basis alone. See ProjectableOverGlue.
 */
template<class B>
concept ProjectionBasis = requires(const B& basis)
{
    typename B::LocalView;
    typename B::LocalView::Tree;
    typename B::LocalView::Tree::FiniteElement;
    requires ScalarLocalFiniteElement<typename B::LocalView::Tree::FiniteElement>;
    { basis.size() } -> std::convertible_to<std::size_t>;
    { basis.localView() } -> std::same_as<typename B::LocalView>;
}
&& (requires { { B::dimension } -> std::convertible_to<int>; }
    || requires { { B::GridView::dimension } -> std::convertible_to<int>; })
&& requires(typename B::LocalView localView, std::size_t i)
{
    { localView.tree() } -> std::same_as<const typename B::LocalView::Tree&>;
    { localView.index(i) } -> std::convertible_to<std::size_t>;
};

/*!
 * \ingroup Concepts
 * \brief A basis whose entities the pattern of a mass matrix can be built from.
 *
 * Satisfied either by a basis over a grid view or by one over a bare entity set.
 */
template<class B>
concept EntityRangeProvider = requires(const B& basis)
{
    requires requires { basis.gridView(); } || requires { basis.entitySet(); };
};

/*!
 * \ingroup Concepts
 * \brief A pair of bases whose local views bind to the entities of a given glue.
 *
 * This is what makes a projection well-formed at the call site rather than deep
 * inside the assembly: the entity types come from the glue, so they cannot be
 * checked against either basis in isolation.
 */
template<class DomainBasis, class TargetBasis, class GlueType>
concept ProjectableOverGlue =
    ProjectionBasis<DomainBasis> && ProjectionBasis<TargetBasis> &&
    requires(typename DomainBasis::LocalView domainLocalView,
             typename TargetBasis::LocalView targetLocalView,
             const typename GlueType::Entity& is)
    {
        domainLocalView.bind(is.domainEntity(0));
        targetLocalView.bind(is.targetEntity(0));
    };

} // end namespace Dumux::Concept

#endif
