// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Point source context of the embedded coupling for problems with point sources
 *        of the element discretization interface
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCECONTEXT_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCECONTEXT_HH

#include <cstddef>
#include <utility>

#include <dumux/common/pointsources.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The context of a point source of the coupling in one element
 *
 * The context is located in an element and carries the identifier of the point source,
 * with which the coupling manager gives access to the data of the coupling, and the weight of the
 * source: the quadrature weight times the integration element of the integration point.
 * A point shared by several elements is given as one context per element with the number of
 * these elements as embeddings.
 */
template<class GlobalPosition>
class PointSourceContext
: public Experimental::PointSourceContext<GlobalPosition>
{
    using ParentType = Experimental::PointSourceContext<GlobalPosition>;
public:
    using Scalar = typename GlobalPosition::value_type;

    PointSourceContext() = default;

    PointSourceContext(GlobalPosition position, std::size_t id, Scalar weight,
                       std::size_t elementIndex, std::size_t embeddings)
    : ParentType(std::move(position), embeddings)
    , id_(id), weight_(weight), elementIndex_(elementIndex)
    {}

    //! The identifier of the point source
    std::size_t id() const
    { return id_; }

    //! The quadrature weight times the integration element
    Scalar weight() const
    { return weight_; }

    //! The index of the element containing the point
    std::size_t elementIndex() const
    { return elementIndex_; }

private:
    std::size_t id_ = 0;
    Scalar weight_ = 0.0;
    std::size_t elementIndex_ = 0;
};

} // end namespace Dumux::EmbeddedCoupling

#endif
