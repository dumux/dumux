// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief A quadrature rule using local refinement to approximate partitioned elements
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_LOCAL_REFINEMENT_QUADRATURE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_LOCAL_REFINEMENT_QUADRATURE_HH

#include <array>
#include <vector>
#include <utility>
#include <algorithm>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief A quadrature rule using local refinement to approximate partitioned elements
 *
 * An indicator function provides a partition of the given geometry. The function returns the partition id
 * for a given position. The algorithm then virtually refines the geometry using local refinemt to
 * approximate the partition boundaries. One integration point per partition is added at the centroid
 * of the partition.
 *
 * This can be used, e.g. for low-order integration of the coupling term
 * in the projection scheme over the interface facets.
 * Each facet may be coupled to multiple 1D elements.
 * See Koch (2021) https://arxiv.org/abs/2106.06358 which also describes this algorithm.
 */
template<class Geometry, class IndicatorFunction>
class LocalRefinementSimplexQuadrature
{
    using IndicatorTriple = std::array<std::size_t, 3>;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using LocalPosition = typename Geometry::LocalCoordinate;
    using Corners = std::array<GlobalPosition, 3>;

    class QuadPoint
    {
        LocalPosition localPos_;
        double weight_;
        std::size_t indicator_;

    public:
        QuadPoint(LocalPosition&& localPos, double weight, std::size_t i)
        : localPos_(weight*std::move(localPos))
        , weight_(weight)
        , indicator_(i)
        {}

        void add(LocalPosition&& localPos, double weight)
        {
            localPos_ += weight*localPos;
            weight_ += weight;
        }

        void finalize()
        {
            localPos_ /= weight_;
        }

        const LocalPosition& position() const { return localPos_; }
        double weight() const { return weight_; }
        std::size_t indicator() const { return indicator_; }
    };

    class Triangle
    {
        Corners corners_;
    public:
        Triangle(const Corners& c) : corners_(c) {}
        Triangle(Corners&& c) : corners_(std::move(c)) {}
        const GlobalPosition& operator[](std::size_t i) const { return corners_[i]; }

        GlobalPosition center() const
        {
            GlobalPosition center(0.0);
            for (const auto& c : corners_)
                center += c;
            center /= corners_.size();
            return center;
        }

        auto volume() const
        {
            const auto ab = corners_[1] - corners_[0];
            const auto ac = corners_[2] - corners_[0];
            return crossProduct(ab, ac).two_norm()*0.5;
        }
    };

public:
    LocalRefinementSimplexQuadrature(const Geometry& geo, std::size_t maxLevel, const IndicatorFunction& i)
    : ind_(i)
    , maxLevel_(maxLevel)
    , geometry_(geo)
    {
        static constexpr int dim = Geometry::mydimension;
        static_assert(dim == 2, "Only triangles are supported so far");

        if (geo.corners() != (dim+1))
            DUNE_THROW(Dune::InvalidStateException, "Only simplex geometries are allowed");

        // evaluate indicator
        const auto tri = Triangle{{ geo.corner(0), geo.corner(1), geo.corner(2) }};
        const auto triple = IndicatorTriple{{ ind_(tri[0]), ind_(tri[1]), ind_(tri[2]) }};
        volume_ = tri.volume();

        // add new sub-qps recursively by adaptive refinement
        addSimplex_(tri, 0, triple, triple);

        for (auto& qp : qps_)
            qp.finalize();
    }

    auto begin() const { return qps_.begin(); }
    auto end() const { return qps_.end(); }
    auto size() const { return qps_.size(); }

private:
    void addSimplex_(const Triangle& tri, std::size_t level, const IndicatorTriple& triple, const IndicatorTriple& childTriple)
    {
        // if all indicator are the same just add the triangle
        if (std::all_of(childTriple.begin(), childTriple.end(), [a0=childTriple[0]](auto a){ return (a == a0); }))
        {
            // the volumes need to sum up to 0.5 (triangle reference element volume)
            auto it = std::find_if(qps_.begin(), qps_.end(), [&](const auto& qp){ return (qp.indicator() == childTriple[0]); });
            if (it != qps_.end())
                it->add(geometry_.local(tri.center()), 0.5*tri.volume()/volume_);
            else
                qps_.emplace_back(geometry_.local(tri.center()), 0.5*tri.volume()/volume_, childTriple[0]);
        }

        // if they are different but the maximum level is reached, split volume in three
        else if (level == maxLevel_)
        {
            for (int i = 0; i < 3; ++i)
            {
                // the volumes need to sum up to 0.5 (triangle reference element volume)
                auto it = std::find_if(qps_.begin(), qps_.end(), [&](const auto& qp){ return (qp.indicator() == childTriple[i]); });
                if (it != qps_.end())
                    it->add(geometry_.local(tri[i]), 0.5*tri.volume()/volume_/3.0);
                else
                    qps_.emplace_back(geometry_.local(tri[i]), 0.5*tri.volume()/volume_/3.0, childTriple[i]);
            }
        }

        // otherwise refine and go recursively over all children
        else
        {
            const auto children = refine(tri);
            for (const auto& c : children)
            {
                // make sure to recompute indicator every time since new indices might come up
                const auto grandChildTriple = IndicatorTriple{{ ind_(c[0]), ind_(c[1]), ind_(c[2]) }};
                addSimplex_(c, level+1, triple, grandChildTriple);
            }
        }
    }

    std::array<Triangle, 4> refine(const Triangle& tri) const
    {
        const auto p01 = 0.5*(tri[0] + tri[1]);
        const auto p12 = 0.5*(tri[1] + tri[2]);
        const auto p20 = 0.5*(tri[2] + tri[0]);

        return {{
            {{ tri[0], p01, p20 }},
            {{ tri[1], p01, p12 }},
            {{ tri[2], p12, p20 }},
            {{ p12, p01, p20 }}
        }};
    }

    const IndicatorFunction &ind_;
    std::size_t maxLevel_;
    const Geometry& geometry_;
    double volume_;

    std::vector<QuadPoint> qps_;
};

} // end namespace Dumux

#endif
