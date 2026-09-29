// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SurfaceRunoff
 * \brief Channel cross-sections: the geometry a one-dimensional reach conveys and stores over.
 *
 * A section is described as functions of the flow depth \f$h\f$ by the wetted area \f$A(h)\f$,
 * the wetted perimeter \f$P(h)\f$, which determines the hydraulic radius \f$R = A/P\f$, and the
 * conveyance \f$K(h) = A R^{2/3}\f$ of Manning's law.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_CROSSSECTION_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_CROSSSECTION_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux::SurfaceRunoff {

/*!
 * \ingroup SurfaceRunoff
 * \brief \f$A R^{2/3}\f$ from an area and a wetted perimeter, zero on a dry section
 */
template<class Scalar>
Scalar conveyanceFromAreaAndPerimeter(const Scalar area, const Scalar perimeter)
{
    using std::pow;
    if (!(area > 0.0) || !(perimeter > 0.0))
        return 0.0;
    return area*pow(area/perimeter, 2.0/3.0);
}

/*!
 * \ingroup SurfaceRunoff
 * \brief A rectangle of fixed width.
 */
template<class Scalar>
class RectangularSection
{
public:
    explicit RectangularSection(const Scalar width)
    : width_(width)
    {
        if (!(width > 0.0))
            DUNE_THROW(Dumux::ParameterException, "Rectangular section needs a positive width");
    }

    Scalar area(const Scalar h) const { using std::max; return max(Scalar(0.0), h)*width_; }
    Scalar topWidth(const Scalar h) const { return h > 0.0 ? width_ : Scalar(0.0); }

    Scalar wettedPerimeter(const Scalar h) const
    {
        using std::max;
        if (!(h > 0.0))
            return 0.0;
        if (std::isinf(width_))
            return std::numeric_limits<Scalar>::infinity();
        return width_ + 2.0*h;
    }

    Scalar conveyance(const Scalar h) const
    {
        using std::max, std::pow;
        if (!(h > 0.0))
            return 0.0;
        if (std::isinf(width_))
            return std::numeric_limits<Scalar>::infinity();
        return conveyanceFromAreaAndPerimeter(area(h), wettedPerimeter(h));
    }

private:
    Scalar width_;
};

/*!
 * \ingroup SurfaceRunoff
 * \brief A trapezoid that stops widening at its bank.
 *
 * Bottom width \f$b\f$ reaching top width \f$t\f$ at the bank depth \f$d\f$. Above the bank, the
 * section continues with vertical walls, leaving the flow over the floodplain to a surrounding
 * two-dimensional domain. \f$t = b\f$ gives a rectangle and \f$b = 0\f$ a triangle.
 */
template<class Scalar>
class TrapezoidalSection
{
public:
    TrapezoidalSection(const Scalar bottomWidth, const Scalar topWidth, const Scalar bankDepth)
    : bottom_(bottomWidth), top_(topWidth), bankDepth_(bankDepth)
    {
        if (!(bankDepth > 0.0))
            DUNE_THROW(Dumux::ParameterException, "Trapezoidal section needs a positive bank depth");
        if (bottomWidth < 0.0 || topWidth < bottomWidth)
            DUNE_THROW(Dumux::ParameterException, "Trapezoidal section needs 0 <= bottom <= top width");

        // length of one sloping bank per unit depth
        const auto offset = 0.5*(top_ - bottom_)/bankDepth_;
        using std::sqrt;
        bankLengthPerDepth_ = sqrt(1.0 + offset*offset);
    }

    Scalar topWidth(const Scalar h) const
    {
        using std::max, std::min;
        if (!(h > 0.0))
            return 0.0;
        return bottom_ + (top_ - bottom_)*min(h, bankDepth_)/bankDepth_;
    }

    Scalar area(const Scalar h) const
    {
        using std::max;
        if (!(h > 0.0))
            return 0.0;
        if (h <= bankDepth_)
            return 0.5*(bottom_ + topWidth(h))*h;
        return 0.5*(bottom_ + top_)*bankDepth_ + top_*(h - bankDepth_);
    }

    Scalar wettedPerimeter(const Scalar h) const
    {
        using std::max;
        if (!(h > 0.0))
            return 0.0;
        if (h <= bankDepth_)
            return bottom_ + 2.0*h*bankLengthPerDepth_;
        return bottom_ + 2.0*bankDepth_*bankLengthPerDepth_ + 2.0*(h - bankDepth_);
    }

    Scalar conveyance(const Scalar h) const
    { return conveyanceFromAreaAndPerimeter(area(h), wettedPerimeter(h)); }

    Scalar bankDepth() const { return bankDepth_; }
    Scalar bankFullArea() const { return area(bankDepth_); }

private:
    Scalar bottom_, top_, bankDepth_, bankLengthPerDepth_;
};

/*!
 * \ingroup SurfaceRunoff
 * \brief A section given as a table of depth against area and conveyance.
 *
 * The conveyance is taken from data, since the conveyance of a compound section is summed over
 * subdivisions with their own roughness and cannot be recovered from a single area and perimeter.
 *
 * Between knots, the area and the conveyance are interpolated linearly. Above the last knot, the
 * area grows with the top width of the last segment and the conveyance with the wide-channel
 * scaling \f$K \propto A^{5/3}\f$.
 */
template<class Scalar>
class TabulatedSection
{
public:
    TabulatedSection(std::vector<Scalar> depth, std::vector<Scalar> area, std::vector<Scalar> conveyance)
    : depth_(std::move(depth)), area_(std::move(area)), conveyance_(std::move(conveyance))
    {
        if (depth_.size() < 2 || depth_.size() != area_.size() || depth_.size() != conveyance_.size())
            DUNE_THROW(Dumux::ParameterException, "Tabulated section needs at least two consistent rows");
        if (!std::is_sorted(depth_.begin(), depth_.end()))
            DUNE_THROW(Dumux::ParameterException, "Tabulated section needs increasing depths");
        if (depth_.front() != 0.0)
            DUNE_THROW(Dumux::ParameterException, "Tabulated section has to start at zero depth");
        if (!std::is_sorted(area_.begin(), area_.end()))
            DUNE_THROW(Dumux::ParameterException, "Tabulated section needs a non-decreasing area");

        // a conveyance decreasing with depth would make kinematic waves travel against the flow
        if (!std::is_sorted(conveyance_.begin(), conveyance_.end()))
            DUNE_THROW(Dumux::ParameterException, "Tabulated section has a decreasing conveyance");
    }

    Scalar area(const Scalar h) const { return interpolate_(h, area_); }
    Scalar conveyance(const Scalar h) const { return interpolate_(h, conveyance_); }

    //! Derivative of the area with respect to the depth
    Scalar topWidth(const Scalar h) const
    {
        if (!(h > 0.0))
            return 0.0;
        const auto i = segment_(h);
        return (area_[i+1] - area_[i])/(depth_[i+1] - depth_[i]);
    }

    Scalar wettedPerimeter(const Scalar h) const
    {
        using std::pow;
        const auto a = area(h), k = conveyance(h);
        if (!(a > 0.0) || !(k > 0.0))
            return 0.0;
        return a/pow(k/a, 1.5); // invert K = A (A/P)^(2/3)
    }

private:
    std::size_t segment_(const Scalar h) const
    {
        const auto it = std::upper_bound(depth_.begin(), depth_.end(), h);
        if (it == depth_.begin())
            return 0;
        const auto i = std::size_t(std::distance(depth_.begin(), it)) - 1;
        return std::min(i, depth_.size() - 2);
    }

    Scalar interpolate_(const Scalar h, const std::vector<Scalar>& y) const
    {
        using std::pow, std::max;
        if (!(h > 0.0))
            return 0.0;

        if (h > depth_.back())
        {
            const auto n = depth_.size() - 1;
            const auto width = (area_[n] - area_[n-1])/(depth_[n] - depth_[n-1]);
            const auto areaAbove = area_[n] + width*(h - depth_.back());
            if (&y == &area_)
                return areaAbove;
            return area_[n] > 0.0 ? conveyance_[n]*pow(areaAbove/area_[n], 5.0/3.0) : Scalar(0.0);
        }

        const auto i = segment_(h);
        const auto w = (h - depth_[i])/(depth_[i+1] - depth_[i]);
        return y[i] + w*(y[i+1] - y[i]);
    }

    std::vector<Scalar> depth_, area_, conveyance_;
};

} // end namespace Dumux::SurfaceRunoff

#endif
