// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Quadrature rules over sub-control volumes and sub-control volume faces
 */
#ifndef DUMUX_DISCRETIZATION_QUADRATURE_RULES_HH
#define DUMUX_DISCRETIZATION_QUADRATURE_RULES_HH

#include <array>
#include <ranges>
#include <dune/common/std/type_traits.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Quadrature point data with weight and interpolation point info
 */
template<class IpData>
class QuadraturePointData
{
    using Scalar = typename IpData::GlobalPosition::value_type;

public:
    QuadraturePointData(Scalar w, IpData ip)
    : weight_(w), ipData_(std::move(ip))
    {}

    const Scalar& weight() const { return weight_; }
    const IpData& ipData() const { return ipData_; }

private:
    Scalar weight_;
    IpData ipData_;
};

namespace QuadratureRules {

/*!
 * \ingroup Discretization
 * \brief Midpoint quadrature rule that uses scv/scvf centers
 */
struct MidpointQuadrature : public Utility::Tag<MidpointQuadrature> {};

/*!
 * \ingroup Discretization
 * \brief Dune-based quadrature rule with specified order
 */
template<int order>
struct DuneQuadrature : public Utility::Tag<DuneQuadrature<order>>
{
    static constexpr int quadratureOrder = order;
};

} // end namespace QuadratureRules

namespace CVFE {

/*!
 * \ingroup Discretization
 * \brief Quadrature rule traits for discretization schemes
 */
template<class GridView,
         class ScvRule = QuadratureRules::MidpointQuadrature,
         class ScvfRule = QuadratureRules::MidpointQuadrature,
         class ElementRule = QuadratureRules::MidpointQuadrature,
         class IntersectionRule = QuadratureRules::MidpointQuadrature>
struct DefaultQuadratureTraits
{
    using ScvQuadratureRule = ScvRule;
    using ScvfQuadratureRule = ScvfRule;
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
};

namespace Detail {

template<class FVElementGeometry>
using DefinesScvQuadratureRuleType = typename FVElementGeometry::ScvQuadratureRule;

template<class FVElementGeometry>
using DefinesScvfQuadratureRuleType = typename FVElementGeometry::ScvfQuadratureRule;

template<class FVElementGeometry>
using DefinesElementQuadratureRuleType = typename FVElementGeometry::ElementQuadratureRule;

template<class FVElementGeometry>
using DefinesIntersectionQuadratureRuleType = typename FVElementGeometry::IntersectionQuadratureRule;

//! Get ScvQuadratureRule type (default is MidpointQuadrature)
template<class FVElementGeometry>
using ScvQuadratureRuleType = Dune::Std::detected_or_t<QuadratureRules::MidpointQuadrature,
                                                       DefinesScvQuadratureRuleType,
                                                       FVElementGeometry>;

//! Get ScvfQuadratureRule type (default is MidpointQuadrature)
template<class FVElementGeometry>
using ScvfQuadratureRuleType = Dune::Std::detected_or_t<QuadratureRules::MidpointQuadrature,
                                                        DefinesScvfQuadratureRuleType,
                                                        FVElementGeometry>;

//! Get ElementQuadratureRule type (default is MidpointQuadrature)
template<class FVElementGeometry>
using ElementQuadratureRuleType = Dune::Std::detected_or_t<QuadratureRules::MidpointQuadrature,
                                                           DefinesElementQuadratureRuleType,
                                                           FVElementGeometry>;

//! Get IntersectionQuadratureRule type (default is MidpointQuadrature)
template<class FVElementGeometry>
using IntersectionQuadratureRuleType = Dune::Std::detected_or_t<QuadratureRules::MidpointQuadrature,
                                                                DefinesIntersectionQuadratureRuleType,
                                                                FVElementGeometry>;

//! Default LocalDofInterpolationPointData type
template<class GridView>
using LocalDofIpData = Dumux::CVFE::LocalDofInterpolationPointData<
        typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate,
        typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate,
        typename IndexTraits<GridView>::LocalIndex
    >;

//! Default InterpolationPointData type
template<class GridView>
using BaseIpData = Dumux::CVFE::InterpolationPointData<
        typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate,
        typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate
    >;

} // end namespace Detail

//! Midpoint quadrature for scv
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::SubControlVolume& scv,
                    QuadratureRules::MidpointQuadrature)
{
    using GridView = typename FVElementGeometry::GridGeometry::GridView;
    using QIpData = IndexedQuadratureInterpolationPointData<Detail::LocalDofIpData<GridView>>;
    const auto& elementGeo = fvGeometry.elementGeometry();
    const auto ipGlobal = scv.center();
    const auto ipLocal = elementGeo.local(ipGlobal);
    const auto localDofIdx = scv.localDofIndex();

    std::array<QuadraturePointData<QIpData>, 1> result{{
        {scv.volume(), QIpData(0, ipLocal, ipGlobal, localDofIdx)}
    }};
    return result;
}

//! Dune quadrature for scv
template<int order, class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::SubControlVolume& scv,
                    QuadratureRules::DuneQuadrature<order>)
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using QIpData = IndexedQuadratureInterpolationPointData<Detail::LocalDofIpData<GridView>>;
    using Extrusion = Extrusion_t<GridGeometry>;

    auto scvGeo = fvGeometry.geometry(scv);
    const auto& elementGeo = fvGeometry.elementGeometry();
    auto quad = Dune::QuadratureRules<typename GridView::ctype, GridView::dimension>::rule(scvGeo.type(), order);
    const auto localDofIdx = scv.localDofIndex();

    return std::views::iota(0u, quad.size())
         | std::views::transform([quad = std::move(quad), scvGeo = std::move(scvGeo), &elementGeo, localDofIdx](const auto idx) {
            const auto& qp = quad[idx];
            const auto ipGlobal = scvGeo.global(qp.position());
            const auto ipLocal = elementGeo.local(ipGlobal);
            return QuadraturePointData<QIpData>{
                qp.weight() * Extrusion::integrationElement(scvGeo, qp.position()),
                QIpData(idx, ipLocal, ipGlobal, localDofIdx)
            };
        });
}

//! Midpoint quadrature for scvf
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::SubControlVolumeFace& scvf,
                    QuadratureRules::MidpointQuadrature)
{
    using GridView = typename FVElementGeometry::GridGeometry::GridView;
    using FaceIpData = FaceInterpolationPointData<Detail::BaseIpData<GridView>, typename IndexTraits<GridView>::LocalIndex>;
    using QFaceIpData = IndexedQuadratureInterpolationPointData<FaceIpData>;
    const auto& elementGeo = fvGeometry.elementGeometry();
    const auto ipGlobal = scvf.ipGlobal();
    const auto ipLocal = elementGeo.local(ipGlobal);

    std::array<QuadraturePointData<QFaceIpData>, 1> result{{
        {scvf.area(), QFaceIpData(0, scvf.unitOuterNormal(), scvf.index(), ipLocal, ipGlobal)}
    }};
    return result;
}

//! Dune quadrature for scvf
template<int order, class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::SubControlVolumeFace& scvf,
                    QuadratureRules::DuneQuadrature<order>)
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using FaceIpData = FaceInterpolationPointData<Detail::BaseIpData<GridView>, typename IndexTraits<GridView>::LocalIndex>;
    using QFaceIpData = IndexedQuadratureInterpolationPointData<FaceIpData>;

    using Extrusion = Extrusion_t<GridGeometry>;

    auto scvfGeo = fvGeometry.geometry(scvf);
    const auto& elementGeo = fvGeometry.elementGeometry();
    auto quad = Dune::QuadratureRules<typename GridView::ctype, GridView::dimension-1>::rule(scvfGeo.type(), order);
    auto normal = scvf.unitOuterNormal();
    const auto scvfIndex = scvf.index();

    return std::views::iota(0u, quad.size())
         | std::views::transform([quad = std::move(quad), scvfGeo = std::move(scvfGeo), &elementGeo, normal = std::move(normal), scvfIndex](const auto idx) {
            const auto& qp = quad[idx];
            const auto ipGlobal = scvfGeo.global(qp.position());
            const auto ipLocal = elementGeo.local(ipGlobal);
            return QuadraturePointData<QFaceIpData>{
                qp.weight() * Extrusion::integrationElement(scvfGeo, qp.position()),
                QFaceIpData(idx, normal, scvfIndex, ipLocal, ipGlobal)
            };
        });
}

//! Midpoint quadrature for element
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::Element& element,
                    QuadratureRules::MidpointQuadrature)
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using QIpData = IndexedQuadratureInterpolationPointData<Detail::BaseIpData<GridView>>;

    const auto& elementGeo = fvGeometry.elementGeometry();
    const auto ipGlobal = elementGeo.center();
    const auto ipLocal = elementGeo.local(ipGlobal);

    std::array<Dumux::QuadraturePointData<QIpData>, 1> result{{
        {elementGeo.volume(), QIpData(0, ipLocal, ipGlobal)}
    }};
    return result;
}

//! Dune quadrature for element
template<int order, class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::Element& element,
                    QuadratureRules::DuneQuadrature<order>)
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using QIpData = IndexedQuadratureInterpolationPointData<Detail::BaseIpData<GridView>>;
    using Extrusion = Extrusion_t<GridGeometry>;

    const auto& elementGeo = fvGeometry.elementGeometry();
    auto quad = Dune::QuadratureRules<typename GridView::ctype, GridView::dimension>::rule(elementGeo.type(), order);

    return std::views::iota(0u, quad.size())
         | std::views::transform([quad = std::move(quad), &elementGeo](const auto idx) {
            const auto& qp = quad[idx];
            const auto ipGlobal = elementGeo.global(qp.position());
            const auto ipLocal = qp.position();
            return Dumux::QuadraturePointData<QIpData>{
                qp.weight() * Extrusion::integrationElement(elementGeo, qp.position()),
                QIpData(idx, ipLocal, ipGlobal)
            };
        });
}

//! Midpoint quadrature for intersection
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::GridGeometry::GridView::Intersection& is,
                    QuadratureRules::MidpointQuadrature)
{
    // per default we don't add any intersection information to the ipData
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using BFlag = Dumux::BoundaryFlag<typename GridView::Grid>;
    using FaceIpData = FEFaceInterpolationPointData<Detail::BaseIpData<GridView>,
                                                    BFlag,
                                                    typename IndexTraits<GridView>::LocalIndex>;
    using QFaceIpData = IndexedQuadratureInterpolationPointData<FaceIpData>;

    const auto& elementGeo = fvGeometry.elementGeometry();
    const auto isGeometry = is.geometry();
    const auto ipGlobal = isGeometry.center();
    const auto ipLocal = elementGeo.local(ipGlobal);

    std::array<Dumux::QuadraturePointData<QFaceIpData>, 1> result{{
        {isGeometry.volume(),
            QFaceIpData(0, is.centerUnitOuterNormal(), BFlag{ is }, is.indexInInside(), ipLocal, ipGlobal)}
    }};
    return result;
}

//! Dune quadrature for intersection
template<int order, class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::GridGeometry::GridView::Intersection& is,
                    QuadratureRules::DuneQuadrature<order>)
{
    // per default we don't add any intersection information to the ipData
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using BFlag = Dumux::BoundaryFlag<typename GridView::Grid>;
    using FaceIpData = FEFaceInterpolationPointData<Detail::BaseIpData<GridView>,
                                                    BFlag,
                                                    typename IndexTraits<GridView>::LocalIndex>;
    using QFaceIpData = IndexedQuadratureInterpolationPointData<FaceIpData>;
    using Extrusion = Extrusion_t<GridGeometry>;

    const auto& elementGeo = fvGeometry.elementGeometry();
    auto isGeometry = is.geometry();
    auto quad = Dune::QuadratureRules<typename GridView::ctype, GridView::dimension-1>::rule(isGeometry.type(), order);
    auto normal = is.centerUnitOuterNormal();
    auto bFlag = BFlag{ is };
    auto index = is.indexInInside();

    return std::views::iota(0u, quad.size())
         | std::views::transform(
            [quad = std::move(quad), isGeometry = std::move(isGeometry), &elementGeo, normal = std::move(normal),
             bFlag = std::move(bFlag), index](const auto idx) {
            const auto& qp = quad[idx];
            const auto ipGlobal = isGeometry.global(qp.position());
            const auto ipLocal = elementGeo.local(ipGlobal);
            return Dumux::QuadraturePointData<QFaceIpData>{
                qp.weight() * Extrusion::integrationElement(isGeometry, qp.position()),
                QFaceIpData(idx, normal, bFlag, index, ipLocal, ipGlobal)
            };
        });
}

//! Generic quadrature rule for scv that uses ScvQuadratureRule type of FVElementGeometry
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::SubControlVolume& scv)
{
    return quadratureRule(fvGeometry, scv, typename Detail::ScvQuadratureRuleType<FVElementGeometry>{});
}

//! Generic quadrature rule for scvf that uses ScvfQuadratureRule type of FVElementGeometry
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::SubControlVolumeFace& scvf)
{
    return quadratureRule(fvGeometry, scvf, typename Detail::ScvfQuadratureRuleType<FVElementGeometry>{});
}

//! Generic quadrature rule for element that uses ElementQuadratureRule type of FVElementGeometry
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::Element& element)
{
    return quadratureRule(fvGeometry, element, typename Detail::ElementQuadratureRuleType<FVElementGeometry>{});
}

//! Generic quadrature rule for intersection that uses IntersectionQuadratureRule type of FVElementGeometry
template<class FVElementGeometry>
auto quadratureRule(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::GridGeometry::GridView::Intersection& intersection)
{
    return quadratureRule(fvGeometry, intersection, typename Detail::IntersectionQuadratureRuleType<FVElementGeometry>{});
}

} // end namespace CVFE

} // end namespace Dumux

#endif
