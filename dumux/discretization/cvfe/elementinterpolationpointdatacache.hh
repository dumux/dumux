// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Interpolation point data at the quadrature points of an element
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_ELEMENT_INTERPOLATION_POINT_DATA_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_ELEMENT_INTERPOLATION_POINT_DATA_CACHE_HH

#include <cstddef>
#include <ranges>
#include <vector>

#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux::Detail::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief Interpolation point data at the quadrature points of the sub-control volume faces of an element
 * \note Updates keep the storage of existing entries.
 */
template<class InterpolationPointData>
class ScvfInterpolationPointDataCache
{
public:
    template<class Problem, class ElementDiscretization, class ElementVariables>
    void update(const Problem& problem,
                const typename ElementDiscretization::Element& element,
                const ElementDiscretization& elemDisc,
                const ElementVariables& elemVars)
    {
        // the first update sizes the storage exactly, later updates reuse it
        if (scvfCache_.empty())
        {
            std::size_t numQps = 0;
            for (const auto& scvf : scvfs(elemDisc))
                numQps += std::ranges::size(Dumux::CVFE::quadratureRule(elemDisc, scvf));
            scvfCache_.reserve(numQps);
        }

        qpsOffset_.resize(elemDisc.numScvf());
        std::size_t numQps = 0;
        for (const auto& scvf : scvfs(elemDisc))
        {
            qpsOffset_[scvf.index()] = numQps;
            for (const auto& qpData : Dumux::CVFE::quadratureRule(elemDisc, scvf))
            {
                if (numQps == scvfCache_.size())
                    scvfCache_.emplace_back();
                scvfCache_[numQps++].update(problem, element, elemDisc, elemVars, qpData.ipData());
            }
        }
        scvfCache_.resize(numQps);
    }

    const InterpolationPointData& scvfCache(std::size_t scvfIdx, std::size_t qpIdx) const
    { return scvfCache_[qpsOffset_[scvfIdx] + qpIdx]; }

    InterpolationPointData& scvfCache(std::size_t scvfIdx, std::size_t qpIdx)
    { return scvfCache_[qpsOffset_[scvfIdx] + qpIdx]; }

private:
    std::vector<InterpolationPointData> scvfCache_;
    std::vector<std::size_t> qpsOffset_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Interpolation point data at the quadrature points of an element and of its boundary faces
 * \note Updates keep the storage of existing entries.
 */
template<class InterpolationPointData>
class ElementAndBoundaryInterpolationPointDataCache
{
public:
    template<class Problem, class ElementDiscretization, class ElementVariables>
    void update(const Problem& problem,
                const typename ElementDiscretization::Element& element,
                const ElementDiscretization& elemDisc,
                const ElementVariables& elemVars)
    {
        const auto elemQuadRule = Dumux::CVFE::quadratureRule(elemDisc, element);
        elementCache_.resize(std::ranges::size(elemQuadRule));
        for (const auto& qpData : elemQuadRule)
            elementCache_[qpData.ipData().qpIndex()].update(problem, element, elemDisc, elemVars, qpData.ipData());

        boundaryFaceCache_.resize(std::ranges::size(boundaryFaces(elemDisc)));
        for (const auto& boundaryFace : boundaryFaces(elemDisc))
        {
            auto& bfCache = boundaryFaceCache_[boundaryFace.index()];
            const auto quadRule = Dumux::CVFE::quadratureRule(elemDisc, boundaryFace);
            bfCache.resize(std::ranges::size(quadRule));
            for (const auto& qpData : quadRule)
                bfCache[qpData.ipData().qpIndex()].update(problem, element, elemDisc, elemVars, qpData.ipData());
        }
    }

    const InterpolationPointData& elementCache(std::size_t qpIdx) const
    { return elementCache_[qpIdx]; }

    InterpolationPointData& elementCache(std::size_t qpIdx)
    { return elementCache_[qpIdx]; }

    const InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx) const
    { return boundaryFaceCache_[bfIdx][qpIdx]; }

    InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx)
    { return boundaryFaceCache_[bfIdx][qpIdx]; }

private:
    std::vector<InterpolationPointData> elementCache_;
    std::vector<std::vector<InterpolationPointData>> boundaryFaceCache_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Interpolation point data at the quadrature points of the sub-control volume faces,
 *        the element and the boundary faces of an element
 */
template<class InterpolationPointData>
class HybridInterpolationPointDataCache
{
public:
    template<class Problem, class ElementDiscretization, class ElementVariables>
    void update(const Problem& problem,
                const typename ElementDiscretization::Element& element,
                const ElementDiscretization& elemDisc,
                const ElementVariables& elemVars)
    {
        scvfCache_.update(problem, element, elemDisc, elemVars);
        elementAndBoundaryCache_.update(problem, element, elemDisc, elemVars);
    }

    const InterpolationPointData& scvfCache(std::size_t scvfIdx, std::size_t qpIdx) const
    { return scvfCache_.scvfCache(scvfIdx, qpIdx); }

    InterpolationPointData& scvfCache(std::size_t scvfIdx, std::size_t qpIdx)
    { return scvfCache_.scvfCache(scvfIdx, qpIdx); }

    const InterpolationPointData& elementCache(std::size_t qpIdx) const
    { return elementAndBoundaryCache_.elementCache(qpIdx); }

    InterpolationPointData& elementCache(std::size_t qpIdx)
    { return elementAndBoundaryCache_.elementCache(qpIdx); }

    const InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx) const
    { return elementAndBoundaryCache_.boundaryFaceCache(bfIdx, qpIdx); }

    InterpolationPointData& boundaryFaceCache(std::size_t bfIdx, std::size_t qpIdx)
    { return elementAndBoundaryCache_.boundaryFaceCache(bfIdx, qpIdx); }

private:
    ScvfInterpolationPointDataCache<InterpolationPointData> scvfCache_;
    ElementAndBoundaryInterpolationPointDataCache<InterpolationPointData> elementAndBoundaryCache_;
};

} // end namespace Dumux::Detail::CVFE

#endif
