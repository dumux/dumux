// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DOXYGEN
#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_GRID_GEOMETRY_DETAIL_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_GRID_GEOMETRY_DETAIL_HH

#include <utility>

namespace Dumux::BoxDfmDetail {

// Add fracture geometries on a bulk grid intersection to storage containers
template<class LocalIndexType = unsigned int,
         class BulkFractureIntersection,
         class BulkGridIndex,
         class GeometryHelper,
         class ScvStorage,
         class ScvfStorage>
auto pushFractureGeometries(const BulkFractureIntersection& is,
                            const BulkGridIndex elementIndex,
                            const GeometryHelper& geometryHelper,
                            ScvStorage& scvStorage,
                            ScvfStorage& scvfStorage) {
    static constexpr int dim = std::decay_t<decltype(is.element)>::Geometry::mydimension;
    static_assert(dim == 2 || dim == 3);

    auto scvLocalIdx = scvStorage.size();
    auto scvfLocalIdx = scvfStorage.size();

    // add fracture scv for each vertex of the intersection
    const auto& refElement = is.referenceElement;
    const auto numCorners = is.intersectionGeometry.corners();
    const auto idxInInside = is.intersection.indexInInside();
    const auto scvOffset = scvStorage.size();

    scvStorage.reserve(scvStorage.size() + numCorners);
    for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
        scvStorage.emplace_back(
            geometryHelper,
            is.intersection,
            is.intersectionGeometry,
            vIdxLocal,
            static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, vIdxLocal, dim)),
            scvLocalIdx++,
            idxInInside,
            elementIndex,
            is.intersectionVertexIndices[vIdxLocal]
        );

    // add fracture scvf for each edge of the intersection in 3d
    unsigned int numScvf = 0;
    if constexpr (dim == 3)
    {
        const auto& faceRefElement = referenceElement(is.intersectionGeometry);
        numScvf = faceRefElement.size(1);
        for (unsigned int edgeIdx = 0; edgeIdx < numScvf; ++edgeIdx)
        {
            // inside/outside scv indices in face local node numbering
            std::vector<LocalIndexType> localScvIndices({
                static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 0, dim-1)),
                static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 1, dim-1))
            });

            // add offset to get the right scv indices
            std::for_each(localScvIndices.begin(),
                          localScvIndices.end(),
                          [scvOffset] (auto& elemLocalIdx) { elemLocalIdx += scvOffset; });

            // add scvf
            scvfStorage.emplace_back(
                geometryHelper,
                is.intersection,
                is.intersectionGeometry,
                edgeIdx,
                scvfLocalIdx++,
                std::move(localScvIndices),
                is.intersection.boundary()
            );
        }
    }

    // dim == 2, intersection is an edge, make 1 scvf
    else
    {
        numScvf = 1;
        // inside/outside scv indices in face local node numbering
        std::vector<LocalIndexType> localScvIndices({0, 1});

        // add offset such that the fracture scvs above are addressed
        std::for_each(localScvIndices.begin(),
                      localScvIndices.end(),
                      [scvOffset] (auto& elemLocalIdx) { elemLocalIdx += scvOffset; });

        // add scvf
        scvfStorage.emplace_back(
            geometryHelper,
            is.intersection,
            is.intersectionGeometry,
            /*idxOnIntersection*/0,
            scvfLocalIdx++,
            std::move(localScvIndices),
            is.intersection.boundary()
        );
    }

    return std::make_pair(numCorners, numScvf);
}

} // end namespace Dumux::BoxDfmDetail

#endif
#endif // DOXYGEN
