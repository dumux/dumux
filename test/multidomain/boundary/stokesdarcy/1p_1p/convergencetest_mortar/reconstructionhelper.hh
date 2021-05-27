// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_ONEP_MORTAR_RECONSTRUCTION_HELPER_HH
#define DUMUX_ONEP_MORTAR_RECONSTRUCTION_HELPER_HH

#include <unordered_map>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/functionspacebasis.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/multidomain/glue.hh>

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
class MortarReconstructionHelper
{
    template<class GridGeometry>
    using IndexType = typename IndexTraits< typename GridGeometry::GridView >::GridIndex;

public:
    template<class GridGeometry>
    using ElementScvfIndexMap = std::unordered_map< IndexType<GridGeometry>,
                                                    std::vector<IndexType<GridGeometry>> >;

    /*!
     * \brief TODO doc me.
     */
    template<class SubDomainGridGeometry, class MortarGridGeometry, class GlueType>
    static void findCoupledElements(const SubDomainGridGeometry& sdGridGeometry,
                                    const MortarGridGeometry& mortarGridGeometry,
                                    const GlueType& glue,
                                    std::vector<std::size_t>& coupledSdElements,
                                    std::vector<std::size_t>& coupledMortarElements)
    {
        // store which elements are coupled
        coupledSdElements.clear();
        coupledMortarElements.clear();
        coupledSdElements.reserve(glue.size());
        coupledMortarElements.reserve(glue.size());
        for (const auto& is : intersections(glue))
        {
            coupledSdElements.push_back( sdGridGeometry.elementMapper().index(is.domainEntity(0)) );
            for (unsigned int nIdx = 0; nIdx < is.numTargetNeighbors(); ++nIdx)
                coupledMortarElements.push_back( mortarGridGeometry.elementMapper().index(is.targetEntity(nIdx)) );
        }
    }

    /*!
     * \brief TODO doc me.
     */
    template<class SDGridGeometry, class MortarGridGeometry>
    static ElementScvfIndexMap<SDGridGeometry>
    findCoupledScvfs(const SDGridGeometry& sdGridGeometry,
                     const MortarGridGeometry& mortarGridGeometry,
                     const std::vector<std::size_t>& coupledSdElements,
                     const std::vector<std::size_t>& coupledMortarElements)
    {
        using ctype = typename SDGridGeometry::GridView::ctype;
        static constexpr ctype eps = 1.5e-7;

        ElementScvfIndexMap<SDGridGeometry> map;
        for (auto eIdxMortar : coupledMortarElements)
        {
            const auto& mortarElement = mortarGridGeometry.element(eIdxMortar);
            const auto& mortarElementGeometry = mortarElement.geometry();
            const auto& mortarElementCenter = mortarElementGeometry.center();
            const auto& mortarElementCorner = mortarElementGeometry.corner(1);
            auto mortarElementEdge = mortarElementCorner - mortarElementGeometry.corner(0);
            mortarElementEdge /= mortarElementEdge.two_norm();

            for (auto eIdxSubDomain : coupledSdElements)
            {
                const auto sdElement = sdGridGeometry.element(eIdxSubDomain);
                auto fvGeometry = localView(sdGridGeometry);
                fvGeometry.bind(sdElement);

                using std::abs;
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // check if face and mortar element are parallel
                    if ( abs(scvf.unitOuterNormal()*mortarElementEdge) > eps )
                        continue;

                    // check if they are in the same plane
                    auto d = scvf.center() - mortarElementCenter;
                    if (d.two_norm() < eps)
                        d = scvf.center() - mortarElementCorner;
                    d /= d.two_norm();

                    if ( abs(scvf.unitOuterNormal()*d) < eps )
                        map[eIdxSubDomain].push_back( scvf.index() );
                }
            }
        }

        // remove duplicates
        for (auto& entry : map)
        {
            auto& v = entry.second;
            std::sort(v.begin(), v.end());
            v.erase( std::unique(v.begin(), v.end()), v.end() );
        }

        return map;
    }

    /*!
     * \brief TODO doc me.
     */
    template<class SDGridGeometry, class MortarGridGeometry>
    static ElementScvfIndexMap<SDGridGeometry>
    findCoupledScvfs(const SDGridGeometry& sdGridGeometry,
                     const MortarGridGeometry& mortarGridGeometry)
    {
        const auto glue = makeGlue(sdGridGeometry, mortarGridGeometry);

        std::vector<std::size_t> coupledSdElements;
        std::vector<std::size_t> coupledMortarElements;
        findCoupledElements(sdGridGeometry, mortarGridGeometry, glue,
                            coupledSdElements, coupledMortarElements);

        return findCoupledScvfs(sdGridGeometry, mortarGridGeometry,
                                coupledSdElements, coupledMortarElements);
    }
};

} // end namespace Dumux

#endif
