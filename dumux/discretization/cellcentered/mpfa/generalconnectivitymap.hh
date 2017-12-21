// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation
 */
#ifndef DUMUX_CC_MPFA_GENERAL_CONNECTIVITY_MAP_HH
#define DUMUX_CC_MPFA_GENERAL_CONNECTIVITY_MAP_HH

#include <vector>
#include <utility>
#include <algorithm>
#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxstencil.hh>

namespace Dumux
{

/*!
 * \ingroup CellCentered
 * \brief General version of the assembly map for cellcentered schemes. To each
 *        cell I we store a list of cells J that are needed to compute the fluxes
 *        in these cells J that depend on cell I. Furthermore, we store for each cell J
 *        a list of scvfs in which cell I is in the stencil, as well as additional scvfs
 *        that are also required to set up the transmissibilities.
 */
template<class TypeTag>
class CCMpfaGeneralConnectivityMap
{
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxStencil = Dumux::FluxStencil<TypeTag>;

    // To each cell "globalI" there will be a list of "globalJ", in which globalI is part
    // of the stencil. We save the scvfs over which fluxes depend on globalI and a list of
    // additional scvfs which are needed temporarily to set up the transmissibilities of the scvfsJ
    struct DataJ
    {
        IndexType globalJ;
        std::vector<IndexType> scvfsJ;
        std::vector<IndexType> additionalScvfs;
    };

    using Map = std::vector<std::vector<DataJ>>;

public:

    /*!
     * \brief Initialize the ConnectivityMap object.
     *
     * \param fvGridGeometry The grid's finite volume geometry.
     */
    void update(const FVGridGeometry& fvGridGeometry)
    {
        map_.resize(fvGridGeometry.gridView().size(0));
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            // We are looking for the elements I, for which this element J is in the flux stencil
            auto globalJ = fvGridGeometry.elementMapper().index(element);

            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            // obtain the data of J in elements I
            std::vector<std::pair<IndexType, std::vector<DataJ>>> dataJForI;

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto& stencil = FluxStencil::stencil(element, fvGeometry, scvf);

                // insert our index in the neighbor stencils of the elements in the flux stencil
                for (auto globalI : stencil)
                {
                    if (globalI == globalJ)
                        continue;

                    auto it = std::find_if(dataJForI.begin(),
                                           dataJForI.end(),
                                           [globalI](const auto& pair) { return pair.first == globalI; });

                    if (it != dataJForI.end())
                    {
                        // get the data J which corresponds to the actual global J
                        // This will be the first entry, as we do so in the else statement (see below)
                        auto& globalJDataJ = it->second[0];

                        // insert actual scvf index in the list of scvfs which couple I and J
                        globalJDataJ.scvfsJ.push_back(scvf.index());

                        // Also, all scvfs connected to a vertex together with the actual scvf
                        // land in the list of additional scvfs. Of that list we will delete those
                        // that are already in the list of scvfsJ later...
                        const auto scvfVectorAtVertex = MpfaHelper::getScvFacesAtVertex(scvf.vertexIndex(), element, fvGeometry);
                        std::vector<IndexType> scvfIndicesAtVertex(scvfVectorAtVertex.size());
                        for (std::size_t i = 0; i < scvfVectorAtVertex.size(); ++i)
                            scvfIndicesAtVertex[i] = scvfVectorAtVertex[i]->index();
                        globalJDataJ.additionalScvfs.insert(globalJDataJ.additionalScvfs.end(),
                                                            scvfIndicesAtVertex.begin(),
                                                            scvfIndicesAtVertex.end());

                        // all the other dofs in the stencil have to appear as "globalJ" to globalI as well
                        for (auto globalJ2 : stencil)
                        {
                            if (globalJ2 == globalJ || globalJ2 == globalI)
                                continue;

                            auto it2 = std::find_if(it->second.begin(),
                                                    it->second.end(),
                                                    [globalJ2](const auto& dataJ) { return dataJ.globalJ == globalJ2; });

                            // if entry for globalJ2 does not exist yet, add globalJ2 to the J-data of globalI
                            // with an empty set of scvfs over which I and J are coupled (i.e. they aren't coupled)
                            if (it2 == it->second.end())
                                it->second.push_back(DataJ({globalJ2, std::vector<IndexType>()}));
                        }
                    }
                    else
                    {
                        // No DataJ for globalI exists yet. Make it and insert data on the actual
                        // global J as first entry in the vector of DataJs belonging to globalI
                        dataJForI.emplace_back(std::make_pair(globalI,
                                                              std::vector<DataJ>({DataJ({globalJ, std::vector<IndexType>({scvf.index()})})})));

                        // Also, all scvfs connected to a vertex together with the actual scvf
                        // land in the list of additional scvfs. Of that list we will delete those
                        // that are already in the list of scvfsJ later...
                        const auto scvfVectorAtVertex = MpfaHelper::getScvFacesAtVertex(scvf.vertexIndex(), element, fvGeometry);
                        std::vector<IndexType> scvfIndicesAtVertex(scvfVectorAtVertex.size());
                        for (unsigned int i = 0; i < scvfVectorAtVertex.size(); ++i)
                            scvfIndicesAtVertex[i] = scvfVectorAtVertex[i]->index();
                        dataJForI.back().second[0].additionalScvfs.insert(dataJForI.back().second[0].additionalScvfs.end(),
                                                                          scvfIndicesAtVertex.begin(),
                                                                          scvfIndicesAtVertex.end());

                        // all the other dofs in the stencil will be "globalJ" to globalI as well
                        for (auto globalJ2 : stencil)
                            if (globalJ2 != globalJ && globalJ2 != globalI)
                                dataJForI.back().second.push_back(DataJ({globalJ2, std::vector<IndexType>()}));
                    }
                }
            }

            // Insert the data into the global map
            for (auto&& pair : dataJForI)
            {
                // obtain the corresponding entry in the map
                auto& dataJVector = map_[pair.first];
                for (auto&& dataJ : pair.second)
                {
                    // delete those additionalScvfs indices that are already in the list of scvfs
                    dataJ.additionalScvfs.erase(std::remove_if(dataJ.additionalScvfs.begin(),
                                                               dataJ.additionalScvfs.end(),
                                                               [&dataJ] (const auto& idx)
                                                               { return MpfaHelper::VectorContainsValue(dataJ.scvfsJ, idx); }),
                                                dataJ.additionalScvfs.end());

                    // if entry for j exists in the map already add scvf and additional scvf indices, create otherwise
                    auto it = std::find_if(dataJVector.begin(),
                                           dataJVector.end(),
                                           [&dataJ](const auto& dataJofMap) { return dataJofMap.globalJ == dataJ.globalJ; });

                    if (it != dataJVector.end())
                    {
                        it->scvfsJ.insert(it->scvfsJ.end(), dataJ.scvfsJ.begin(), dataJ.scvfsJ.end());
                        it->additionalScvfs.insert(it->additionalScvfs.end(), dataJ.additionalScvfs.begin(), dataJ.additionalScvfs.end());
                    }
                    else
                        dataJVector.emplace_back(std::move(dataJ));
                }
            }
        }
    }

    const std::vector<DataJ>& operator[] (const IndexType globalI) const
    { return map_[globalI]; }

private:
    Map map_;
};
} // end namespace Dumux

#endif
