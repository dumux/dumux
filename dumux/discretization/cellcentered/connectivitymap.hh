// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCDiscretization
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation. This is used for
 *        finite-volume schemes with symmetric sparsity pattern in the global matrix.
 */
#ifndef DUMUX_CC_CONNECTIVITY_MAP_HH
#define DUMUX_CC_CONNECTIVITY_MAP_HH

#include <vector>
#include <utility>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/fluxstencil.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief A simple version of the connectivity map for cellcentered schemes.
 *        This implementation works for schemes in which for a given cell I only
 *        those cells J have to be prepared in whose stencil the cell I appears.
 *        This means that for the flux calculations in the cells J (in order to compute
 *        the derivatives with respect to cell I), we do not need data on any additional cells J
 *        to compute these fluxes. The same holds for scvfs in the cells J, i.e. we need only those
 *        scvfs in the cells J in which the cell I is in the stencil.
 */
template<class GridGeometry>
class CCSimpleConnectivityMap
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using FluxStencil = Dumux::FluxStencil<FVElementGeometry>;
    static constexpr int maxElemStencilSize = GridGeometry::maxElementStencilSize;

    struct DataJ
    {
        GridIndexType globalJ;
        typename FluxStencil::ScvfStencilIForJ scvfsJ;
        // A list of additional scvfs is needed for compatibility
        // reasons with more complex connectivity maps (see mpfa)
        typename FluxStencil::ScvfStencilIForJ additionalScvfs;
    };

    using Map = std::vector<std::vector<DataJ>>;

public:

    /*!
     * \brief Initialize the ConnectivityMap object.
     *
     * \param gridGeometry The grid's finite volume geometry.
     */
    void update(const GridGeometry& gridGeometry)
    {
        map_.clear();
        map_.resize(gridGeometry.gridView().size(0));

        // container to store for each element J the elements I which have J in their flux stencil
        Dune::ReservedVector<std::pair<GridIndexType, DataJ>, maxElemStencilSize> dataJForI;
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            // We are looking for the elements I, for which this element J is in the flux stencil
            const auto globalJ = gridGeometry.elementMapper().index(element);
            fvGeometry.bindElement(element);

            // obtain the data of J in elements I
            dataJForI.clear();

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto& stencil = FluxStencil::stencil(element, fvGeometry, scvf);

                // insert our index in the neighbor stencils of the elements in the flux stencil
                for (auto globalI : stencil)
                {
                    if (globalI == globalJ)
                        continue;

                    auto it = std::find_if(dataJForI.begin(), dataJForI.end(),
                                           [globalI](const auto& pair) { return pair.first == globalI; });

                    if (it != dataJForI.end())
                        it->second.scvfsJ.push_back(scvf.index());
                    else
                    {
                        if (dataJForI.size() > maxElemStencilSize - 1)
                            DUNE_THROW(Dune::InvalidStateException, "Maximum admissible stencil size (" << maxElemStencilSize-1
                                                                     << ") is surpassed (" << dataJForI.size() << "). "
                                                                     << "Please adjust the GridGeometry traits accordingly!");

                        dataJForI.push_back(std::make_pair(globalI, DataJ({globalJ, {scvf.index()}, {}})));
                    }
                }
            }

            for (auto&& pair : dataJForI)
                map_[pair.first].emplace_back(std::move(pair.second));
        }
    }

    const std::vector<DataJ>& operator[] (const GridIndexType globalI) const
    { return map_[globalI]; }

private:
    Map map_;
};

} // end namespace Dumux

#endif
