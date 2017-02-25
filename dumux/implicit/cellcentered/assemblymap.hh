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
 *        that contribute to the derivative calculation. This is used for
 *        finite-volume schemes with symmetric sparsity pattern in the global matrix.
 */
#ifndef DUMUX_CC_ASSEMBLY_MAP_HH
#define DUMUX_CC_ASSEMBLY_MAP_HH

#include <dune/istl/bcrsmatrix.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup CellCentered
 * \brief A simple version of the assembly map for cellcentered schemes.
 *        This implementation works for schemes in which for a given cell I only
 *        those cells J have to be prepared in whose stencil the cell I appears.
 *        This means that for the flux calculations in the cells J (in order to compute
 *        the derivatives with respect to cell I), we do not need data on any additional cells J
 *        to compute these fluxes. The same holds for scvfs in the cells J, i.e. we need only those
 *        scvfs in the cells J in which the cell I is in the stencil.
 */
template<class TypeTag>
class CCSimpleAssemblyMap
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using IndexType = typename GridView::IndexSet::IndexType;

    struct DataJ
    {
        IndexType globalJ;
        std::vector<IndexType> scvfsJ;
        // A list of additional scvfs is needed for compatibility
        // reasons with more complex assembly maps (see mpfa)
        std::vector<IndexType> additionalScvfs;
    };

    using Map = std::vector<std::vector<DataJ>>;

public:

    /*!
     * \brief Initialize the AssemblyMap object.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(const Problem& problem)
    {
        map_.clear();
        map_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            // We are looking for the elements I, for which this element J is in the flux stencil
            auto globalJ = problem.elementMapper().index(element);

            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element);

            // obtain the data of J in elements I
            std::vector<std::pair<IndexType, DataJ>> dataJForI;

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                FluxVariables fluxVars;
                const auto& stencil = fluxVars.computeStencil(problem, element, fvGeometry, scvf);

                // insert our index in the neighbor stencils of the elements in the flux stencil
                for (auto globalI : stencil)
                {
                    if (globalI == globalJ)
                        continue;

                    auto it = std::find_if(dataJForI.begin(),
                                           dataJForI.end(),
                                           [globalI](const auto& pair) { return pair.first == globalI; });

                    if (it != dataJForI.end())
                        it->second.scvfsJ.push_back(scvf.index());
                    else
                        dataJForI.emplace_back(std::make_pair(globalI, DataJ({globalJ,
                                                                              std::vector<IndexType>({scvf.index()}),
                                                                              std::vector<IndexType>()})));
                }
            }

            for (auto&& pair : dataJForI)
                map_[pair.first].emplace_back(std::move(pair.second));
        }
    }

    const std::vector<DataJ>& operator[] (const IndexType globalI) const
    { return map_[globalI]; }

private:
    Map map_;
};

} // end namespace Dumux

#endif
