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
#ifndef DUMUX_CC_ASSEMBLY_MAP_HH
#define DUMUX_CC_ASSEMBLY_MAP_HH

#include <dumux/implicit/properties.hh>


namespace Dumux
{

template<class TypeTag>
class CCAssemblyMap
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using IndexType = typename GridView::IndexSet::IndexType;
    using Map = std::vector< std::vector< std::vector<IndexType> > >;

public:

    /*!
     * \brief Initialize the AssemblyMap object.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(const Problem& problem)
    {
        map_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            // get a local finite volume geometry object that is bindable
            auto fvGeometryJ = localView(problem.model().globalFvGeometry());

            auto globalI = problem.elementMapper().index(element);
            const auto& neighborStencil = problem.model().stencils(element).neighborStencil();

            map_[globalI].reserve(neighborStencil.size());
            for (auto globalJ : neighborStencil)
            {
                const auto& elementJ = fvGeometryJ.globalFvGeometry().element(globalJ);

                // find the flux vars needed for the calculation of the flux into element
                std::vector<IndexType> fluxVarIndices;

                // only non-ghost neighbors (J) have to be considered, derivatives from non-ghost to ghost dofs
                // are assembled when assembling the ghost element (I)
                if (elementJ.partitionType() != Dune::GhostEntity)
                {
                    fvGeometryJ.bindElement(elementJ);
                    for (auto&& scvFaceJ : scvfs(fvGeometryJ))
                    {
                        auto fluxVarsIdx = scvFaceJ.index();

                        // if globalI is in flux var stencil, add to list
                        FluxVariables fluxVars;
                        const auto fluxStencil = fluxVars.computeStencil(problem, elementJ, fvGeometryJ, scvFaceJ);

                        for (auto globalIdx : fluxStencil)
                        {
                            if (globalIdx == globalI)
                            {
                                fluxVarIndices.push_back(fluxVarsIdx);
                                break;
                            }
                        }
                    }
                }
                map_[globalI].emplace_back(std::move(fluxVarIndices));
            }
        }
    }

    //! Some implementations of the assembly map can be solution-dependent and need to be updated
    void update() {}

    const std::vector<std::vector<IndexType>>& operator [] (const IndexType globalI) const
    { return map_[globalI]; }

private:
    Map map_;
};

}

#endif
