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
 * \brief Implements the notion of stencils for vertex-centered models
 */
#ifndef DUMUX_BOX_STENCILS_HH
#define DUMUX_BOX_STENCILS_HH

#include <set>
#include <dumux/implicit/box/properties.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class BoxStencilsVector;

/*!
 * \brief Element-related stencils
 */
template<class TypeTag>
class BoxElementStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;

    static const int dim = GridView::dimension;
public:
    void update(const Problem& problem, const Element& element)
    {
        for(int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
            elementStencil_.push_back(problem.vertexMapper().subIndex(element, vIdxLocal, dim));
    }

    //! The full element stencil (all element this element is interacting with)
    const Stencil& elementStencil() const
    {
        return elementStencil_;
    }

private:
    Stencil elementStencil_;
};

/*!
 * \brief Vertex-related stencils
 */
template<class TypeTag>
class BoxVertexStencils
{
    friend class BoxStencilsVector<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;
public:
    //! The full vertex stencil (all vertices this vertex is interacting with)
    const Stencil& vertexStencil() const
    {
        return vertexStencil_;
    }

    //! The scv indices connected to a vertex
    const Stencil& vertexScvs() const
    {
        return vertexScvs_;
    }

    //! The element indices adjacent to the vertex
    const Stencil& elementIndices() const
    {
        return elementIndices_;
    }

private:
    //! The full vertex stencil (all vertices this vertex is interacting with)
    Stencil& vertexStencil()
    {
        return vertexStencil_;
    }

    Stencil& vertexScvs()
    {
        return vertexScvs_;
    }

    Stencil& elementIndices()
    {
        return elementIndices_;
    }

    Stencil vertexStencil_;
    Stencil vertexScvs_;
    Stencil elementIndices_;
};

/*!
 * \ingroup BoxModel
 * \brief The global stencil container class
 */
template<class TypeTag>
class BoxStencilsVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int dim = GridView::dimension;

public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        elementStencils_.resize(problem.gridView().size(0));
        vertexStencils_.resize(problem.gridView().size(dim));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            elementStencils_[eIdx].update(problem, element);

            // bind the FvGeometry to the element before using it
            problem.model().fvGeometries_().bindElement(element);
            const auto& fvGeometry = problem.model().fvGeometries(element);

            for (const auto& scv : fvGeometry.scvs())
            {
                auto vIdxGlobal = scv.dofIndex();
                vertexStencils_[vIdxGlobal].vertexScvs().push_back(scv.index());
                vertexStencils_[vIdxGlobal].elementIndices().push_back(eIdx);

                for (const auto& scvJ : fvGeometry.scvs())
                {
                    auto vIdxGlobalJ = scvJ.dofIndex();

                    // make sure that vertex j is in the neighbor set
                    // of vertex i and vice-versa
                    vertexStencils_[vIdxGlobal].vertexStencil().push_back(vIdxGlobalJ);
                    vertexStencils_[vIdxGlobalJ].vertexStencil().push_back(vIdxGlobal);
                }
            }
        }

        // The vertex indices in the stencils have to be made unique
        for (auto& vertexStencil : vertexStencils_)
        {
            auto& vertStencil = vertexStencil.vertexStencil();
            std::sort(vertStencil.begin(), vertStencil.end());
            vertStencil.erase(std::unique(vertStencil.begin(), vertStencil.end()), vertStencil.end());
        }
    }

    //! overload for elements
    template <class Entity>
    typename std::enable_if<Entity::codimension == 0, const BoxElementStencils<TypeTag>&>::type
    get(const Entity& entity) const
    {
        return elementStencils_[problemPtr_->elementMapper().index(entity)];
    }

    //! overload for vertices
    template <class Entity>
    typename std::enable_if<Entity::codimension == Entity::dimension, const BoxVertexStencils<TypeTag>&>::type
    get(const Entity& entity) const
    {
        return vertexStencils_[problemPtr_->vertexMapper().index(entity)];
    }

private:
    std::vector<BoxElementStencils<TypeTag>> elementStencils_;
    std::vector<BoxVertexStencils<TypeTag>> vertexStencils_;
    const Problem* problemPtr_;
};

} // end namespace

#endif
