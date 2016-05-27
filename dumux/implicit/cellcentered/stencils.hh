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
 * \brief Implements the notion of stencils for cell-centered models
 */
#ifndef DUMUX_CC_STENCILS_HH
#define DUMUX_CC_STENCILS_HH

#include <set>
#include <dumux/implicit/cellcentered/properties.hh>

namespace Dumux
{

/*!
 * \brief Element-related stencils
 */
template<class TypeTag>
class CCElementStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;
public:
    void update(const Problem& problem, const Element& element)
    {
        elementStencil_.clear();

        const auto& fvGeometry = problem.model().fvGeometries(element);
        // loop over sub control faces
        for (const auto& scvf : fvGeometry.scvfs())
        {
            FluxVariables fluxVars;
            const auto& stencil = fluxVars.computeStencil(problem, scvf);

            elementStencil_.insert(elementStencil_.end(), stencil.begin(), stencil.end());
        }
        // make values in elementstencil unique
        std::sort(elementStencil_.begin(), elementStencil_.end());
        elementStencil_.erase(std::unique(elementStencil_.begin(), elementStencil_.end()), elementStencil_.end());

        auto globalI = problem.elementMapper().index(element);
        neighborStencil_ = elementStencil_;

        //remove the element itself and possible ghost neighbors from the stencil
        auto pred = [&problem, globalI](const int i) -> bool
        {
            if (i == globalI)
                return true;
            if (problem.model().fvGeometries().element(i).partitionType() == Dune::GhostEntity)
                return true;
            return false;
        };
        neighborStencil_.erase(std::remove_if(neighborStencil_.begin(), neighborStencil_.end(), pred), neighborStencil_.end());
    }

    //! The full element stencil (all element this element is interacting with)
    const Stencil& elementStencil() const
    {
        return elementStencil_;
    }

    //! The full element stencil without this element
    const Stencil& neighborStencil() const
    {
        return neighborStencil_;
    }

private:
    Stencil elementStencil_;
    Stencil neighborStencil_;
};

/*!
 * \ingroup CCModel
 * \brief The global stencil container class
 */
template<class TypeTag>
class CCStencilsVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using StencilType = CCElementStencils<TypeTag>;

public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        elementStencils_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            // bind the FvGeometry to the element before using it
            problem.model().fvGeometries_().bind(element);

            auto eIdx = problem.elementMapper().index(element);
            elementStencils_[eIdx].update(problem, element);
        }
    }

    //! overload for elements
    template <class Entity>
    typename std::enable_if<Entity::codimension == 0, const StencilType&>::type
    get(const Entity& entity) const
    {
        return elementStencils_[problemPtr_->elementMapper().index(entity)];
    }

private:
    std::vector<CCElementStencils<TypeTag>> elementStencils_;
    const Problem* problemPtr_;
};

} // end namespace

#endif
