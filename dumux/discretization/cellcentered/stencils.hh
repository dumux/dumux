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
#ifndef DUMUX_DISCRETIZATION_CC_STENCILS_HH
#define DUMUX_DISCRETIZATION_CC_STENCILS_HH

#include <dumux/implicit/cellcentered/properties.hh>

namespace Dumux
{

/*!
 * \brief Element-related stencils for symmetric cc methods.
 */
template<class TypeTag>
class CCSymmetricElementStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;
public:
    //! Update the stencil.
    void update(const Problem& problem,
                const Element& element)
    {
        // restrict the FvGeometry locally and bind to the element
        auto fvGeometry = localView(problem.model().globalFvGeometry());
        fvGeometry.bindElement(element); // for TPFA bind element is enough

        elementStencil_.clear();

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            const auto& stencil = fluxVars.computeStencil(problem, element, fvGeometry, scvf);
            elementStencil_.insert(elementStencil_.end(), stencil.begin(), stencil.end());
        }
        // make values in elementstencil unique
        std::sort(elementStencil_.begin(), elementStencil_.end());
        elementStencil_.erase(std::unique(elementStencil_.begin(), elementStencil_.end()), elementStencil_.end());

        auto globalI = problem.elementMapper().index(element);
        neighborStencil_ = elementStencil_;

        // remove the element itself and possible ghost neighbors from the neighbor stencil
        neighborStencil_.erase(std::remove_if(neighborStencil_.begin(), neighborStencil_.end(),
                                             [globalI](int i){ return (i == globalI); }),
                               neighborStencil_.end());
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
 * \brief The global stencil container class for symmetric cc methods
 */
template<class TypeTag>
class CCSymmetricStencilsVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using StencilType = CCSymmetricElementStencils<TypeTag>;

public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        elementStencils_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
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
    std::vector<StencilType> elementStencils_;
    const Problem* problemPtr_;
};

//! Forward declaration of the global stencil class for nonsymmetric cc methods
template<class TypeTag> class CCNonSymmetricStencilsVector;

/*!
 * \brief Element-related stencils for non-symmetric cc methods
 */
template<class TypeTag>
class CCNonSymmetricElementStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;
public:
    //! Update the stencil. We expect a bound fvGeometry
    void update(const Problem& problem,
                const Element& element,
                CCNonSymmetricStencilsVector<TypeTag>& stencilsVector)
    {
        // the element index of our own here
        auto eIdx = problem.elementMapper().index(element);
        elementStencil_.clear();

        // restrict the FvGeometry locally and bind to the element
        auto fvGeometry = localView(problem.model().globalFvGeometry());
        fvGeometry.bindElement(element);

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            const auto& stencil = fluxVars.computeStencil(problem, element, fvGeometry, scvf);

            // insert stencil into the element stencil
            elementStencil_.insert(elementStencil_.end(), stencil.begin(), stencil.end());

            // insert our index in the neighbor stencils of the elements in the flux stencil
            for (auto globalI : stencil)
            {
                if (globalI != eIdx)
                    stencilsVector[globalI].addNeighbor(eIdx);
            }
        }

        // make values in element stencil unique
        std::sort(elementStencil_.begin(), elementStencil_.end());
        elementStencil_.erase(std::unique(elementStencil_.begin(), elementStencil_.end()), elementStencil_.end());
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

    void addNeighbor(const IndexType nIdx)
    {
        neighborStencil_.push_back(nIdx);
    }

    void makeNeighborStencilUnique()
    {
        // make values in neighbor stencil unique
        std::sort(neighborStencil_.begin(), neighborStencil_.end());
        neighborStencil_.erase(std::unique(neighborStencil_.begin(), neighborStencil_.end()), neighborStencil_.end());
    }

private:
    Stencil elementStencil_;
    Stencil neighborStencil_;
};

/*!
 * \ingroup CCModel
 * \brief The global stencil container class for non-symmetric cc methods
 */
template<class TypeTag>
class CCNonSymmetricStencilsVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using StencilType = CCNonSymmetricElementStencils<TypeTag>;

    // The stencil type of an element has to access bracket operator
    friend StencilType;

public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        elementStencils_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            elementStencils_[eIdx].update(problem, element, *this);
        }

        for (auto&& stencil : elementStencils_)
            stencil.makeNeighborStencilUnique();
    }

    //! overload for elements
    template <class Entity>
    typename std::enable_if<Entity::codimension == 0, const StencilType&>::type
    get(const Entity& entity) const
    {
        return elementStencils_[problemPtr_->elementMapper().index(entity)];
    }

private:
    StencilType& operator[] (const IndexType globalIdx)
    { return elementStencils_[globalIdx]; }

    std::vector<StencilType> elementStencils_;
    const Problem* problemPtr_;
};

} // end namespace

#endif
