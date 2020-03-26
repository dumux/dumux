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
#ifndef DUMUX_VARIABLECLASS_HH
#define DUMUX_VARIABLECLASS_HH

#include "properties.hh"


// for  parallelization
//#include <dumux/parallel/elementhandles.hh>

/**
 * @file
 * @brief  Base class holding the variables for sequential models.
 */

namespace Dumux
{

/*!
 * \ingroup Sequential
 */
//! Base class holding the variables and discretized data for sequential models.
/*!
 * Stores global information and variables that are common for all sequential models and also functions needed to access these variables.
 * Can be directly used for a single phase model.
 *
 * @tparam TypeTag The Type Tag
 *
 */
template<class TypeTag>
class VariableClass
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    enum
    {
        dim = GridView::dimension,
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;

    using VertexMapper = typename SolutionTypes::VertexMapper;
    using ElementMapper = typename SolutionTypes::ElementMapper;

public:
    using CellDataVector = typename std::vector<CellData>;

private:
    const GridView gridView_;
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;
    CellDataVector cellDataVector_;

public:
    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    VariableClass(const GridView& gridView) :
        gridView_(gridView),
        elementMapper_(gridView, Dune::mcmgElementLayout()),
        vertexMapper_(gridView, Dune::mcmgVertexLayout())
    {}

    //! Initializes the variable class
    /*! Method initializes the cellData vector.
     * Should be called from problem init()
     */
    void initialize()
    {
        elementMapper_.update();
        vertexMapper_.update();
        cellDataVector_.resize(gridView_.size(0));
    }

    //! Resizes sequential variable vectors
    /*! Method that change the size of the vectors for h-adaptive simulations.
     *
     *\param size Size of the current (refined and coarsened) grid
     */
    void adaptVariableSize(const int size)
    {
        cellDataVector_.resize(size);
    }

    //! Return the vector holding all cell data
    CellDataVector& cellDataGlobal()
    {
        return cellDataVector_;
    }

    const CellDataVector& cellDataGlobal() const
    {
        assert(cellDataVector_.size() == gridView_.size(0));

        return cellDataVector_;
    }

    //! Return the cell data of a specific cell
    CellData& cellData(const int idx)
    {
        assert(cellDataVector_.size() == gridView_.size(0));

        return cellDataVector_[idx];
    }

    const CellData& cellData(const int idx) const
    {
        assert(cellDataVector_.size() == gridView_.size(0));

        return cellDataVector_[idx];
    }

    //! Get index of element (codim 0 entity)
    /*! Get index of element (codim 0 entity).
     * @param element codim 0 entity
     * \return element index
     */
    int index(const Element& element) const
    {
        return elementMapper_.index(element);
    }

    //! Get index of vertex (codim dim entity)
    /*! Get index of vertex (codim dim entity).
     * @param vertex codim dim entity
     * \return vertex index
     */
    int index(const Vertex& vertex) const
    {
        return vertexMapper_.index(vertex);
    }

    //!Return gridView
    const GridView& gridView() const
    {
        return gridView_;
    }

    //! Return mapper for elements (for adaptive grids)
    ElementMapper& elementMapper()
    {
        return elementMapper_;
    }
    //! Return mapper for elements (for static grids)
    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }
    //! Return mapper for vertices (for adaptive grids)
    VertexMapper& vertexMapper()
    {
        return vertexMapper_;
    }
    //! Return mapper for vertices (for static grids)
    const VertexMapper& vertexMapper() const
    {
        return vertexMapper_;
    }
};
}
#endif
