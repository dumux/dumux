// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_VARIABLECLASS_HH
#define DUMUX_VARIABLECLASS_HH

#include "decoupledproperties.hh"

// for  parallelization
#include <dumux/parallel/elementhandles.hh>

/**
 * @file
 * @brief  Base class holding the variables for sequential models.
 * @author Markus Wolff
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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dim>::Entity Vertex;

    typedef typename SolutionTypes::VertexMapper VertexMapper;
    typedef typename SolutionTypes::ElementMapper ElementMapper;

public:
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;//!<type for vector of scalars
    typedef typename std::vector <CellData> CellDataVector;

private:
    const GridView& gridView_;
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;
    const int codim_;
    CellDataVector cellDataVector_;
protected:
    ScalarSolutionType primaryVariablesVector_[numEq];

public:
    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(const GridView& gridView, int codim = 0) :
        gridView_(gridView), elementMapper_(gridView), vertexMapper_(gridView), codim_(codim)
    {
        for (int i = 0; i < numEq; i++)
        {
            primaryVariablesVector_[i].resize(gridView.size(codim));
            primaryVariablesVector_[i] = 0;
        }
        cellDataVector_.resize(gridView.size(codim));
    }

    // serialization methods
    //! Function needed for restart option.
    template<class Restarter>
    void serialize(Restarter &res)
    {
        res.template serializeEntities<0> (*this, gridView_);
    }

    //! Function needed for restart option.
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities<0> (*this, gridView_);
    }

    //! Function needed for restart option.
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = elementMapper_.map(element);
        outstream << primaryVariablesVector_[0][globalIdx][0];
        for (int i = 1; i < numEq; i++)
        {
            outstream <<"   "<< primaryVariablesVector_[i][globalIdx][0];
        }
    }

    //! Function needed for restart option.
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = elementMapper_.map(element);
        instream >> primaryVariablesVector_[0][globalIdx][0];
        for (int i = 1; i < numEq; i++)
        {
            instream >> primaryVariablesVector_[i][globalIdx][0];
        }
    }

    //! Resizes decoupled variable vectors
    /*! Method that change the size of the vectors for h-adaptive simulations.
     *
     *\param size Size of the current (refined and coarsened) grid
     */
    void adaptVariableSize(int size)
        {
        for (int i = 0; i < numEq; i++)
        {
            primaryVariablesVector_[i].resize(size);
            primaryVariablesVector_[i] = 0;
        }
            cellDataVector_.resize(size);
        }

//    void communicatePressure()
//    {
//#if HAVE_MPI
//        ElementHandleAssign<typename ScalarSolutionType::block_type, ScalarSolutionType, ElementMapper> elementHandle(primaryVariables_, elementMapper_);
//        gridView_.communicate(elementHandle,
//                              Dune::InteriorBorder_All_Interface,
//                              Dune::ForwardCommunication);
//#endif
//    }

    //! Return pressure vector
    const ScalarSolutionType& primaryVariablesGlobal(int eqIdx) const
    {
        return primaryVariablesVector_[eqIdx];
    }

    ScalarSolutionType& primaryVariablesGlobal(int eqIdx)
    {
        return primaryVariablesVector_[eqIdx];
    }

    //! Return saturation vector
    CellDataVector& cellDataGlobal()
    {
        return cellDataVector_;
    }

    const CellDataVector& cellDataGlobal() const
    {
        return cellDataVector_;
    }

    //! Return saturation vector
    CellData& cellData(int idx)
    {
        return cellDataVector_[idx];
    }

    const CellData& cellData(int idx) const
    {
        return cellDataVector_[idx];
    }

    //! Get index of element (codim 0 entity)
    /*! Get index of element (codim 0 entity).
     * @param element codim 0 entity
     * \return element index
     */
    int index(const Element& element) const
    {
        return elementMapper_.map(element);
    }

    //! Get index of vertex (codim dim entity)
    /*! Get index of vertex (codim dim entity).
     * @param vertex codim dim entity
     * \return vertex index
     */
    int index(const Vertex& vertex) const
    {
        return vertexMapper_.map(vertex);
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
