// $Id$
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

//#define HACK_SINTEF_RESPROP

#include <dune/istl/bvector.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include "decoupledproperties.hh"

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld, numPhase = GET_PROP_VALUE(TypeTag, PTAG(
                NumPhases))
    };

    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename SolutionTypes::VertexMapper VertexMapper;
    typedef typename SolutionTypes::ElementMapper ElementMapper;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
public:
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;//!<type for vector of scalars
    typedef typename SolutionTypes::PhaseProperty PhasePropertyType;//!<type for vector of phase properties
    typedef typename SolutionTypes::FluidProperty FluidPropertyType;//!<type for vector of fluid properties
    typedef typename SolutionTypes::PhasePropertyElemFace PhasePropertyElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef typename SolutionTypes::DimVecElemFace DimVecElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars

private:
    const GridView& gridView_;
    const ElementMapper elementMapper_;
    const VertexMapper vertexMapper_;
    const int gridSize_;

    const int codim_;

    ScalarSolutionType pressure_;
    DimVecElemFaceType velocity_;

    PhasePropertyElemFaceType potential_;

public:
    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */

    VariableClass(const GridView& gridView, Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        gridView_(gridView), elementMapper_(gridView), vertexMapper_(gridView), gridSize_(gridView_.size(0)), codim_(0)
    {
        initializeGlobalVariables(initialVel);
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(const GridView& gridView, int codim, Dune::FieldVector<Scalar,
            dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        gridView_(gridView), elementMapper_(gridView), vertexMapper_(gridView), gridSize_(gridView_.size(codim)), codim_(codim)
    {
        initializeGlobalVariables(initialVel);
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
        outstream << pressure_[globalIdx];
    }

    //! Function needed for restart option.
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = elementMapper_.map(element);
        instream >> pressure_[globalIdx];
    }

    //! initializes the potential differences stored at the element faces.
    void initializePotentials(Dune::FieldVector<Scalar, dim>& initialPot)
    {
        if (initialPot.two_norm())
        {
            // compute update vector
            ElementIterator eItEnd = gridView_.template end<0> ();
            for (ElementIterator eIt = gridView_.template begin<0> (); eIt != eItEnd; ++eIt)
            {
                // cell index
                int globalIdxI = elementMapper_.map(*eIt);

                // run through all intersections with neighbors and boundary
                IntersectionIterator isItEnd = gridView_.iend(*eIt);
                for (IntersectionIterator isIt = gridView_.ibegin(*eIt); isIt != isItEnd; ++isIt)
                {
                    // local number of facet
                    int indexInInside = isIt->indexInInside();

                    Dune::FieldVector<Scalar, dimWorld> unitOuterNormal = isIt->centerUnitOuterNormal();

                    for (int i = 0; i < numPhase; i++) {potential_[globalIdxI][indexInInside][i] = initialPot * unitOuterNormal;}
                }
            }
        }
        else
        {
            potential_ = Dune::FieldVector<Scalar, numPhase> (0);
        }
        return;
    }

private:
    void initializeGlobalVariables(Dune::FieldVector<Scalar, dim>& initialVel)
    {
        //resize to grid size
        pressure_.resize(gridSize_);
        velocity_.resize(gridSize_);//depends on pressure
        potential_.resize(gridSize_);//depends on pressure

        //initialise variables
        pressure_ = 0;
        velocity_ = initialVel;
        initializePotentials(initialVel);
    }

    //Write saturation and pressure into file
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        if (codim_ == 0)
        {
            ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (this->gridSize());

            *pressure = this->pressure();

            writer.addCellData(pressure, "pressure");
        }
        if (codim_ == dim)
        {
            ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (this->gridSize());

            *pressure = this->pressure();

            writer.addVertexData(pressure, "pressure");
        }

        return;
    }
public:

    //! Return pressure vector
    const ScalarSolutionType& pressure() const
    {
        return pressure_;
    }

    ScalarSolutionType& pressure()
    {
        return pressure_;
    }

    //! Return velocity vector
    const DimVecElemFaceType& velocity() const
    {
        return velocity_;
    }

    DimVecElemFaceType& velocity()
    {
        return velocity_;
    }

    //! Return vector of wetting phase potential gradients
    Dune::FieldVector<Scalar, numPhase>& potential(int Idx1, int Idx2)
    {
        return potential_[Idx1][Idx2];
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

    //!Return the number of data elements
    int gridSize() const
    {
        return gridSize_;
    }

    //!Return gridView
    const GridView& gridView() const
    {
        return gridView_;
    }

    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }

    //! Get pressure
    /*! evaluate pressure at given element
     @param element entity of codim 0
     \return value of pressure
     */
    const Dune::FieldVector<Scalar, 1>& pressElement(const Element& element) const
    {
        return pressure_[elementMapper_.map(element)];
    }

    //! Get velocity at given element face
    /*! evaluate velocity at given location
     @param element entity of codim 0
     @param indexInInside index in reference element
     \return vector of velocity
     */
    Dune::FieldVector<Scalar, dim>& velocityElementFace(const Element& element, const int indexInInside)
    {
        int elemId = elementMapper_.map(element);

        return (velocity_[elemId][indexInInside]);
    }

    const Dune::FieldVector<Scalar, dim>& velocityElementFace(const Element& element, const int indexInInside) const
    {
        int elemId = elementMapper_.map(element);

        return (velocity_[elemId][indexInInside]);
    }
};
}
#endif
