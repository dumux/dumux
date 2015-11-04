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
 *
 * \brief Provides a grid creator which a regular grid made of
 *        quadrilaterals.
 */
#ifndef DUMUX_CAKE_GRID_CREATOR_HH
#define DUMUX_CAKE_GRID_CREATOR_HH

//#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/basicproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

template <class TypeTag>
class LinearSpacer
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    LinearSpacer(Scalar offSet, Scalar length, Scalar noOfPoints)
    : offSet_(offSet), length_(length), noOfPoints_(noOfPoints)
    {}

    void createSpacing() const
    {
        DUNE_THROW(Dune::NotImplemented, "Linear spacer not implemented yet.");
    }

private:
    Scalar offSet_;
    Scalar length_;
    Scalar noOfPoints_;
};

template <class TypeTag>
class PowerLawSpacer
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    PowerLawSpacer(Scalar offSet, Scalar length, Scalar noOfPoints)
    : offSet_(offSet), length_(length), noOfPoints_(noOfPoints)
    {
        power_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PowerLawSpacer, Power);
    }

    void createSpacing() const
    {
        DUNE_THROW(Dune::NotImplemented, "Power law spacer not implemented yet.");
    }

private:
    Scalar power_;
    Scalar offSet_;
    Scalar length_;
    Scalar noOfPoints_;
};

/*!
 * \brief Provides a grid creator which a regular grid made of
 *        quadrilaterals.
 *
 * A quadirlateral is a line segment in 1D, a rectangle in 2D and a
 * cube in 3D.
 */
template <class TypeTag,
          class SpacerRadius = Dumux::LinearSpacer<TypeTag>>
class CakeGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    typedef std::shared_ptr<Grid> GridPointer;
    typedef Dumux::LinearSpacer<TypeTag> SpacerDefault;

    enum { dim = Grid::dimension,
           dimWorld = Grid::dimensionworld
         };

    enum
    {
            radiusIdx = 0,
            angleIdx = 1,
            zIdx = 2
    };

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        if (dim == 2)
        {
            DUNE_THROW(Dune::NotImplemented, "The CakeGridCreator is not implemented for 1D and 2D.");
        }
        gridPtr() = CakeGridCreator<TypeTag, SpacerRadius>::createCakeGrid();
    }

    static std::shared_ptr<Grid> createCakeGrid()
    {
        // The required parameters
        typedef Dune::FieldVector<typename Grid::ctype, dimWorld> GlobalPosition;
        const GlobalPosition polarCoorMin =
                GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, GlobalPosition, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), PolarCoorMin);
        const GlobalPosition polarCoorMax =
                GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, GlobalPosition, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), PolarCoorMax);

        // The optional parameters (they have a default)
        typedef std::array<unsigned int, dim> CellArray;
        CellArray cells;
        std::fill(cells.begin(), cells.end(), 1);
        try { cells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Cells); }
        catch (Dumux::ParameterException &e) { }

        SpacerRadius spacerR(polarCoorMin[radiusIdx], polarCoorMax[radiusIdx], cells[radiusIdx]);
        SpacerDefault spacerAngle(polarCoorMin[angleIdx], polarCoorMax[angleIdx], cells[angleIdx]);
        SpacerDefault spacerZ(polarCoorMin[zIdx], polarCoorMax[zIdx], cells[zIdx]);
        GridFactory gf;
        Dune::FieldVector<Scalar, dim> pos;

        gf.insertVertex({0, 0, 0});
        gf.insertVertex({1, 0, 0});
        gf.insertVertex({0, 1, 0});
        gf.insertVertex({1, 1, 0});
        gf.insertVertex({0, 0, 1});
        gf.insertVertex({1, 0, 1});
        gf.insertVertex({0, 1, 1});
        gf.insertVertex({1, 1, 1});

        Dune::GeometryType type;
        type.makeHexahedron();
        std::vector<unsigned int> vid({0, 1, 2, 3, 4, 5, 6, 7});
        gf.insertElement(type, vid);
        return std::shared_ptr<Grid>(gf.createGrid());
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        gridPtr()->loadBalance();
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    static GridPointer &gridPtr()
    {
        static GridPointer cakeGrid;
        return cakeGrid;
    }
};

}

#endif
