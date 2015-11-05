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
#include <dumux/common/propertysystem.hh>

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
template <class TypeTag>
class CakeGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    typedef std::shared_ptr<Grid> GridPointer;
    typedef Dumux::LinearSpacer<TypeTag> SpacerDefault;
    typedef typename Dumux::PowerLawSpacer<TypeTag> PowerLawSpacer;

    enum { dim = Grid::dimension,
           dimWorld = Grid::dimensionworld
         };

    enum
    {
            radiusIdx = 0,
            angleIdx = 1,
            zIdx = 2
    };
    typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;
    typedef std::array<unsigned int, dim> CellArray;
    typedef Dune::FieldVector<typename Grid::ctype, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "The CakeGridCreator is not implemented for 1D and 2D.");
        }
        VectorField coordinates;;
        createVectors(coordinates);
        gridPtr() = CakeGridCreator<TypeTag>::createCakeGrid(coordinates);
    }

    static void createVectors(VectorField &coordinates)
    {
        //Set the coordinate vector
        coordinates.resize(1);
        coordinates[0] = {0, 0, 0};
        coordinates.resize(2);
        coordinates[1] = {1, 0, 0};
        coordinates.resize(3);
        coordinates[2] = {0, 1, 0};
        coordinates.resize(4);
        coordinates[3] = {1, 1, 0};
        coordinates.resize(5);
        coordinates[4] = {0, 0, 1};
        coordinates.resize(6);
        coordinates[5] = {1, 0, 1};
        coordinates.resize(7);
        coordinates[6] = {0, 1, 1};
        coordinates.resize(8);
        coordinates[7] = {1, 1, 1};

        // The required parameters
        GlobalPosition polarCoorMin;
        GlobalPosition polarCoorMax;
        CellArray cells;

        try {
        polarCoorMin =
                GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, GlobalPosition, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), PolarCoorMin);
        polarCoorMax =
                GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, GlobalPosition, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), PolarCoorMax);

        std::fill(cells.begin(), cells.end(), 1);
        cells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Cells);
        }
        catch (Dumux::ParameterException &e) {
            DUNE_THROW(Dumux::ParameterException, "Please supply the mandatory parameters: "
                                          << "PolarCoorMin, PolarCoorMax or Cells");
        }

        //Check the input file whether UsePowerLawSpacerRadius is specified, if yes use type PowerLawSpacer else default type
        try {
            if(GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, bool, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), UsePowerLawSpacerRadius)){
                PowerLawSpacer spacerR(polarCoorMin[radiusIdx], polarCoorMax[radiusIdx], cells[radiusIdx]);
            }
            else
                SpacerDefault spacerR(polarCoorMin[radiusIdx], polarCoorMax[radiusIdx], cells[radiusIdx]);
        }
        catch (...)
        {
            SpacerDefault spacerR(polarCoorMin[radiusIdx], polarCoorMax[radiusIdx], cells[radiusIdx]);
        }
        SpacerDefault spacerAngle(polarCoorMin[angleIdx], polarCoorMax[angleIdx], cells[angleIdx]);
        SpacerDefault spacerZ(polarCoorMin[zIdx], polarCoorMax[zIdx], cells[zIdx]);

    }

    static std::shared_ptr<Grid> createCakeGrid(VectorField &coordinates)
    {

        GridFactory gf;
        Dune::FieldVector<Scalar, dim> pos;

        gf.insertVertex(coordinates[0]);
        gf.insertVertex(coordinates[1]);
        gf.insertVertex(coordinates[2]);
        gf.insertVertex(coordinates[3]);
        gf.insertVertex(coordinates[4]);
        gf.insertVertex(coordinates[5]);
        gf.insertVertex(coordinates[6]);
        gf.insertVertex(coordinates[7]);

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
