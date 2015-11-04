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
    LinearSpacer(Scalar noOfPoints, Scalar length, Scalar offSet)
    : noOfPoints_(noOfPoints), length_(length), offSet_(offSet)
    {

    }

    void createSpacing() const
    {
        DUNE_THROW(Dune::NotImplemented, "Linear spacer not implemented yet.");
    }

private:
    Scalar noOfPoints_;
    Scalar length_;
    Scalar offSet_;
};

template <class TypeTag>
class PowerLawSpacer
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    PowerLawSpacer(Scalar noOfPoints, Scalar length, Scalar offSet)
    : noOfPoints_(noOfPoints), length_(length), offSet_(offSet)
    {
        power_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PowerLawSpacer, Power);
    }


    void createSpacing() const
    {
        DUNE_THROW(Dune::NotImplemented, "Power law spacer not implemented yet.");
    }

private:
    Scalar power_;
    Scalar noOfPoints_;
    Scalar length_;
    Scalar offSet_;
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


    enum { dim = Grid::dimension };

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

//        //Get the necessary parameters from the parameter file
//        Scalar outerRadius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, OuterRadius);
//        Scalar wellRadius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, WellRadius);
//        int noOfCellsR = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NoOfCellsR);
//        Scalar angle = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, Angle);
//        int noOfSlice = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NoOfSlice);
//        Scalar height = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, Height);
//        Scalar bOR = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, BOR);
//        int noOfCellsHeight = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NoOfCellsHeight);


//        SpacerRadius spacerR(noOfCellsR, outerRadius, wellRadius);
//        SpacerDefault spacerAngle(noOfSlice, angle, 0.0);
//        SpacerDefault spacerZ(noOfCellsHeight, height, bOR);



//        Dune::shared_ptr<Grid> gridp(gf.createGrid());
        gridPtr() = CakeGridCreator<TypeTag, SpacerRadius>::createCakeGrid();
//        gridPtr() = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cellRes);
    }

    static std::shared_ptr<Grid> createCakeGrid()
    {
        GridFactory gf;
        Dune::FieldVector<Scalar, dim> pos;

        pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

        Dune::GeometryType type;
        type.makeHexahedron();
        std::vector<unsigned int> vid(8);
        vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3;
        vid[4] = 4; vid[5] = 5; vid[6] = 6; vid[7] = 7;
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
