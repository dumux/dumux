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
/*
 * author: Larissa de Vries
 * last update: 05.12.2016
 */

#ifndef DUMUX_CAKE_GRID_CREATOR_HH
#define DUMUX_CAKE_GRID_CREATOR_HH

#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/common/propertysystem.hh>
#include <vector>
#include <iostream>
#include <cmath>

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
    typedef std::vector<Scalar> vector;
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
        std::cout << "makeGrid() starts..." << std::endl;
        vector dR;
        vector dA;
        vector dZ;
        createVectors(dR, dA, dZ);
        gridPtr() = CakeGridCreator<TypeTag>::createCakeGrid(dR, dA, dZ);
        std::cout << "makeGrid() ends..." << std::endl;
    }

    static void createVectors(vector &dR, vector &dA, vector &dZ)
    {
        GlobalPosition polarCoorMin;
        GlobalPosition polarCoorMax;
        CellArray cells;

       try{
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


        if (polarCoorMin[0] > polarCoorMax[0]
            || polarCoorMin[1] > polarCoorMax[1]
            || polarCoorMin[2] > polarCoorMax[2]){
            std::cout << "Error: Input vectors are not ascending." << std::endl;
        }

        std::cout << "createVectors() starts..." << std::endl;

        // save array entries as doubles
        double dr = polarCoorMax[0] - polarCoorMin[0];
        double da = polarCoorMax[1] - polarCoorMin[1];
        double dz = polarCoorMax[2] - polarCoorMin[2];
        double cellR = cells[0];
        double cellA = cells[1];
        double cellZ = cells[2];

        // calculate step size in r, a and z direction
        double stepsizeR = (dr/cellR);
        double stepsizeA = (da/cellA);
        double stepsizeZ = (dz/cellZ);

        dR.resize(cellR+1, 0.0);
        dA.resize(cellA+1, 0.0);
        dZ.resize(cellZ+1, 0.0);

        dR[0] = polarCoorMin[0];
        for (int i = 1; i <= cellR; ++i){
            dR[i] = dR[0]+i*stepsizeR;
        }

        dA[0] = polarCoorMin[1];
        for (int j = 1; j <= cellA; ++j){
            dA[j] = (dA[0]+j*stepsizeA)/180*M_PI;
        }

        dZ[0] = polarCoorMin[2];
        for (int k = 1; k <= cellZ; ++k){
            dZ[k] = dZ[0]+k*stepsizeZ;
        }

        std::cout << "createVectors() ends..." << std::endl;

    }

    static std::shared_ptr<Grid> createCakeGrid(vector &dR, vector &dA, vector &dZ)
    {
        /* std::cout << dR[0] << std::endl;
        std::cout << dA[0] << std::endl;
        std::cout << dZ[0] << std::endl;

        std::cout << dR[1] << std::endl;
        std::cout << dA[1] << std::endl;
        std::cout << dZ[1] << std::endl; */

        std::cout << "createCakeGrid() starts..." << std::endl;
        GridFactory gf;
        Dune::GeometryType type;
        type.makeHexahedron();

        std::cout << "create nodes starts..." << std::endl;

        //create nodes
        for (int j = 0; j <= dA.size() - 1; ++j){
            for (int l = 0; l <= dZ.size() - 1; ++l){
                for (int i = 0; i <= dR.size()- 1; ++i){
                    // Get wellRadius
                    Scalar wellRadius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, WellRadius);
                    // transformation
                    double dX = cos(dA[j])*wellRadius + cos(dA[j])*dR[i];
                    double dY = sin(dA[j])*wellRadius + sin(dA[j])*dR[i];
                    Dune::FieldVector <double, 3> v = {dX, dY, dZ[l]};
                    std::cout << "Coordinates of : " << dX << " " << dY << " " << dZ[l] << std::endl;
                    gf.insertVertex(v);

                }
            }

        }



        std::cout << "create nodes ends..." << std::endl;



        // assign nodes
        int z = 0;
        int t = 0;
        std::cout << "assign nodes starts..." << std::endl;
        for (int j = 0; j < dA.size() - 1; ++j){
            for (int l = 0; l < dZ.size() - 1; ++l){
                if (j < dA.size() - 2){
                    for (int i = 0; i < dR.size() - 1; ++i){
                        std::cout << "assign nodes not 360 starts..." << std::endl;
                        std::vector<unsigned int> vid({z, z+1, z+dR.size()*dZ.size(), z+dR.size()*dZ.size()+1, z+dR.size(), z+dR.size()+1, z+dR.size()*dZ.size()+dR.size(), z+dR.size()*dZ.size()+dR.size()+1});

                        for (int i = 0; i < vid.size(); ++i)
                        {
                            std::cout << "vid = " << vid[i] << std::endl;
                        }

                        gf.insertElement(type, vid);
                        std::cout << "assign nodes not 360 ends..." << std::endl;

                        z = z+1;
                    }
                    z = z+1;
                    std::cout << "assign nodes ends..." << std::endl;

                }
                // TODO else wird immer im letzten Schritt betreten...
                else {
                    // assign nodes for 360°-cake
                    std::cout << "assign nodes 360 starts..." << std::endl;
                    for (int i = 0; i < dR.size() - 1; ++i){
                        // z = z + 1;
                        std::vector<unsigned int> vid({z, z+1, t, t+1, z+dR.size(), z+dR.size()+1, t+dR.size(), t+dR.size()+1});
                        for (int i = 0; i < vid.size(); ++i)
                        {
                            std::cout << "vid = " << vid[i] << std::endl;
                        }
                        t = t + 1;
                        gf.insertElement(type, vid);
                        }
                    t = t + 1;
                    }
                    std::cout << "assign nodes 360 ends..." << std::endl;

            }

            z = z + dR.size();
        }
        // assign nodes ends...
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
