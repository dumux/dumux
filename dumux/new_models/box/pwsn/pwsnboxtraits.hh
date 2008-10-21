/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_PWSN_BOX_TRAITS_HH
#define DUMUX_PWSN_BOX_TRAITS_HH

#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc used for
 *        the PwSnBoxModel.
 */
namespace Dune
{
    /*!
     * \brief Specify the shape functions, operator assemblers, etc used for
     *        the PwSnBoxModel.
     */
    template<class Scalar, class Grid>
    class PwSnBoxTraits
    {
    private:
        typedef typename Grid::ctype     CoordScalar;

        enum {
            GridDim = Grid::dimension,
            WorldDim = Grid::dimensionworld,

            SIZE = Dune::LagrangeShapeFunctionSetContainer<CoordScalar,
                                                           Scalar,
                                                           Grid::dimension>::maxsize
        };
    
    public:
        enum {
            NumUnknowns = 2,
            PwIndex = 0,
            SnIndex = 1
        };
        
        //! A vector of all primary variables at a point
        typedef Dune::FieldVector<Scalar, NumUnknowns>  UnknownsVector;
        //! boundary condition vector
        typedef Dune::FieldVector<Dune::BoundaryConditions::Flags, 
                                  NumUnknowns> BoundaryTypeVector;
        
        
        //! The finite volume dual cells of a finite element cell
        typedef Dune::FVElementGeometry<Grid>  FVElementGeometry;
        
        /*!
         * \brief Represents a local solution.
         *
         * A solution vector is attached at each vertex of the cell.
         */
        struct LocalFunction
        {
            UnknownsVector atSubContVol[SIZE];
        };
        
        /*!
         * \brief Data which is attached to each vertex of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct CachedSubContVolData
        {
            Scalar Sw;
            Scalar pC;
            Scalar pN;
            
            UnknownsVector mobility;  //Vector with the number of phases
        };
        
        /*!
         * \brief Cached data for the each vertex of the cell.
         */
        struct CachedCellData
        {
            CachedSubContVolData  atSubContVol[SIZE];
        };

        //! The function which represents a solution for a fixed time
        //! step. We use first-order vertex centered FE polynomials.
        typedef Dune::LeafP1Function<Grid, 
                                     Scalar, 
                                     NumUnknowns>   SpatialFunction;
        
        //! The OperatorAssembler which assembles the global stiffness
        //! matrix
        typedef Dune::LeafP1OperatorAssembler<Grid,
                                              Scalar,
                                              NumUnknowns>  JacobianAssembler;

        //! a single of shape function used for the BoxFunction inside
        //! cells.
        typedef Dune::LagrangeShapeFunctionSetContainer<CoordScalar,
                                                        Scalar,
                                                        GridDim> ShapeFunctionSetContainer;
        
        //! The actual shape functions which are being used. If a 
        //! grid only contains simplices or tetrahedra, it is more
        //! efficent to use LagrangeShapeFunctions::p1cube, or
        //! LagrangeShapeFunctions::p1simplex
        //!
        //! TODO: Use specialization to take advantage of simplex grids
        //!       and structured grids.
        static const ShapeFunctionSetContainer &shapeFunctions()
            { 
                return Dune::LagrangeShapeFunctions<CoordScalar, Scalar, GridDim>::general;
            }
    };
}

#endif
