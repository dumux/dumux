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
#ifndef DUMUX_2P2C_BOX_TRAITS_HH
#define DUMUX_2P2C_BOX_TRAITS_HH

#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc used for
 *        the TwoPTwoCBoxModel.
 */
namespace Dune
{
    /*!
     * \brief Specify the shape functions, operator assemblers, etc used for
     *        the TwoPTwoCBoxModel.
     */
    template<class Scalar, class Grid>
    class TwoPTwoCBoxTraits
    {
    private:
        typedef typename Grid::ctype     CoordScalar;

        enum {
            GridDim = Grid::dimension,
            WorldDim = Grid::dimensionworld,
        };
    
    public:
        enum {
            PrimaryVariables = 2,

            NumPhases     = 2,
            NumComponents = 2
            
            // Phase indices
            WettingIndex = 0,
            NonwettingIndex = 1,

            // Bit masks for the phase state. If a bit is set, the
            // respective fluid is present in a phase.
            WettingMask    = 0x01,
            NonwettingMask = 0x02
        };
        
        //! A vector of all primary variables at a point
        typedef Dune::FieldVector<Scalar, NumUnknowns>  UnknownsVector;
        //! boundary condition vector
        typedef Dune::FieldVector<Dune::BoundaryConditions::Flags, 
                                  NumUnknowns> BoundaryTypeVector;
        
        
        //! The finite volume cell segments of a finite element cell
        typedef Dune::FVElementGeometry<Grid>  FVElementGeometry;

        /*!
         * \brief A single of shape function used for the BoxFunction
         *        inside cells.
         */
        typedef Dune::LagrangeShapeFunctionSetContainer<CoordScalar,
                                                        Scalar,
                                                        GridDim> ShapeFunctionSetContainer;
        
        /*!
         * \brief The actual shape functions which are being used. 
         * 
         * If a grid only contains simplices or tetrahedra, it is more
         * efficent to use LagrangeShapeFunctions::p1cube, or
         * LagrangeShapeFunctions::p1simplex
         *
         * TODO: Use specialization to take advantage of simplex grids
         *       and structured grids.
        */
        static const ShapeFunctionSetContainer &shapeFunctions()
            { 
                return Dune::LagrangeShapeFunctions<CoordScalar, Scalar, GridDim>::general;
            }

    private:        
        enum {
            _maxDOF = ShapeFunctionSetContainer::maxsize
        }

    public:
        
        /*!
         * \brief Represents a local solution.
         *
         * A solution vector is attached at each vertex of the cell.
         */
        struct LocalFunction
        {
            UnknownsVector atSubContVol[_maxDOF];
        };
        
        /*!
         * \brief Data which is attached to each vertex of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct CachedSubContVolData
        {
            Scalar satN;
            Scalar satW;
            Scalar pW;
            Scalar pC;
            Scalar pN;
            Scalar temperature;
            UnknownsVector mobility;  //Vector with the number of phases
            UnknownsVector density;
            FieldMatrix<Scalar, NumComponents, NumPhases> massfrac;
            int phasestate;
        };
        
        /*!
         * \brief Cached data for the each vertex of the cell.
         */
        struct CachedCellData
        {
            CachedSubContVolData  atSubContVol[_maxDOF];
        };

        /*!
         * \brief The function which represents a solution for a fixed
         *        time step. We use first-order vertex centered FE
         *        polynomials.
         */
        typedef Dune::LeafP1Function<Grid, 
                                     Scalar, 
                                     NumUnknowns>   SpatialFunction;
        
        /*!
         * \brief The OperatorAssembler which assembles the global
         *        stiffness matrix.
         */
        typedef Dune::LeafP1OperatorAssembler<Grid,
                                              Scalar,
                                              NumUnknowns>  JacobianAssembler;
    };
}

#endif
