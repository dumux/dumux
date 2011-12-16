// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
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

/*!
 * \ingroup MPFA2p
 * \ingroup Properties
 * \file
 *
 * \brief Properties for the MPFA-O method.
 */
#ifndef DUMUX_MPFAOPROPERTIES_HH
#define DUMUX_MPFAOPROPERTIES_HH

// dumux environment
#include <dumux/common/basicproperties.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/alugrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

namespace Dumux
{
/*!
 * \brief Indices denoting the different grid types.
 */
struct GridTypes
{
public:
    //
    static const int any = 0;
    //SGrid
    static const int sGrid = 1;
    //YaspGrid
    static const int yaspGrid = 2;
    //UGGrid
    static const int ugGrid = 3;
    //ALUGrid
    static const int aluGrid = 4;
};

template<class Grid, int dim>
struct GridImp
{
    static const int imp = GridTypes::any;
};

template<int dim>
struct GridImp<Dune::YaspGrid<dim>, dim>
{
    static const int imp = GridTypes::yaspGrid;
};

template<int dim>
struct GridImp<Dune::SGrid<dim, dim>, dim>
{
    static const int imp = GridTypes::sGrid;
};

#if HAVE_UG
template<int dim>
struct GridImp<Dune::UGGrid<dim>, dim>
{
    static const int imp = GridTypes::ugGrid;
};
#endif

namespace Properties
{
NEW_TYPE_TAG(MPFAProperties);

NEW_PROP_TAG( GridTypeIndices );
NEW_PROP_TAG( GridImplementation ); //returns kind of grid implementation

//forward declaration
NEW_PROP_TAG( Grid );

SET_PROP(MPFAProperties, GridImplementation)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
public:
    static const int value = GridImp<Grid, Grid::dimension>::imp;
};

SET_TYPE_PROP(MPFAProperties, GridTypeIndices, GridTypes);

}
}// end of Dune namespace
#endif
