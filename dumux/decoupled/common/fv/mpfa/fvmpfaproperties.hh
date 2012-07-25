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
 * \ingroup IMPET
 * \ingroup IMPETProperties
 * \file
 *
 * \brief Properties for a MPFA method.
 */
#ifndef DUMUX_FVMPFAPROPERTIES_HH
#define DUMUX_FVMPFAPROPERTIES_HH

// dumux environment
#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid/2d/alugrid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

namespace Dumux
{
/*!
 *
 *
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
//! \cond \private
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

#if HAVE_ALUGRID
template<int dim>
struct GridImp<Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>, dim>
{
    static const int imp = GridTypes::aluGrid;
};
#endif

#if HAVE_UG
template<int dim>
struct GridImp<Dune::UGGrid<dim>, dim>
{
    static const int imp = GridTypes::ugGrid;
};
#endif
//! \endcond

namespace Properties
{
//! Basic Type tag for MFPA models
NEW_TYPE_TAG(MPFAProperties);

NEW_PROP_TAG( GridTypeIndices );//!< The grid type indices to decide which grid is used
NEW_PROP_TAG( GridImplementation ); //!< Gives kind of grid implementation in form of a GridType
NEW_PROP_TAG( EnableComplexLStencil ); //!< Enable use of four different L-shapes instead of 2 (3-d)
NEW_PROP_TAG( MPFATransmissibilityCriterionThreshold ); //!< Threshold for transmissibility choice
NEW_PROP_TAG( MPFATransmissibilityCriterion ); //!< Choose transmissibility criterion
}
}

namespace Dumux
{
namespace Properties
{

//! \cond \private
SET_PROP(MPFAProperties, GridImplementation)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
public:
    static const int value = GridImp<Grid, Grid::dimension>::imp;
};
//! \endcond

//! Set grid type indices
SET_TYPE_PROP(MPFAProperties, GridTypeIndices, GridTypes);

//! Allow use of all available L-shapes
SET_BOOL_PROP(MPFAProperties, EnableComplexLStencil, true);

//! Allow use of all available L-shapes
SET_SCALAR_PROP(MPFAProperties, MPFATransmissibilityCriterionThreshold, 1e-8);

//! Set standard criterion
SET_INT_PROP(MPFAProperties, MPFATransmissibilityCriterion, 0);
}
}// end of Dune namespace
#endif
