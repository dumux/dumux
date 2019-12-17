// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup IMPET
 * \ingroup IMPETProperties mpfa
 * \file
 *
 * \brief Properties for a MPFA method.
 */
#ifndef DUMUX_FVMPFAPROPERTIES_HH
#define DUMUX_FVMPFAPROPERTIES_HH

// dumux environment
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
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

#if HAVE_DUNE_ALUGRID
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
namespace TTag {
struct MPFAProperties {};
}

template<class TypeTag, class MyTypeTag>
struct  GridTypeIndices  { using type = UndefinedProperty; };//!< The grid type indices to decide which grid is used
template<class TypeTag, class MyTypeTag>
struct  GridImplementation  { using type = UndefinedProperty; }; //!< Gives kind of grid implementation in form of a GridType
template<class TypeTag, class MyTypeTag>
struct  MPFAInteractionVolume  { using type = UndefinedProperty; }; //!< Type of the data container for one interaction volume
template<class TypeTag, class MyTypeTag>
struct  MPFAInteractionVolumeContainer  { using type = UndefinedProperty; }; //!< Type of the data container for one interaction volume
template<class TypeTag, class MyTypeTag>
struct  MPFATransmissibilityCalculator  { using type = UndefinedProperty; }; //!< Type defining the transmissibility calculation
}
}

namespace Dumux
{
namespace Properties
{

//! \cond \private
template<class TypeTag>
struct GridImplementation<TypeTag, TTag::MPFAProperties>
{
private:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
public:
    static const int value = GridImp<Grid, Grid::dimension>::imp;
};
//! \endcond

//! Set grid type indices
template<class TypeTag>
struct GridTypeIndices<TypeTag, TTag::MPFAProperties> { using type = GridTypes; };
}
}// end of Dune namespace
#endif
