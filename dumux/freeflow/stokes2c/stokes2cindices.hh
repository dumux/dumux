// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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

/*!
 * \file
 *
 * \brief Defines the indices required for the compositional Stokes BOX model.
 */
#ifndef DUMUX_STOKES2C_INDICES_HH
#define DUMUX_STOKES2C_INDICES_HH

#include <dumux/freeflow/stokes/stokesindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxIndices
 * \brief The common indices for the compositional stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct Stokes2cCommonIndices : public StokesCommonIndices<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Phase indices
    static const int lPhaseIdx = FluidSystem::lPhaseIdx; //!< Index of the liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int lCompIdx = 0; //!< Index of the liquid's primary component
    static const int gCompIdx = 1; //!< Index of the gas' primary component

    // equation and primary variable indices
    static const int dim = StokesCommonIndices<TypeTag>::dim;
    static const int transportIdx = PVOffset + dim+1; //! The index for the transport equation.
    static const int massOrMoleFracIndex = transportIdx; //! The index for the mass or mole fraction in primary variable vectors.

    static const int phaseIdx = gPhaseIdx; //!< Index of the gas phase (required to use the same fluid system in coupled models)
};
} // end namespace

#endif
