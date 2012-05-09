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
 * \brief Defines the indices required for the compositional Stokes box model.
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
 * \brief The common indices for the compositional Stokes box model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct Stokes2cCommonIndices : public StokesCommonIndices<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Phase index
    static const int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx); //!< Index of the employed phase in case of a two-phase fluidsystem (set by default to nPhase)

    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< \deprecated use phaseIdx instead, only one phase employed
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< \deprecated use phaseIdx instead, only one phase employed
    static const int lPhaseIdx = wPhaseIdx; //!< \deprecated use phaseIdx instead, only one phase employed
    static const int gPhaseIdx = nPhaseIdx; //!< \deprecated use phaseIdx instead, only one phase employed

    // Component indices
    static const int phaseCompIdx = phaseIdx; //!< The index of the main component of the considered phase
    static const int transportCompIdx = (unsigned int)(1-phaseIdx); //!< The index of the transported (minor) component; ASSUMES phase indices of 0 and 1

    static const int comp1Idx = 0; //!< \deprecated Index of the wetting's primary component
    static const int comp0Idx = 1; //!< \deprecated Index of the non-wetting's primary component
    static const int lCompIdx = transportCompIdx; //!< \deprecated use transportComp instead
    static const int gCompIdx = phaseCompIdx; //!< \deprecated use phaseCompIdx instead

    // equation and primary variable indices
    static const int dim = StokesCommonIndices<TypeTag>::dim;
    static const int transportEqIdx = PVOffset + dim+1; //!< The index for the transport equation
    static const int transportIdx = transportEqIdx; //!< \deprecated use transportEqIdx instead

    static const int massOrMoleFracIdx = transportEqIdx; //!< The index of the mass or mole fraction of the transported component in primary variable vectors
    static const int massOrMoleFracIndex = massOrMoleFracIdx; //!< \deprecated use massOrMoleFracIdx instead

};
} // end namespace

#endif
