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
 * \brief Primary Variables for the Staggered Grid models
 */
#ifndef DUMUX_STAGGERED_PRIMARYVARIABLES_HH
#define DUMUX_STAGGERED_PRIMARYVARIABLES_HH

#include "properties.hh"
#include <dune/istl/multitypeblockvector.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesModel
 * \brief This class inherits from DUNE's MultiTypeBlockVector and provides a specific [] operator for convenience
 */
template<class TypeTag,
         class CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables),
         class FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables)>
class StaggeredPrimaryVariables : public Dune::MultiTypeBlockVector<CellCenterPrimaryVariables, FacePrimaryVariables>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices =  typename GET_PROP_TYPE(TypeTag, Indices);
    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using PrimaryVariables = Dune::MultiTypeBlockVector<CellCenterPrimaryVariables, FacePrimaryVariables>;

public:
    StaggeredPrimaryVariables() = default;

    StaggeredPrimaryVariables(const Scalar value)
    {
        (*this)[cellCenterIdx] = value;
        (*this)[faceIdx] = value;
    }

    const auto& operator [](const unsigned int i) const
    {
        if(i < Indices::faceOffset)
            return PrimaryVariables::operator[](cellCenterIdx)[i];
        else
            return PrimaryVariables::operator[](faceIdx)[i - Indices::faceOffset];
    }

    auto& operator [](const unsigned int i)
    {
        if(i < Indices::faceOffset)
            return PrimaryVariables::operator[](cellCenterIdx)[i];
        else
            return PrimaryVariables::operator[](faceIdx)[i - Indices::faceOffset];
    }
};

} // namespace Dumux

#endif
