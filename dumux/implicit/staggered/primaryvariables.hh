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

    using ParentType = Dune::MultiTypeBlockVector<CellCenterPrimaryVariables, FacePrimaryVariables>;

public:
    StaggeredPrimaryVariables() = default;

     /*!
     * \brief Constructor to initialize all entries with the same value
     *
     * \param value The value
     */
    StaggeredPrimaryVariables(const Scalar value) noexcept
    {
        (*this)[cellCenterIdx] = value;
        (*this)[faceIdx] = value;
    }

     /*!
     * \brief Constructor to initialize the cellcenter and face primary values with given values
     *
     * \param ccPriVars The cellcenter primary variables used for initialization
     * \param facePriVars The face primary variables used for initialization
     */
    StaggeredPrimaryVariables(CellCenterPrimaryVariables&& ccPriVars, FacePrimaryVariables&& facePriVars) noexcept
    {
        (*this)[cellCenterIdx] = std::move(ccPriVars);
        (*this)[faceIdx] = std::move(facePriVars);
    }

     /*!
     * \brief Operator overload which allows to automatically access the "right" priVars vector via pvIdx.
     *        const version
     * \note: the ParentType (DUNE multitypeblockvector) [] operator has to be visible (using ...)
     *
     * \param pvIdx The global index of the primary variable
     */
    using ParentType::operator [];
    const Scalar& operator [](const unsigned int pvIdx) const
    {
        if(pvIdx < Indices::faceOffset)
            return ParentType::operator[](cellCenterIdx)[pvIdx];
        else
            return ParentType::operator[](faceIdx)[pvIdx - Indices::faceOffset];
    }

     /*!
     * \brief Operator overload which allows to automatically access the "right" priVars vector via pvIdx
     *        non-const version
     *
     * \param pvIdx The global index of the primary variable
     */
    Scalar& operator [](const unsigned int pvIdx)
    {
        if(pvIdx < Indices::faceOffset)
            return ParentType::operator[](cellCenterIdx)[pvIdx];
        else
            return ParentType::operator[](faceIdx)[pvIdx - Indices::faceOffset];
    }
};

} // namespace Dumux

#endif
