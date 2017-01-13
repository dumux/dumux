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
        (*this)[cellCenterIdx] = std::forward<decltype(ccPriVars)>(ccPriVars);
        (*this)[faceIdx] = std::forward<decltype(facePriVars)>(facePriVars);
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

/*!
 * \brief This class generates a range [a,b) which can be used in a for loop, e.g.
 *        for(auto i : range(3) { ... i=0,1,2 }   or
 *        for(auto i : range(5, 8) { ... i=5,6,7 }
 *        see: https://en.wikipedia.org/wiki/Generator_(computer_programming)
 *        TODO: should this be moved to the math.hh header?
 */
class range
{
public:
    // constructors
    range(int end): last_(end), iter_(0) {}
    range(int begin, int end): last_(end), iter_(begin) {}
    range() = delete;

    // Iterable functions
     const range& begin() const { return *this; }
     const range& end() const { return *this; }

    // Iterator functions
    bool operator!=(const range&) const { return iter_ < last_; }
    void operator++() { ++iter_; }
    int operator*() const { return iter_; }
private:
    int last_;
    int iter_;
};

/*!
* \brief Class which provides two ranges of indices (cc and face)
*        cc: for(auto i : PriVarIndices(cellCenterIdx)) { ... }
*        face: for(auto i : PriVarIndices(faceIdx)) { ... }
*/
template<class TypeTag>
class PriVarIndices : public range
{
    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    using cellCenterIdx = typename DofTypeIndices::CellCenterIdx;
    using faceIdx = typename DofTypeIndices::FaceIdx;

    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
    static constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);

public:
    PriVarIndices(cellCenterIdx) : range(0, numEqCellCenter) {}
    PriVarIndices(faceIdx) : range(numEqCellCenter, numEqCellCenter + numEqFace) {}
};

} // namespace Dumux

#endif
