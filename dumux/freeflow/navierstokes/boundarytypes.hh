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
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesBoundaryTypes
 */
#ifndef FREEFLOW_NAVIERSTOKES_BOUNDARY_TYPES_HH
#define FREEFLOW_NAVIERSTOKES_BOUNDARY_TYPES_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Class to specify the type of a boundary condition for the Navier-Stokes model.
 */
template <int numEq>
class NavierStokesBoundaryTypes : public BoundaryTypes<numEq>
{
    using ParentType = BoundaryTypes<numEq>;

public:
    NavierStokesBoundaryTypes()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
            resetEq(eqIdx);
    }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(const int eqIdx)
    {
        ParentType::resetEq(eqIdx);

        boundaryInfo_[eqIdx].isSymmetry = false;
        boundaryInfo_[eqIdx].isBeaversJoseph = false;
    }

    /*!
     * \brief Sets a symmetry boundary condition for all equations
     */
    void setAllSymmetry()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
        {
            resetEq(eqIdx);
            boundaryInfo_[eqIdx].isSymmetry = true;
        }
    }

    /*!
     * \brief Returns true if the there is a symmetry boundary condition
     */
    bool isSymmetry() const
    { return boundaryInfo_[0].isSymmetry; }

    /*!
     * \brief  Prevent setting all boundary conditions to Dirichlet.
     */
    template<class T = void>
    void setAllDirichlet()
    {
        static_assert(AlwaysFalse<T>::value, "Setting all boundary types to Dirichlet not permitted!");
    }

    /*!
     * \brief  Prevent setting all boundary conditions to Neumann.
     */
    template<class T = void>
    void setAllNeumann()
    {
        static_assert(AlwaysFalse<T>::value, "Setting all boundary types to Neumann not permitted!");
    }

    /*!
     * \brief Set a boundary condition for a single equation to
     *        Beavers-Joseph(-Saffmann) (special case of Dirichlet b.c.).
     */
    void setBeaversJoseph(const int eqIdx)
    {
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].isBeaversJoseph = true;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        Beavers-Joseph(-Saffman) boundary condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isBeaversJoseph(const int eqIdx) const
    { return boundaryInfo_[eqIdx].isBeaversJoseph; }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        Beavers-Joseph(-Saffman) boundary condition.
     */
    bool hasBeaversJoseph() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isBeaversJoseph)
                return true;
        return false;
    }

protected:
    //! use bitfields to minimize the size
    struct NavierStokesBoundaryInfo
    {
        bool isSymmetry : 1;
        bool isBeaversJoseph : 1;
    };

    std::array<NavierStokesBoundaryInfo, numEq> boundaryInfo_;
};

} // end namespace Dumux

#endif
