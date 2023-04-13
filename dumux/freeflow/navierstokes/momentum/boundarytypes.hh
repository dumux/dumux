// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesBoundaryTypes
 */
#ifndef FREEFLOW_NAVIERSTOKES_MOMENTUM_BOUNDARY_TYPES_HH
#define FREEFLOW_NAVIERSTOKES_MOMENTUM_BOUNDARY_TYPES_HH

#include <dumux/common/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Class to specify the type of a boundary condition for the Navier-Stokes model.
 */
template <int size>
class NavierStokesMomentumBoundaryTypes : public BoundaryTypes<size>
{
    using ParentType = BoundaryTypes<size>;

public:
    NavierStokesMomentumBoundaryTypes()
    {
        for (int eqIdx=0; eqIdx < size; ++eqIdx)
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
        for (int eqIdx=0; eqIdx < size; ++eqIdx)
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
        for (int i = 0; i < size; ++i)
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

    std::array<NavierStokesBoundaryInfo, size> boundaryInfo_;
};

} // end namespace Dumux

#endif
