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
 * \brief Class to specify the type of a boundary for multidomain problems.
 */
#ifndef MULTIDOMAIN_BOUNDARY_TYPES_HH
#define MULTIDOMAIN_BOUNDARY_TYPES_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup BC
 * \brief Class to specify the type of a boundary for multidomain problems.
 */
template <int numEq>
class MultidomainBoundaryTypes : public BoundaryTypes<numEq>
{
    typedef BoundaryTypes<numEq> ParentType;

public:
    MultidomainBoundaryTypes()
    { reset(); }

    /*!
     * \brief Reset the boundary types.
     *
     * After this method no equations will be disabled and neither
     * neumann nor dirichlet conditions will be evaluated. This
     * corrosponds to a neumann zero boundary.
     */
    void reset()
    {
        ParentType::reset();
        for (int i=0; i < numEq; ++i) {
            boundaryCouplingInfo_[i].visited = 0;
            boundaryCouplingInfo_[i].isCouplingDirichlet = 0;
            boundaryCouplingInfo_[i].isCouplingNeumann = 0;
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling Dirichlet.
     */
    void setAllCouplingDirichlet()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingDirichlet(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling Neumann.
     */
    void setAllCouplingNeumann()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingNeumann(eqIdx);
        }
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling inflow.
     */
    void setCouplingDirichlet(int eqIdx)
    {
        boundaryCouplingInfo_[eqIdx].visited = 1;
        boundaryCouplingInfo_[eqIdx].isCouplingDirichlet = 1;
        boundaryCouplingInfo_[eqIdx].isCouplingNeumann = 0;

        Valgrind::SetDefined(boundaryCouplingInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling outflow.
     */
    void setCouplingNeumann(int eqIdx)
    {
        boundaryCouplingInfo_[eqIdx].visited = 1;
        boundaryCouplingInfo_[eqIdx].isCouplingDirichlet = 0;
        boundaryCouplingInfo_[eqIdx].isCouplingNeumann = 1;

        Valgrind::SetDefined(boundaryCouplingInfo_[eqIdx]);
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        Dirichlet coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isCouplingDirichlet(unsigned eqIdx) const
    { return boundaryCouplingInfo_[eqIdx].isCouplingDirichlet; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Dirichlet coupling condition.
     */
    bool hasCouplingDirichlet() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryCouplingInfo_[i].isCouplingDirichlet)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        Neumann coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isCouplingNeumann(unsigned eqIdx) const
    { return boundaryCouplingInfo_[eqIdx].isCouplingNeumann; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Neumann coupling condition.
     */
    bool hasCouplingNeumann() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryCouplingInfo_[i].isCouplingNeumann)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        Mortar coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isCouplingMortar(unsigned eqIdx) const
    {
        return false;
    }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Mortar coupling condition.
     */
    bool hasCouplingMortar() const
    {
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isCoupling(unsigned eqIdx) const
    {
        return boundaryCouplingInfo_[eqIdx].isCouplingDirichlet
               || boundaryCouplingInfo_[eqIdx].isCouplingNeumann;
    }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        coupling condition.
     */
    bool hasCoupling() const
    {
        for (int i = 0; i < numEq; ++i)
            if (isCoupling(i))
                return true;
        return false;
    }

private:
    // this is a bitfield structure!
    struct __attribute__((__packed__)) {
        unsigned char visited : 1;
        unsigned char isCouplingDirichlet : 1;
        unsigned char isCouplingNeumann : 1;
    } boundaryCouplingInfo_[numEq];
};

}

#endif
