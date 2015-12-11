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
 * \brief Class to specify the type of a boundary.
 */
#ifndef BOUNDARY_TYPES_HH
#define BOUNDARY_TYPES_HH

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup BC
 * \brief Class to specify the type of a boundary.
 */
template <int numEq>
class BoundaryTypes
{
public:
    BoundaryTypes()
    { reset(); }

    /*!
     * \brief Reset the boundary types.
     *
     * After this method no equations will be disabled and neither
     * Neumann nor Dirichlet conditions will be evaluated. This
     * corrosponds to a Neumann zero boundary.
     */
    void reset()
    {
        for (int i=0; i < numEq; ++i) {
            boundaryInfo_[i].visited = 0;

            boundaryInfo_[i].isDirichlet = 0;
            boundaryInfo_[i].isNeumann = 0;
            boundaryInfo_[i].isOutflow = 0;
            boundaryInfo_[i].isCouplingInflow = 0;
            boundaryInfo_[i].isCouplingOutflow = 0;
            boundaryInfo_[i].isMortarCoupling = 0;

            eq2pvIdx_[i] = i;
            pv2eqIdx_[i] = i;
        }
    }

    /*!
     * \brief Returns true if the boundary types for a given equation
     *        has been specified.
     *
     * \param eqIdx The index of the equation
     */
    bool isSet(int eqIdx) const
    { return boundaryInfo_[eqIdx].visited; }

    /*!
     * \brief Make sure the boundary conditions are well-posed.
     *
     * If they are not, an assertation fails and the program aborts!
     * (if the NDEBUG macro is not defined)
     */
    void checkWellPosed() const
    {
#ifndef NDEBUG
        for (int i=0; i < numEq; ++i)
            // if this fails, at least one condition is missing.
            assert(boundaryInfo_[i].visited);
#endif
    }

    /*!
     * \brief Set all boundary conditions to Neumann.
     */
    void setAllNeumann()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setNeumann(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to Dirichlet.
     */
    void setAllDirichlet()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx)
        {
            setDirichlet(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to Neumann.
     */
    void setAllOutflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setOutflow(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling inflow.
     */
    DUNE_DEPRECATED_MSG("setAllCouplingInflow() is deprecated")
    void setAllCouplingInflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            boundaryInfo_[eqIdx].visited = 1;
            boundaryInfo_[eqIdx].isDirichlet = 0;
            boundaryInfo_[eqIdx].isNeumann = 0;
            boundaryInfo_[eqIdx].isOutflow = 0;
            boundaryInfo_[eqIdx].isCouplingInflow = 1;
            boundaryInfo_[eqIdx].isCouplingOutflow = 0;
            boundaryInfo_[eqIdx].isMortarCoupling = 0;

            Valgrind::SetDefined(boundaryInfo_[eqIdx]);
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling outflow.
     */
    DUNE_DEPRECATED_MSG("setAllCouplingOutflow() is deprecated")
    void setAllCouplingOutflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            boundaryInfo_[eqIdx].visited = 1;
            boundaryInfo_[eqIdx].isDirichlet = 0;
            boundaryInfo_[eqIdx].isNeumann = 0;
            boundaryInfo_[eqIdx].isOutflow = 0;
            boundaryInfo_[eqIdx].isCouplingInflow = 0;
            boundaryInfo_[eqIdx].isCouplingOutflow = 1;
            boundaryInfo_[eqIdx].isMortarCoupling = 0;

            Valgrind::SetDefined(boundaryInfo_[eqIdx]);
        }
    }

    /*!
     * \brief Set all boundary conditions to mortar coupling.
     */
    DUNE_DEPRECATED_MSG("setAllMortarCoupling() is deprecated")
    void setAllMortarCoupling()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            boundaryInfo_[eqIdx].visited = 1;
            boundaryInfo_[eqIdx].isDirichlet = 0;
            boundaryInfo_[eqIdx].isNeumann = 0;
            boundaryInfo_[eqIdx].isOutflow = 0;
            boundaryInfo_[eqIdx].isCouplingInflow = 0;
            boundaryInfo_[eqIdx].isCouplingOutflow = 0;
            boundaryInfo_[eqIdx].isMortarCoupling = 1;

            Valgrind::SetDefined(boundaryInfo_[eqIdx]);
        }
    }

    /*!
     * \brief Set a Neumann boundary condition for a single a single
     *        equation.
     *
     * \param eqIdx The index of the equation
     */
    void setNeumann(int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 0;
        boundaryInfo_[eqIdx].isNeumann = 1;
        boundaryInfo_[eqIdx].isOutflow = 0;
        boundaryInfo_[eqIdx].isCouplingInflow = 0;
        boundaryInfo_[eqIdx].isCouplingOutflow = 0;
        boundaryInfo_[eqIdx].isMortarCoupling = 0;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a Dirichlet boundary condition for a single primary
     *        variable
     *
     * \param pvIdx The index of the primary variable for which the
     *              Dirichlet condition should apply.
     * \param eqIdx The index of the equation which should used to set
     *              the Dirichlet condition
     */
    void setDirichlet(int pvIdx, int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 1;
        boundaryInfo_[eqIdx].isNeumann = 0;
        boundaryInfo_[eqIdx].isOutflow = 0;
        boundaryInfo_[eqIdx].isCouplingInflow = 0;
        boundaryInfo_[eqIdx].isCouplingOutflow = 0;
        boundaryInfo_[eqIdx].isMortarCoupling = 0;

        // update the equation <-> primary variable mapping
        eq2pvIdx_[eqIdx] = pvIdx;
        pv2eqIdx_[pvIdx] = eqIdx;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a Neumann boundary condition for a single a single
     *        equation.
     *
     * \param eqIdx The index of the equation on which the outflow
     *              condition applies.
     */
    void setOutflow(int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 0;
        boundaryInfo_[eqIdx].isNeumann = 0;
        boundaryInfo_[eqIdx].isOutflow = 1;
        boundaryInfo_[eqIdx].isCouplingInflow = 0;
        boundaryInfo_[eqIdx].isCouplingOutflow = 0;
        boundaryInfo_[eqIdx].isMortarCoupling = 0;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling inflow.
     */
    DUNE_DEPRECATED_MSG("setCouplingInflow() is deprecated")
    void setCouplingInflow(int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 0;
        boundaryInfo_[eqIdx].isNeumann = 0;
        boundaryInfo_[eqIdx].isOutflow = 0;
        boundaryInfo_[eqIdx].isCouplingInflow = 1;
        boundaryInfo_[eqIdx].isCouplingOutflow = 0;
        boundaryInfo_[eqIdx].isMortarCoupling = 0;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling outflow.
     */
    DUNE_DEPRECATED_MSG("setCouplingOutflow() is deprecated")
    void setCouplingOutflow(int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 0;
        boundaryInfo_[eqIdx].isNeumann = 0;
        boundaryInfo_[eqIdx].isOutflow = 0;
        boundaryInfo_[eqIdx].isCouplingInflow = 0;
        boundaryInfo_[eqIdx].isCouplingOutflow = 1;
        boundaryInfo_[eqIdx].isMortarCoupling = 0;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to mortar coupling.
     */
    DUNE_DEPRECATED_MSG("setMortarCoupling() is deprecated")
    void setMortarCoupling(int eqIdx)
    {
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 0;
        boundaryInfo_[eqIdx].isNeumann = 0;
        boundaryInfo_[eqIdx].isOutflow = 0;
        boundaryInfo_[eqIdx].isCouplingInflow = 0;
        boundaryInfo_[eqIdx].isCouplingOutflow = 0;
        boundaryInfo_[eqIdx].isMortarCoupling = 1;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a Dirichlet boundary condition for a single primary
     *        variable.
     *
     * Depending on the discretization, setting the Dirichlet condition
     * will replace the balance equation with index equal to pvIdx.
     *
     * \param pvIdx The index of the primary variable inside a
     *              PrimaryVariables object.
     */
    void setDirichlet(int pvIdx)
    {
        setDirichlet(pvIdx, pvIdx);
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        Dirichlet condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isDirichlet(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isDirichlet; }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        Dirichlet condition.
     */
    bool hasDirichlet() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isDirichlet)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        Neumann condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isNeumann(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isNeumann; }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        Neumann condition.
     */
    bool hasNeumann() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isNeumann)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        outflow condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isOutflow(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isOutflow; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        outflow condition.
     */
    bool hasOutflow() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isOutflow)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        inflow coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    DUNE_DEPRECATED_MSG("isCouplingInflow() is deprecated")
    bool isCouplingInflow(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isCouplingInflow; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        inflow coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasCouplingInflow() is deprecated")
    bool hasCouplingInflow() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingInflow)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        outflow coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    DUNE_DEPRECATED_MSG("isCouplingOutflow() is deprecated")
    bool isCouplingOutflow(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isCouplingOutflow; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        outflow coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasCouplingOutflow() is deprecated")
    bool hasCouplingOutflow() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingOutflow)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        mortar coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    DUNE_DEPRECATED_MSG("isMortarCoupling() is deprecated")
    bool isMortarCoupling(unsigned eqIdx) const
    {
        return boundaryInfo_[eqIdx].isMortarCoupling;
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        mortar coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasMortarCoupling() is deprecated")
    bool hasMortarCoupling() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isMortarCoupling)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        coupling condition.
     *
     * \param eqIdx The index of the equation
     *
     * This only used for correcting the pressure in the Stokes equation.
     */
    bool isCoupling(unsigned eqIdx) const
    {
        return false;
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        coupling condition.
     */
    bool hasCoupling() const
    {
        return false;
    }

    /*!
     * \brief Returns the index of the equation which should be used
     *        for the Dirichlet condition of the pvIdx's primary
     *        variable.
     *
     * \param pvIdx The index of the primary variable which is be set
     *              by the Dirichlet condition.
     */
    unsigned dirichletToEqIndex(unsigned pvIdx) const
    { return pv2eqIdx_[pvIdx]; }

    /*!
     * \brief Returns the index of the primary variable which should
     *        be used for the Dirichlet condition given an equation
     *        index.
     *
     * \param eqIdx The index of the equation which is used to set
     *              the Dirichlet condition.
     */
    unsigned eqToDirichletIndex(unsigned eqIdx) const
    { return eq2pvIdx_[eqIdx]; }

private:
    // this is a bitfield structure!
    struct __attribute__((__packed__)) {
        unsigned char visited : 1;
        unsigned char isDirichlet : 1;
        unsigned char isNeumann : 1;
        unsigned char isOutflow : 1;
        unsigned char isCouplingInflow : 1;
        unsigned char isCouplingOutflow : 1;
        unsigned char isMortarCoupling : 1;
    } boundaryInfo_[numEq];

    unsigned char eq2pvIdx_[numEq];
    unsigned char pv2eqIdx_[numEq];
};

}

#endif
