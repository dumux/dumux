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
     * \brief Reset the boundary types for all equations.
     *
     * After this method no equations will be disabled and neither
     * Neumann nor Dirichlet conditions will be evaluated. This
     * corresponds to a Neumann zero boundary.
     */
    void reset()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
        {
            resetEq(eqIdx);
        }
    }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(int eqIdx)
    {
          boundaryInfo_[eqIdx].visited = 0;

          boundaryInfo_[eqIdx].isDirichlet = 0;
          boundaryInfo_[eqIdx].isNeumann = 0;
          boundaryInfo_[eqIdx].isOutflow = 0;
          boundaryInfo_[eqIdx].isCouplingDirichlet = 0;
          boundaryInfo_[eqIdx].isCouplingNeumann = 0;
          boundaryInfo_[eqIdx].isCouplingMortar = 0;

          eq2pvIdx_[eqIdx] = eqIdx;
          pv2eqIdx_[eqIdx] = eqIdx;
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
     * If they are not, an assertion fails and the program aborts!
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
     * \brief Set all boundary conditions to Dirichlet-like coupling
     */
    void setAllCouplingDirichlet()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingDirichlet(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to Neumann-like coupling.
     */
    void setAllCouplingNeumann()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingNeumann(eqIdx);
        }
    }
    /*!
     * \brief Set all boundary conditions to mortar coupling.
     */
    void setAllCouplingMortar()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingDirichlet(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling inflow.
     */
    DUNE_DEPRECATED_MSG("setAllCouplingInflow() is deprecated. Use setAllCouplingNeumann() instead.")
    void setAllCouplingInflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingNeumann(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to coupling outflow.
     */
    DUNE_DEPRECATED_MSG("setAllCouplingOutflow() is deprecated. Use setAllCouplingDirichlet() instead.")
    void setAllCouplingOutflow()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingDirichlet(eqIdx);
        }
    }

    /*!
     * \brief Set all boundary conditions to mortar coupling.
     */
    DUNE_DEPRECATED_MSG("setAllMortarCoupling() is deprecated. Use setAllCouplingMortar() instead.")
    void setAllMortarCoupling()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            setCouplingMortar(eqIdx);
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
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isNeumann = 1;

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
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isDirichlet = 1;

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
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isOutflow = 1;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to
     *        a Dirichlet-like coupling condition.
     */
    void setCouplingDirichlet(int eqIdx)
    {
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isCouplingDirichlet = 1;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to
     *        a Neumann-like coupling condition.
     */
    void setCouplingNeumann(int eqIdx)
    {
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isCouplingNeumann = 1;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to
     *        a mortar coupling condition.
     */
    void setCouplingMortar(int eqIdx)
    {
        resetEq(eqIdx);
        boundaryInfo_[eqIdx].visited = 1;
        boundaryInfo_[eqIdx].isCouplingMortar = 1;

        Valgrind::SetDefined(boundaryInfo_[eqIdx]);
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling inflow.
     */
    DUNE_DEPRECATED_MSG("setCouplingInflow() is deprecated. Use setCouplingNeumann() instead.")
    void setCouplingInflow(int eqIdx)
    {
        setCouplingNeumann(eqIdx);
    }

    /*!
     * \brief Set a boundary condition for a single equation to coupling outflow.
     */
    DUNE_DEPRECATED_MSG("setCouplingOutflow() is deprecated. Use setCouplingDirichlet() instead.")
    void setCouplingOutflow(int eqIdx)
    {
        setCouplingDirichlet(eqIdx);
    }

    /*!
     * \brief Set a boundary condition for a single equation to mortar coupling.
     */
    DUNE_DEPRECATED_MSG("setMortarCoupling() is deprecated. Use setCouplingMortar() instead.")
    void setMortarCoupling(int eqIdx)
    {
        setCouplingMortar(eqIdx);
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
    DUNE_DEPRECATED_MSG("isCouplingInflow() is deprecated. Use isCouplingNeumann() instead.")
    bool  flow(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isCouplingNeumann; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        inflow coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasCouplingInflow() is deprecated. Use hasCouplingNeumann() instead.")
    bool hasCouplingInflow() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingNeumann)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        outflow coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    DUNE_DEPRECATED_MSG("isCouplingOutflow() is deprecated. Use isCouplingDirichlet() instead.")
    bool isCouplingOutflow(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isCouplingDirichlet; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        outflow coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasCouplingOutflow() is deprecated. Use hasCouplingDirichlet() instead.")
    bool hasCouplingOutflow() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingDirichlet)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        mortar coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    DUNE_DEPRECATED_MSG("isMortarCoupling() is deprecated. Use isCouplingMortar() instead.")
    bool isMortarCoupling(unsigned eqIdx) const
    {
        return boundaryInfo_[eqIdx].isCouplingMortar;
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        mortar coupling condition.
     */
    DUNE_DEPRECATED_MSG("hasMortarCoupling() is deprecated. Use hasCouplingMortar() instead.")
    bool hasMortarCoupling() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingMortar)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if an equation is used to specify an
     *        Dirichlet coupling condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isCouplingDirichlet(unsigned eqIdx) const
    { return boundaryInfo_[eqIdx].isCouplingDirichlet; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Dirichlet coupling condition.
     */
    bool hasCouplingDirichlet() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingDirichlet)
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
    { return boundaryInfo_[eqIdx].isCouplingNeumann; }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Neumann coupling condition.
     */
    bool hasCouplingNeumann() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingNeumann)
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
        return boundaryInfo_[eqIdx].isCouplingMortar;
    }

    /*!
     * \brief Returns true if some equation is used to specify an
     *        Mortar coupling condition.
     */
    bool hasCouplingMortar() const
    {
        for (int i = 0; i < numEq; ++i)
            if (boundaryInfo_[i].isCouplingMortar)
                return true;
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
        return boundaryInfo_[eqIdx].isCouplingDirichlet
               || boundaryInfo_[eqIdx].isCouplingNeumann
               || boundaryInfo_[eqIdx].isCouplingMortar;
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        coupling condition.
     */
    bool hasCoupling() const
    {
        for (int i = 0; i < numEq; ++i)
            if (isCoupling(i))
                return true;
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

protected:
    // this is a bitfield structure!
    struct __attribute__((__packed__)) {
        unsigned char visited : 1;
        unsigned char isDirichlet : 1;
        unsigned char isNeumann : 1;
        unsigned char isOutflow : 1;
        unsigned char isCouplingDirichlet : 1;
        unsigned char isCouplingNeumann : 1;
        unsigned char isCouplingMortar : 1;
    } boundaryInfo_[numEq];

    unsigned char eq2pvIdx_[numEq];
    unsigned char pv2eqIdx_[numEq];
};

}

#endif
