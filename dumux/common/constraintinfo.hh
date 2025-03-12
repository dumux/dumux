// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Class to specify information related to constraints
 */
#ifndef DUMUX_CONSTRAINT_INFO_HH
#define DUMUX_CONSTRAINT_INFO_HH

#include <algorithm>
#include <array>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief Class to specify information related to constraints
 */
template <int numEq>
class ConstraintInfo
{
public:
    ConstraintInfo()
    { reset(); }

    //! we might have a constraint for each equation
    static constexpr int size()
    { return numEq; }

    /*!
     * \brief Reset for all equations.
     */
    void reset()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
            resetEq(eqIdx);
    }

    /*!
     * \brief Reset for one equation.
     */
    void resetEq(int eqIdx)
    {
        isConstraint_[eqIdx] = false;
        eq2pvIdx_[eqIdx] = eqIdx;
        pv2eqIdx_[eqIdx] = eqIdx;
    }

    /*!
     * \brief Set all as constraints.
     */
    void setAll()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx)
            set(eqIdx);
    }

    /*!
     * \brief Set a constraint condition for a single primary
     *        variable
     *
     * \param pvIdx The index of the primary variable for which the
     *              constraint condition should apply.
     * \param eqIdx The index of the equation which should be used to set
     *              the constraint condition
     */
    void set(int pvIdx, int eqIdx)
    {
        resetEq(eqIdx);
        isConstraint_[eqIdx] = true;

        // update the equation <-> primary variable mapping
        eq2pvIdx_[eqIdx] = pvIdx;
        pv2eqIdx_[pvIdx] = eqIdx;
    }

    /*!
     * \brief Set a constraint condition for a single primary
     *        variable.
     *
     * \param pvIdx The index of the primary variable inside a
     *              PrimaryVariables object.
     */
    void set(int pvIdx)
    { set(pvIdx, pvIdx); }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        constraint condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isConstraint(unsigned eqIdx) const
    { return isConstraint_[eqIdx]; }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        constraint condition.
     */
    bool hasConstraint() const
    {
        return std::any_of(isConstraint_.begin(),
                           isConstraint_.end(),
                           [](const auto val){ return val; }
                           );
    }

    /*!
     * \brief Returns the index of the equation which should be used
     *        for the constraint condition of the pvIdx's primary
     *        variable.
     *
     * \param pvIdx The index of the primary variable which is be set
     *              by the constraint condition.
     */
    unsigned priVarToEqIndex(unsigned pvIdx) const
    { return pv2eqIdx_[pvIdx]; }

    /*!
     * \brief Returns the index of the primary variable which should
     *        be used for the constraint condition given an equation
     *        index.
     *
     * \param eqIdx The index of the equation which is used to set
     *              the constraint condition.
     */
    unsigned eqToPriVarIndex(unsigned eqIdx) const
    { return eq2pvIdx_[eqIdx]; }

protected:
    std::array<bool, numEq> isConstraint_;
    std::array<unsigned int, numEq> eq2pvIdx_, pv2eqIdx_;
};

} // end namespace Dumux

#endif
