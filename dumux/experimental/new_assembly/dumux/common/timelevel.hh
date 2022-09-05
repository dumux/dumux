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
 * \ingroup Common
 * \copydoc Dumux::TimeLevel
 */
#ifndef DUMUX_COMMON_TIME_LEVEL_HH
#define DUMUX_COMMON_TIME_LEVEL_HH

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Class that represents a time level, possibly within a time integration step.
 */
template<typename S>
class TimeLevel
{
public:
    using Scalar = S;

    /*!
     * \brief Construct a time level with information on an ongoing time step.
     * \param curTime The current time level
     * \param prevTime The previous time level
     * \param dtFraction The fraction of the ongoing time step this level corresponds to.
     */
    TimeLevel(Scalar curTime, Scalar prevTime, Scalar dtFraction)
    : curTime_(curTime)
    , prevTime_(prevTime)
    , timeStepFraction_(dtFraction)
    {}

    /*!
     * \brief Construct a time level.
     * \param curTime The current time level
     */
    explicit TimeLevel(Scalar curTime)
    : TimeLevel(curTime, curTime, 1.0)
    {}

    //! Return the current time
    Scalar current() const { return curTime_; }
    //! Return the time at the beginning of the time step
    Scalar previous() const { return prevTime_; }
    //! Return the fraction of the ongoing time step this level corresponds to
    Scalar timeStepFraction() const { return timeStepFraction_; }

private:
    Scalar curTime_;
    Scalar prevTime_;
    Scalar timeStepFraction_;
};

} // namespace Dumux

#endif
