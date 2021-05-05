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
 * \ingroup TimeStepping
 * \brief Class that represents a time level during time integration.
 */
#ifndef DUMUX_TIMESTEPPING_TIME_LEVEL_HH
#define DUMUX_TIMESTEPPING_TIME_LEVEL_HH

namespace Dumux::Experimental {

/*!
 * \brief Class that represents a time level during time integration.
 */
template<class Scalar>
class TimeLevel
{
public:

    /*!
     * \brief Construct a time level with a time.
     * \note This can be used in contexts outside of time integration,
     *       where no information on a previous time or time step size is needed.
     */
    explicit TimeLevel(Scalar curTime)
    : curTime_(curTime)
    , prevTime_(curTime)
    , timeStepFraction_(1.0)
    {}

    /*!
     * \brief Construct a time level with information on an ongoing time step.
     * \param curTime The current time level
     * \param prevTime The previous time level
     * \param dtFraction The fraction of a time step this level corresponds to.
     * \note Within a time integration step, several time levels might occur
     *       when multi-stage methods are used. The argument dtFraction allows
     *       for determining the time that will be reached at the end of the
     *       time integration step.
     */
    TimeLevel(Scalar curTime, Scalar prevTime, Scalar dtFraction)
    : curTime_(curTime)
    , prevTime_(prevTime)
    , timeStepFraction_(dtFraction)
    {}

    //! Return the current time
    Scalar current() const { return curTime_; }
    //! Return the time at the beginning of time integration
    Scalar previous() const { return prevTime_; }
    //! Return the fraction of the time step this level corresponds to
    Scalar timeStepFraction() const { return timeStepFraction_; }

private:
    Scalar curTime_;
    Scalar prevTime_;
    Scalar timeStepFraction_;
};

} // end namespace Dumux::Experimental

#endif
