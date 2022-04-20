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
 * \brief Base classes for realizing the pull-variant of the observer pattern.
 */
#ifndef DUMUX_COMMON_OBSERVER_HH
#define DUMUX_COMMON_OBSERVER_HH

#include <vector>
#include <memory>
#include <algorithm>

namespace Dumux {

//! Interface for observing classes
struct Observer
{
    virtual ~Observer() = default;
    virtual void update() = 0;
};

using ObserverPointer = std::shared_ptr<Observer>;

//! Interface for classes that can be observed
class Observable
{
public:
    virtual ~Observable() = default;

    //! Register a new observer
    void addObserver(const ObserverPointer& observer)
    {
        if (!observesAlready_(observer))
            observers_.push_back(observer);
    }

protected:
    //! Allows parent classes to notify all observers after a change
    void notifyAllObservers_() const
    {
        std::for_each(observers_.begin(),
                      observers_.end(),
                      [] (auto& obs) { obs->update(); });
    }

private:
    bool observesAlready_(const ObserverPointer& observer) const
    {
        return std::any_of(observers_.begin(),
                           observers_.end(),
                           [&] (const auto& obs) { return obs == observer; });
    }

    std::vector<ObserverPointer> observers_;
};

} // end namespace Dumux
#endif
