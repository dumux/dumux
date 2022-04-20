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

#include <list>
#include <memory>
#include <algorithm>

namespace Dumux {

//! Interface for observing classes
class Observer
{
    friend class Observers;
    virtual void update() = 0;
public:
    virtual ~Observer() = default;
};

//! A list of observers to be used by an observee
class Observers
{
public:
    //! attach a new observer
    void attach(Observer* o)
    {
        if (!observing_(o))
            observers_.push_back(o);
    }

    //! detach a given observer
    void detach(Observer* o)
    {
        observers_.remove(o);
    }

    //! notify all observers that the subject has been updated
    void notifyAll() const
    {
        std::for_each(
            observers_.begin(), observers_.end(),
            [&](Observer* obs) { obs->update(); }
        );
    }

    bool empty() const
    {
        return observers_.empty();
    }

private:
    bool observing_(Observer* o) const
    {
        return std::any_of(
            observers_.begin(), observers_.end(),
            [=](Observer* obs){ return obs == o; }
        );
    }

    std::list<Observer*> observers_;
};

//! Interface for classes that can be observed
class Observee
{
public:
    Observee()
    : observers_(std::make_unique<Observers>())
    {}

    ~Observee() noexcept(true)
    {
        if (!observers_->empty())
        {
            std::cerr << "The list of observers is not empty! "
                      << "Detach all observers before destroying the observee."
                      << std::endl;

            std::exit(1);
        }
    }

    Observers& observers() const
    {
        return *observers_;
    }

private:
    std::unique_ptr<Observers> observers_;
};

} // end namespace Dumux

#endif
