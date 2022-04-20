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

// forward declaration
template<class Subject>
class Observers;

//! Interface for observing classes
template<class Subject>
class Observer
{
    friend class Observers<Subject>;
    virtual void update_(const Subject&) = 0;
public:
    virtual ~Observer() = default;
};

//! A list of observers to be used by an observee
template<class Subject>
class Observers
{
    using Obs = Observer<Subject>;
public:
    //! attach a new observer
    void attach(Obs* o)
    {
        if (!observing_(o))
            observers_.push_back(o);
    }

    //! detach a given observer
    void detach(Obs* o)
    {
        observers_.remove(o);
    }

    //! notify all observers that the subject has been updated
    void notifyAll(const Subject& subject) const
    {
        std::for_each(
            observers_.begin(), observers_.end(),
            [&](Obs* obs) { obs->update_(subject); }
        );
    }

    bool empty() const
    {
        return observers_.empty();
    }

private:
    bool observing_(Obs* o) const
    {
        return std::any_of(
            observers_.begin(), observers_.end(),
            [=](Obs* obs){ return obs == o; }
        );
    }

    std::list<Obs*> observers_;
};

//! Interface for classes that can be observed
template<class Impl>
class Observee
{
    using Obss = Observers<Impl>;
    using Obs = Observer<Impl>;
public:
    Observee()
    : observers_(std::make_unique<Obss>())
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

    void attach(Obs* o) const
    {
       observers_->attach(o);
    }

    //! detach a given observer
    void detach(Obs* o) const
    {
       observers_->detach(o);
    }

protected:
    //! notify all observers that the subject has been updated
    void notifyAllObservers_(const Impl& impl) const
    {
       observers_->notifyAll(impl);
    }

private:
    std::unique_ptr<Obss> observers_;
};

} // end namespace Dumux

#endif
