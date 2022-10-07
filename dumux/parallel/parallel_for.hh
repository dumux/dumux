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
 * \ingroup Parallel
 * \brief Parallel for loop (multithreading)
 */

#ifndef DUMUX_PARALLEL_PARALLEL_FOR_HH
#define DUMUX_PARALLEL_PARALLEL_FOR_HH

#include <dumux/parallel/multithreading.hh>

#if HAVE_CPP_PARALLEL_ALGORITHMS
#include <algorithm>
#include <execution>
#include <dune/common/rangeutilities.hh>
#endif

#if HAVE_TBB
#include <tbb/parallel_for.h>
#endif

#if HAVE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

// contents of the detail namespace might change
// any time without prior notice (do not use directly)
#ifndef DOXYGEN // hide from doxygen
namespace Dumux::Detail {

// This should be specialized for different ExecutionBackends
template<class FunctorType, class ExecutionBackend>
class ParallelFor;


// Serial backend implementation
template<class FunctorType>
class ParallelFor<FunctorType, Multithreading::ExecutionBackends::Serial>
{
public:
    ParallelFor(const std::size_t count, const FunctorType& functor)
    : functor_(functor), count_(count) {}

    void execute() const
    {
        for (std::size_t i = 0; i < count_; ++i)
            functor_(i);
    }

private:
    FunctorType functor_;
    std::size_t count_;
};

#if HAVE_CPP_PARALLEL_ALGORITHMS
// C++ parallel algorithms backend implementation
template<class FunctorType>
class ParallelFor<FunctorType, Multithreading::ExecutionBackends::Cpp>
{
public:
    ParallelFor(const std::size_t count, const FunctorType& functor)
    : functor_(functor), range_(count) {}

    void execute() const
    {
        std::for_each(std::execution::par_unseq, range_.begin(), range_.end(), functor_);
    }

private:
    FunctorType functor_;
    Dune::IntegralRange<std::size_t> range_;
};
#endif


#if HAVE_TBB
// TBB backend implementation
template<class FunctorType>
class ParallelFor<FunctorType, Multithreading::ExecutionBackends::TBB>
{
public:
    ParallelFor(const std::size_t count, const FunctorType& functor)
    : functor_(functor), count_(count) {}

    void execute() const
    {
        tbb::parallel_for(std::size_t{0}, count_, [&](const std::size_t i){ functor_(i); });
    }

private:
    FunctorType functor_;
    std::size_t count_;
};
#endif // HAVE_TBB

#if HAVE_KOKKOS
// Kokkos backend implementation
template<class FunctorType>
class ParallelFor<FunctorType, Multithreading::ExecutionBackends::Kokkos>
{
public:
    ParallelFor(const std::size_t count, const FunctorType& functor)
    : functor_(functor), count_(count) {}

    void execute() const
    {
        Kokkos::parallel_for(count_, [&](const std::size_t i){ functor_(i); });
    }

private:
    FunctorType functor_;
    std::size_t count_;
};
#endif // HAVE_KOKKOS


#if HAVE_OPENMP
// OpenMP backend implementation
template<class FunctorType>
class ParallelFor<FunctorType, Multithreading::ExecutionBackends::OpenMP>
{
public:
    ParallelFor(const std::size_t count, const FunctorType& functor)
    : functor_(functor), count_(count) {}

    void execute() const
    {
        #pragma omp parallel for
        for (std::size_t i = 0; i < count_; ++i)
            functor_(i);
    }

private:
    FunctorType functor_;
    std::size_t count_;
};
#endif // HAVE_OPENMP


} // end namespace Detail
#endif // DOXYGEN


namespace Dumux {

/*!
 * \ingroup Parallel
 * \brief A parallel for loop (multithreading)
 * \param count the number of work tasks to perform
 * \param functor functor executed for each task (get task number as argument)
 */
template<class FunctorType>
inline void parallelFor(const std::size_t count, const FunctorType& functor)
{
    using ExecutionBackend = Detail::Multithreading::ExecutionBackend;
    Detail::ParallelFor<FunctorType, ExecutionBackend> action(count, functor);
    action.execute();
}

} // end namespace Dumux

#endif
