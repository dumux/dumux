// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup Parallel
 * \brief Parallel for loop (multithreading)
 */

#ifndef DUMUX_PARALLEL_PARALLEL_FOR_HH
#define DUMUX_PARALLEL_PARALLEL_FOR_HH

#include <dumux/parallel/multithreading.hh>

#if DUMUX_HAVE_CPP_PARALLEL_ALGORITHMS
#include <algorithm>
#include <execution>
#include <dune/common/rangeutilities.hh>
#endif

#if HAVE_TBB
#include <tbb/parallel_for.h>
#endif

#if DUMUX_HAVE_KOKKOS
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

#if DUMUX_HAVE_CPP_PARALLEL_ALGORITHMS
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

#if DUMUX_HAVE_KOKKOS
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
#endif // DUMUX_HAVE_KOKKOS


#if DUMUX_HAVE_OPENMP
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
#endif // DUMUX_HAVE_OPENMP


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
