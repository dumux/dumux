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
 *
 * \brief A preconditioner for overlapping matrices and vectors
 */
#ifndef DUMUX_OVERLAPPING_PRECONDITIONER_HH
#define DUMUX_OVERLAPPING_PRECONDITIONER_HH

#warning This file is deprecated and will be removed after Dumux 2.9

#include <dumux/common/exceptions.hh>
#include <dune/istl/preconditioners.hh>

#include "overlappingscalarproduct.hh"

namespace Dumux {

template <class SeqPreCond, class Overlap>
class OverlappingPreconditioner :
    public Dune::Preconditioner<typename SeqPreCond::domain_type,
                                typename SeqPreCond::range_type>
{
public:
    typedef typename SeqPreCond::domain_type domain_type;
    typedef typename SeqPreCond::range_type range_type;

    enum { category = Dune::SolverCategory::overlapping };

    OverlappingPreconditioner(SeqPreCond &seqPreCond, const Overlap &overlap)
        : seqPreCond_(seqPreCond), overlap_(&overlap)
    {
    }

    void pre(domain_type &x, range_type &y)
    {
        seqPreCond_.pre(x, y);

/*
        // communicate the results on the overlap
        x.syncAverage();
        y.syncAverage();
*/
    };

    void apply(domain_type &x, const range_type &d)
    {
#if HAVE_MPI
        if (overlap_->peerSet().size() > 0) {
            // set the residual and right hand side on the front to zero
            range_type dd(d);
            dd.resetFront();

            // make sure that all processes react the same if the
            // sequential preconditioner on one process throws an
            // exception
            short success;
            try {
                // execute the sequential preconditioner
                seqPreCond_.apply(x, dd);
                short localSuccess = 1;
                MPI_Allreduce(&localSuccess, // source buffer
                              &success, // destination buffer
                              1, // number of objects in buffers
                              MPI_SHORT, // data type
                              MPI_MIN, // operation
                              MPI_COMM_WORLD); // communicator
            }
            catch (const Dune::Exception &e) {
                std::cout << "Process " << overlap_->myRank()
                          << " threw exception in sequential preconditioner: " << e.what() << "\n";
                short localSuccess = 0;
                MPI_Allreduce(&localSuccess, // source buffer
                              &success, // destination buffer
                              1, // number of objects in buffers
                              MPI_SHORT, // data type
                              MPI_MIN, // operation
                              MPI_COMM_WORLD); // communicator
            }

            if (success)
                x.syncAverage();
            else
                DUNE_THROW(NumericalProblem,
                           "Preconditioner threw an exception on some process.");
        }
        else
#endif // HAVE_MPI
            seqPreCond_.apply(x, d);
    };

    void post(domain_type &x)
    {
        seqPreCond_.post(x);
    };

private:
    SeqPreCond seqPreCond_;
    const Overlap *overlap_;
};

} // namespace Dumux

#endif
