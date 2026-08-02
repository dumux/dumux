// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief A compressed subdomain-local view of a global residual
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_RESIDUAL_VIEW_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_RESIDUAL_VIEW_HH

#include <cstddef>
#include <vector>

#include <dune/istl/bvector.hh>

#include <dumux/nonlinear/preconditioning/partition.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief A compressed, subdomain-owned residual keyed by global degree-of-freedom index
 *
 * Entries addressed outside \f$ N_i \f$ are discarded, mirroring SubdomainMatrixView's row handling:
 * an assembly loop extended to the fringe evaluates residuals that are not wanted.
 */
template<class PV>
class SubdomainResidualView
{
public:
    using PrimaryVariables = PV;
    using LocalVector = Dune::BlockVector<PrimaryVariables>;
    using Index = Partition::Index;

    SubdomainResidualView(const Partition& partition, SubdomainIndex i)
    : rowOf_(partition.numDofs(), Partition::noIndex)
    {
        for (Index k = 0; k < partition.dofs(i).size(); ++k)
            rowOf_[partition.dofs(i)[k]] = k;

        vector_.resize(partition.dofs(i).size());
        vector_ = 0.0;
    }

    PrimaryVariables& operator[] (Index globalI)
    {
        const auto local = rowOf_[globalI];
        return local == Partition::noIndex ? sink_ : vector_[local];
    }

    const PrimaryVariables& operator[] (Index globalI) const
    {
        const auto local = rowOf_[globalI];
        return local == Partition::noIndex ? sink_ : vector_[local];
    }

    std::size_t size() const { return vector_.size(); }
    void setToZero() { vector_ = 0.0; }

    const LocalVector& vector() const { return vector_; }
    LocalVector& vector() { return vector_; }

private:
    LocalVector vector_;
    std::vector<Index> rowOf_;
    mutable PrimaryVariables sink_;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
