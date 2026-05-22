// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Constraints related to Dirichlet boundaries
 */
#ifndef DUMUX_DIRICHLET_CONSTRAINTS_HH
#define DUMUX_DIRICHLET_CONSTRAINTS_HH

#include <concepts>
#include <unordered_map>

#include <dumux/common/indextraits.hh>

namespace Dumux {

template<class Info, class Values>
struct ConstraintData
{
    using ConstraintInfo = Info;
    using ConstraintValues = Values;

    const ConstraintInfo& constraintInfo() const
    { return info_; }

    const ConstraintValues& values() const
    { return values_; }

    ConstraintInfo info_;
    ConstraintValues values_;
};

template<class DirichletConstraintInfo, class DirichletValues, class IndexType>
struct DirichletConstraintData : public ConstraintData<DirichletConstraintInfo, DirichletValues>
{
    using GridIndexType = IndexType;

    GridIndexType dofIndex() const
    { return dofIdx_; }

    GridIndexType dofIdx_;
};

} // end namespace Dumux

#endif
