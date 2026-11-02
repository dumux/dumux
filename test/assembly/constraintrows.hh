// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Detection of constraint rows in assembled Jacobians
 */
#ifndef DUMUX_TEST_ASSEMBLY_CONSTRAINT_ROWS_HH
#define DUMUX_TEST_ASSEMBLY_CONSTRAINT_ROWS_HH

#include <cstddef>

namespace Dumux::Test {

/*!
 * \brief Whether the matrix row i is the row of a constraint, i.e. the identity in its diagonal block
 *        and, for the constraint of a periodic degree of freedom, the negative identity in the block of
 *        the degree of freedom it is mapped to
 * \param periodic Only report constraints of periodic degrees of freedom
 */
template<class MatrixRow>
bool isConstraintRow(const MatrixRow& row, std::size_t i, bool periodic = false)
{
    bool hasNegativeIdentity = false;
    for (auto col = row.begin(); col != row.end(); ++col)
    {
        bool isZero = true, isIdentity = true, isNegativeIdentity = true;
        for (int r = 0; r < int(col->N()); ++r)
            for (int c = 0; c < int(col->M()); ++c)
            {
                const auto value = (*col)[r][c];
                const double diag = r == c ? 1.0 : 0.0;
                isZero = isZero && value == 0.0;
                isIdentity = isIdentity && value == diag;
                isNegativeIdentity = isNegativeIdentity && value == -diag;
            }

        if (col.index() == i && !isIdentity)
            return false;
        if (col.index() != i && !isZero)
        {
            if (!isNegativeIdentity || hasNegativeIdentity)
                return false;
            hasNegativeIdentity = true;
        }
    }
    return !periodic || hasNegativeIdentity;
}

} // end namespace Dumux::Test

#endif
