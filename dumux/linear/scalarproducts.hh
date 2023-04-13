// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Classes to compute scalar products
 */
#ifndef DUMUX_LINEAR_SCALAR_PRODUCTS_HH
#define DUMUX_LINEAR_SCALAR_PRODUCTS_HH

#include <array>
#include <memory>
#include <algorithm>

#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/solvercategory.hh>
#include <dune/istl/scalarproducts.hh>

namespace Dumux {

/*!
 * \brief A scalar product for multi-type vectors
 *
 * Consistent vectors in interior and border are assumed
 * \tparam X The type of the sequential vector to use for the left hand side,
 * e.g. Dune::MultiTypeBlockVector or another type fulfilling its interface
 * \tparam C The type of the communication object
 * This must either be Dune::OwnerOverlapCopyCommunication or a type
 * implementing the same interface
 *
 * Dune::OwnerOverlapCopyCommunication can represent a overlapping or
 * a non-overlapping decomposition. This class allows to use different
 * types of decompositions for each sub-domain of the vector.
 */
template<class X, class C>
class ParallelMultiTypeScalarProduct : public Dune::ScalarProduct<X>
{
    static constexpr std::size_t numSubDomains = X::size();
public:
    // public aliases for a dune-istl like interface
    using domain_type = X;
    using field_type = typename X::field_type;
    using real_type = typename Dune::FieldTraits<field_type>::real_type;
    using communication_type = C;

    ParallelMultiTypeScalarProduct (const std::array<std::shared_ptr<const communication_type>, numSubDomains>& comms)
    : comms_(comms)
    {}

    /*!
     * \brief Dot product of two vectors
     * It is assumed that the vectors are consistent on the interior+border partition
     * According to Blatt and Bastian (2009)
     * https://doi.org/10.1504/IJCSE.2008.021112 they only have to be in a
     * "valid representation" (i.e. all dofs owned by the process have the same value as the global vector)
     */
    field_type dot (const X& x, const X& y) const override
    {
        field_type result = 0.0;

        // use the communicators for the subdomain scalar products and sum up
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i)
        {
            field_type subResult = 0.0;
            comms_[i]->dot(x[i], y[i], subResult);
            result += subResult;
        });

        return result;
    }

    /*!
     * \brief compute 2-norm of a right-hand side vector
     */
    real_type norm (const X& x) const override
    {
        using std::sqrt;
        return sqrt(dot(x, x));
    }

    /*!
     * \brief category of the scalar product
     *
     * see Dune::SolverCategory::Category
     *
     * as we have potentially several categories choose overlapping
     * if there is no clear category
     * This is part of a check mechanism and the category has to
     * match with the linear operator and preconditioner when used
     * in a parallel solver.
     */
    Dune::SolverCategory::Category category() const override
    {
        if (std::all_of(comms_.begin(), comms_.end(),
            [](auto& c){ return c->category() == Dune::SolverCategory::sequential; }
        ))
            return Dune::SolverCategory::sequential;
        else if (std::all_of(comms_.begin(), comms_.end(),
            [](auto& c){ return c->category() == Dune::SolverCategory::nonoverlapping; }
        ))
            return Dune::SolverCategory::nonoverlapping;
        else
            return Dune::SolverCategory::overlapping;
    }

private:
    std::array<std::shared_ptr<const communication_type>, numSubDomains> comms_;
};

} // end namespace Dumux

#endif
