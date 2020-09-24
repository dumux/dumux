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
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the Brooks Corey constitutive relations.
 */
#ifndef DUMUX_BROOKS_COREY_PARAMS_HH
#define DUMUX_BROOKS_COREY_PARAMS_HH

#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *        for the Brooks Corey constitutive relations.
 * \see BrooksCorey
 */
// "Use new material laws! Removal after 3.3")
template <class ScalarT>
class BrooksCoreyParams
{
public:
    using Scalar = ScalarT;

    BrooksCoreyParams() = default;

    BrooksCoreyParams(Scalar pe, Scalar lambda)
    : pe_(pe), lambda_(lambda)
    {}

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        return Dune::FloatCmp::eq(pe_, otherParams.pe(), /*eps*/1e-6*pe_)
               && Dune::FloatCmp::eq(lambda_, otherParams.lambda(), /*eps*/1e-6*lambda_);
    }

    /*!
     * \brief Returns the entry pressure in \f$\mathrm{[Pa]}\f$
     */
    Scalar pe() const
    { return pe_; }

    /*!
     * \brief Set the entry pressure in \f$\mathrm{[Pa]}\f$]
     */
    void setPe(Scalar v)
    { pe_ = v; }


    /*!
     * \brief Returns the lambda shape parameter \f$\mathrm{[-]}\f$
     */
    Scalar lambda() const
    { return lambda_; }

    /*!
     * \brief Set the lambda shape parameter \f$\mathrm{[-]}\f$
     */
    void setLambda(Scalar v)
    { lambda_ = v; }

private:
    Scalar pe_;
    Scalar lambda_;
};
} // namespace Dumux

#endif
