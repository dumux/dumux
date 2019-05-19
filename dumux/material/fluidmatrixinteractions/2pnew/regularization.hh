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
 * \brief Implementation of the regularized version of the van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_DEFAULT_REGULARIZATION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_DEFAULT_REGULARIZATION_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/optional.hh>

namespace Dumux {
namespace FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The regularization base class that conatins all the interfaces
 * \note Depending on the application not all interfaces have to be implemented
 *       by a custom regularization. The user will be notified with a runtime warning.
 */
template<class Scalar>
class TwoPRegularization
{
public:
    using Return = OptionalFloatingPointNumber<Scalar>;

    //! Init does nothing per default
    template<typename... Args> void init(Args&&...) {}

    //! No parameters per default
    template<class S> struct Params {};

    /*!
     * \brief The regularized capillary pressure-saturation curve
     */
    Return pc(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized partial derivative of the capillary pressure w.r.t. the saturation
     */
    Return dpc_dsw(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    Return sw(const Scalar pc) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized partial derivative of the saturation to the capillary pressure
     */
    Return dsw_dpc(const Scalar pc) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized relative permeability for the wetting phase
     */
    Return krw(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    Return dkrw_dsw(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized relative permeability for the non-wetting phase
     */
    Return krn(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }

    /*!
     * \brief The regularized derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    Return dkrn_dsw(const Scalar sw) const
    {
        DUNE_THROW(Dune::NotImplemented, "The regularization does not implement this function!");
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A tag to turn off regularization and it's overhead
 */
template<class Scalar>
struct NoTwoPRegularization : public TwoPRegularization<Scalar>
{};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The default regularization policy
 * \todo TODO: Implement the default old-style regularization
 */
template<class Scalar>
class TwoPDefaultRegularization : public TwoPRegularization<Scalar>
{};

} // end namespace FluidMatrix
} // end namespace Dumux

#endif
