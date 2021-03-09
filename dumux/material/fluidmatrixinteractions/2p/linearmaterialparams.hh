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
 * \brief   Parameters for the linear capillary pressure and
 *          relative permeability <-> saturation relations
 */
#ifndef LINEAR_MATERIAL_PARAMS_HH
#define LINEAR_MATERIAL_PARAMS_HH

// TODO Deprecated. Remove after 3.3

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Reference implementation of params for the linear material
 *        law.
 */
template<class ScalarT>
class LinearMaterialParams
{
public:
    using Scalar = ScalarT;

    LinearMaterialParams()
    {}

    LinearMaterialParams(Scalar entryPc, Scalar maxPc)
    {
        setEntryPc(entryPc);
        setMaxPc(maxPc);
    }


    /*!
     * \brief Return the entry pressure for the linear material law in \f$\mathrm{[Pa]}\f$.
     *
     * The entry pressure is reached at \f$\mathrm{\overline{S}_w = 1}\f$
     */
    Scalar entryPc() const
    { return entryPc_; }

    /*!
     * \brief Set the entry pressure for the linear material law in \f$\mathrm{[Pa]}\f$.
     *
     * The entry pressure is reached at \f$\mathrm{\overline{S}_w = 1}\f$
     */
    void setEntryPc(Scalar v)
    { entryPc_ = v; }

    /*!
     * \brief Return the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    Scalar maxPc() const
    { return maxPc_; }

    /*!
     * \brief Set the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    void setMaxPc(Scalar v)
    { maxPc_ = v; }

private:
    Scalar entryPc_;
    Scalar maxPc_;
};
} // namespace Dumux

#endif
