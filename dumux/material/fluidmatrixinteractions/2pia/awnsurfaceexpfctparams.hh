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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the parameters for a function relating volume specific interfacial area to capillary pressure and saturation.
 * This function is exponential.
 */
#ifndef AWN_SURFACE_EXP_FCT_PARAMS_HH
#define AWN_SURFACE_EXP_FCT_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief implementation of interfacial area surface params
 */
template<class ScalarT>
class AwnSurfaceExpFctParams
{
public:
    using Scalar = ScalarT;

    AwnSurfaceExpFctParams()
    {}

    AwnSurfaceExpFctParams(const Scalar a1, const Scalar a2, const Scalar a3)
    {
        setA1(a1);
        setA2(a2);
        setA3(a3);
    }

    /*!
         * \brief Return the \f$\mathrm{a_{1}}\f$ shape parameter of awn surface.
         */
        const Scalar a1() const
        { return a1_; }
        /*!
         * \brief Return the \f$\mathrm{a_{2}}\f$ shape parameter of awn surface.
         */
        const Scalar a2() const
        { return a2_; }
        /*!
         * \brief Return the \f$\mathrm{a_{3}}\f$ shape parameter of awn surface.
         */
        const Scalar a3() const
        { return a3_; }
        /*!
         * \brief Return the \f$\mathrm{S_{wr}}\f$ shape parameter of awn surface.
         */
        const Scalar Swr() const
        { return Swr_; }

    /*!
     * \brief Set the \f$\mathrm{a_{1}}\f$ shape parameter.
     */
    void setA1(const Scalar v)
    { a1_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{2}}\f$ shape parameter.
     */
    void setA2(const Scalar v)
    { a2_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{3}}\f$ shape parameter.
     */
    void setA3(const Scalar v)
    { a3_ = v; }

    /*!
     * \brief Set the \f$\mathrm{S_{wr}}\f$ for the surface.
     */
    void setSwr(const Scalar v)
    { Swr_ = v; }

private:
    Scalar a1_;
    Scalar a2_;
    Scalar a3_;
    Scalar Swr_;
};
} // namespace Dumux

#endif
