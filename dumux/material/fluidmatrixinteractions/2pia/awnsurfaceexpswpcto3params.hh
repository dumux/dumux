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
 * This is of third order in pc.
 */
#ifndef AWN_SURFACE_EXP_SW_PC_TO_3_PARAMS
#define AWN_SURFACE_EXP_SW_PC_TO_3_PARAMS

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of interfacial area surface params
 */
template<class ScalarT>
class AwnSurfaceExpSwPcTo3Params
{
public:
    using Scalar = ScalarT;

    AwnSurfaceExpSwPcTo3Params()
    {}

    AwnSurfaceExpSwPcTo3Params(const Scalar a1,
                               const Scalar a2,
                               const Scalar a3)
    {
        setA1(a1);
        setA2(a2);
        setA3(a3);
    }

     /*!
     * \brief Return the \f$\mathrm{[a_{1}}\f$ shape parameter of awn surface.
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

private:
    Scalar Swr_;
    Scalar Snr_;
    Scalar a1_;
    Scalar a2_;
    Scalar a3_;
};
} // namespace Dumux

#endif
