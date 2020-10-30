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
 * This parametrization is a second order polynomial.
 */
#ifndef AWN_SURFACE_POLYNOMIAL_2ND_ORDER_PARAMS_HH
#define AWN_SURFACE_POLYNOMIAL_2ND_ORDER_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of interfacial area surface params
 */
template<class ScalarT>
class AwnSurfacePolynomial2ndOrderParams
{
public:
    using Scalar = ScalarT;

    AwnSurfacePolynomial2ndOrderParams()
    {}

    AwnSurfacePolynomial2ndOrderParams(Scalar a00, Scalar a10, Scalar a20, Scalar a11, Scalar a01, Scalar a02)
    {
        setA00(a00);
        setA10(a10);
        setA20(a20);
        setA11(a11);
        setA01(a01);
        setA02(a02);
    }

    /*!
     * \brief Return the \f$\mathrm{a_{00}}\f$ shape parameter of awn surface.
     */
    const Scalar a00() const
    { return a00_; }

    /*!
     * \brief Return the \f$\mathrm{a_{10}}\f$ shape parameter of awn surface.
     */
    const Scalar a10() const
    { return a10_; }

    /*!
     * \brief Return the \f$\mathrm{a_{20}}\f$ shape parameter of awn surface.
     */
    const Scalar a20() const
    { return a20_; }

    /*!
     * \brief Return the \f$\mathrm{a_{11}}\f$ shape parameter of awn surface.
     */
    const Scalar a11() const
    { return a11_; }

    /*!
     * \brief Return the \f$\mathrm{a_{01}}\f$ shape parameter of awn surface.
     */
    const Scalar a01() const
    { return a01_; }

    /*!
     * \brief Return the \f$\mathrm{a_{02}}\f$ shape parameter of awn surface.
     */
    const Scalar a02() const
    { return a02_; }

    /*!
     * \brief Set the \f$\mathrm{a_{00}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA00(const Scalar v)
    { a00_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{10}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA10(const Scalar v)
    { a10_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{20}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA20(const Scalar v)
    { a20_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{11}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA11(const Scalar v)
    { a11_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{01}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA01(const Scalar v)
    { a01_ = v; }

    /*!
     * \brief Set the \f$\mathrm{a_{02}}\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setA02(const Scalar v)
    { a02_ = v; }

private:
    Scalar a00_;
    Scalar a10_;
    Scalar a20_;
    Scalar a11_;
    Scalar a01_;
    Scalar a02_;
};
} // namespace Dumux

#endif
