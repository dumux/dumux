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
 * \ingroup InputOutput
 * \brief Interface for plotting the two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_IO_PLOT_PC_KR_SW_HH
#define DUMUX_IO_PLOT_PC_KR_SW_HH

#include <string>
#include <tuple>
#include <algorithm>
#include <dumux/common/math.hh>

namespace Dumux {

namespace Detail {
template<class Function, class Range>
Range evalFunctionForRange(const Function& f, const Range& range)
{
    Range result = range;
    std::transform(range.begin(), range.end(), result.begin(), [&](auto x){ return f(x); });
    return result;
}
} // end namespace Detail

/*!
 * \ingroup InputOutput
 * \brief sample the pc-sw curve
 */
template<class PcKrSw, class V>
auto samplePcSw(const PcKrSw& curve, const V& sw)
{
    return Detail::evalFunctionForRange([&](const auto s){ return curve.pc(s); }, sw);
}

/*!
 * \ingroup InputOutput
 * \brief sample the pc-sw curve derivative wrt sw
 */
template<class PcKrSw, class V>
auto samplePcSwDerivative(const PcKrSw& curve, const V& sw)
{
    return Detail::evalFunctionForRange([&](const auto s){ return curve.dpc_dsw(s); }, sw);
}

/*!
 * \ingroup InputOutput
 * \brief sample the sw-pc curve derivative wrt pc
 */
template<class PcKrSw, class V>
auto samplePcSwInverseDerivative(const PcKrSw& curve, const V& pc)
{
    return Detail::evalFunctionForRange([&](const auto p){ return curve.dsw_dpc(p); }, pc);
}

/*!
 * \ingroup InputOutput
 * \brief sample sw-pc curve but return the log10 of the capillary pressure
 */
template<class PcKrSw, class V>
auto sampleLog10PcSw(const PcKrSw& curve, const V& sw)
{
    return Detail::evalFunctionForRange([&](const auto s){ using std::log10; return log10(curve.pc(s)); }, sw);
}

/*!
 * \ingroup InputOutput
 * \brief sample krw-sw and krn-sw curves
 */
template<class PcKrSw, class V>
auto sampleRelPerms(const PcKrSw& curve, const V& sw)
{
    return std::make_tuple(
        Detail::evalFunctionForRange([&](const auto s){ return curve.krw(s); }, sw),
        Detail::evalFunctionForRange([&](const auto s){ return curve.krn(s); }, sw)
    );
}

/*!
 * \ingroup InputOutput
 * \brief sample the derivatives of the krw-sw and krn-sw curves
 */
template<class PcKrSw, class V>
auto sampleRelPermDerivatives(const PcKrSw& curve, const V& sw)
{
    return std::make_tuple(
        Detail::evalFunctionForRange([&](const auto s){ return curve.dkrw_dsw(s); }, sw),
        Detail::evalFunctionForRange([&](const auto s){ return curve.dkrn_dsw(s); }, sw)
    );
}

// forward declaration
template<class S> class GnuplotInterface;

namespace Detail {
template<class S, class V>
void addDataSetToGnuplot(GnuplotInterface<S>& gnuplot,
                         const V& x, const V& y,
                         const std::string& curveName,
                         const std::string& curveOptions,
                         const std::string& xLabel,
                         const std::string& yLabel)
{
    gnuplot.setXlabel(xLabel);
    gnuplot.setYlabel(yLabel);
    gnuplot.addDataSetToPlot(x, y, curveName, curveOptions);
}
} // end namespace Detail

/*!
 * \ingroup InputOutput
 * \brief Helper functions related to gnuplot
 */
namespace Gnuplot {

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addPcSw(GnuplotInterface<S>& gnuplot, const V& sw, const V& pc,
             const std::string& curveName = "pc-sw",
             const std::string& curveOptions = "w l",
             const std::string& xLabel = "wetting phase saturation [-]",
             const std::string& yLabel = "capillary pressure [Pa]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, pc, curveName, curveOptions, xLabel, yLabel);
}

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addPcSwDerivative(GnuplotInterface<S>& gnuplot, const V& sw, const V& dpc_dsw,
                       const std::string& curveName = "dpc-dsw",
                       const std::string& curveOptions = "w l",
                       const std::string& xLabel = "wetting phase saturation [-]",
                       const std::string& yLabel = "derivative of capillary pressure [Pa]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, dpc_dsw, curveName, curveOptions, xLabel, yLabel);
}

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addPcSwInverseDerivative(GnuplotInterface<S>& gnuplot, const V& sw, const V& dpc_dsw,
                              const std::string& curveName = "dsw-dpc",
                              const std::string& curveOptions = "w l",
                              const std::string& xLabel = "capillary pressure [Pa]",
                              const std::string& yLabel = "derivative of saturation [1/Pa]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, dpc_dsw, curveName, curveOptions, xLabel, yLabel);
}

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addLog10PcSw(GnuplotInterface<S>& gnuplot, const V& sw, const V& log10pc,
                  const std::string& curveName = "log10-pc-sw",
                  const std::string& curveOptions = "w l",
                  const std::string& xLabel = "wetting phase saturation [-]",
                  const std::string& yLabel = "log10 of capillary pressure [Pa]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, log10pc, curveName, curveOptions, xLabel, yLabel);
}

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addRelPerms(GnuplotInterface<S>& gnuplot, const V& sw, const V& krw, const V& krn,
                 const std::string& curveName = "relperm",
                 const std::string& curveOptions = "w l",
                 const std::string& xLabel = "wetting phase saturation [-]",
                 const std::string& yLabel = "relative permeability [-]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, krw, curveName + "_krw", curveOptions, xLabel, yLabel);
    Detail::addDataSetToGnuplot(gnuplot, sw, krn, curveName + "_krn", curveOptions, xLabel, yLabel);
}

/*!
 * \ingroup InputOutput
 * \brief Convenience function for adding material law quantities to gnuplot
 */
template<class S, class V>
void addRelPermDerivatives(GnuplotInterface<S>& gnuplot, const V& sw, const V& krw, const V& krn,
                           const std::string& curveName = "relperm_dsw",
                           const std::string& curveOptions = "w l",
                           const std::string& xLabel = "wetting phase saturation [-]",
                           const std::string& yLabel = "derivative of the relative permeability [-]")
{
    Detail::addDataSetToGnuplot(gnuplot, sw, krw, curveName + "_krw", curveOptions, xLabel, yLabel);
    Detail::addDataSetToGnuplot(gnuplot, sw, krn, curveName + "_krn", curveOptions, xLabel, yLabel);
}

} // end namespace Gnuplot
} // end namespace Dumux

#endif
