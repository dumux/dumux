// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Implementation of the regularized version of the van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef DUMUX_WEBB_REGULARISATION_HH
#define DUMUX_WEBB_REGULARISATION_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

/*!
 * \brief Parameters that are necessary for the \em regularization of
 *        the VanGenuchten material law after Webb (2000) for low saturations.
 *        For high saturation we use a spline with linear continuation.
 * \ingroup fluidmatrixinteractionsparams
 */
template<class Scalar>
class WebbRegularizationParams
{
public:
    WebbRegularizationParams()
    {
        setOvenDryPc(1.0e9);
        setHighSweThreshold_pc(0.99);
    }

    ////////////////////////////////////////////////////////
    //! Setter methods to be called by the user
    ////////////////////////////////////////////////////////

    //! Set the capillary pressure of the oven dry material
    //! If you change this you have to call initRegularizationParams again
    void setOvenDryPc(const Scalar pc)
    { ovenDryPc_ = pc; }

    //! We set high saturation threshold (effective saturation)
    //! If you change this you have to call initRegularizationParams again
    void setHighSweThreshold_pc(const Scalar swe)
    { highSweThreshold_ = swe; }

    //! Just for compatibility, this is chosen by setting oven dry pc to fit the curve
    void setLowSweThreshold_pc(const Scalar swe)
    { }

    //! This is a free functions, being a friend we can access private members
    template<class Params>
    friend void initRegularizationParams(Params& params)
    {
        using MaterialLaw = typename Params::MaterialLaw;
        using EffToAbsPolicy = typename Params::EffToAbsPolicy;

        static constexpr Scalar min = 1e-40;
        static constexpr Scalar eps = 1e-16;

        ////////////////////////////////////////////////////////////
        //! Precompute regularization for small wetting saturations
        ////////////////////////////////////////////////////////////

        // compute the threshold saturation and capillary pressure
        // start quite close to swr to walk downwards with the Newton
        Scalar sw_star = params.swr() + 0.001*(1 - params.snr() - params.swr());
        Scalar swe_star = EffToAbsPolicy::swToSwe(params, sw_star);
        Scalar pc_star = MaterialLaw::pcRaw(params, swe_star);
        Scalar dpcdsw_star = MaterialLaw::dpc_dsweRaw(params, swe_star)*EffToAbsPolicy::dswe_dsw(params) / (std::log(10)*pc_star);

        const auto res = [](Scalar sw, Scalar pc, Scalar slope, Scalar ovenDryPc)
                         { return std::log10(ovenDryPc) - std::log10(pc) + sw*slope; };

        Scalar residual = res(sw_star, pc_star, dpcdsw_star, params.ovenDryPc());

        int iterations = 0;
        while (std::abs(residual) > 1e-5)
        {
            auto defl_sw = sw_star + eps;
            auto defl_swe = EffToAbsPolicy::swToSwe(params, defl_sw);
            auto defl_pc = MaterialLaw::pcRaw(params, defl_swe);
            auto defl_slope = MaterialLaw::dpc_dsweRaw(params, defl_swe)*EffToAbsPolicy::dswe_dsw(params) / (std::log(10)*defl_pc);
            const auto deflected = res(defl_sw, defl_pc, defl_slope, params.ovenDryPc());
            const auto derivative = (deflected-residual)/eps;

            // make sure we don't get out of physical bounds
            sw_star = std::min(1.0-min, std::max(min, sw_star - residual/derivative));
            swe_star = std::min(1.0-min, std::max(min, EffToAbsPolicy::swToSwe(params, sw_star)));
            sw_star = EffToAbsPolicy::sweToSw(params, swe_star);
            pc_star = MaterialLaw::pcRaw(params, swe_star);
            dpcdsw_star = MaterialLaw::dpc_dsweRaw(params, swe_star)*EffToAbsPolicy::dswe_dsw(params) / (std::log(10)*pc_star);

            residual = res(sw_star, pc_star, dpcdsw_star, params.ovenDryPc());

            ++iterations;

            if (iterations > 200)
                DUNE_THROW(Dune::SystemError, "Didn't succeed in regularizing the pc-sw curve using Webb's method.");
        }

        std::cout << "Regularizing pc-sw curve at sw = " << sw_star
                  << " (took " << iterations << " Newton iterations, residual = " << std::abs(residual) << ")." << std::endl;

        params.setHighPcThreshold(pc_star);
        params.setLowSwThreshold(sw_star);
        params.setLowSwSlope(dpcdsw_star);

        ////////////////////////////////////////////////////////////
        //! Precompute regularization for high wetting saturations
        ////////////////////////////////////////////////////////////

        // convert the effective saturation to an absolute one
        // which we actually need later
        const Scalar highSweTh = params.highSweThreshold();
        params.setHighSwThreshold(EffToAbsPolicy::sweToSw(params, highSweTh));
        params.setLowPcThreshold(MaterialLaw::pcRaw(params, highSweTh));
        //! TODO this only works for VanGenuchten
        params.setHighSwSlope(-2*params.lowPcThreshold()/(1.0 - highSweTh));
        const Scalar mTh = MaterialLaw::dpc_dsweRaw(params, highSweTh);
        params.setHighSwSpline(highSweTh, 1.0,
                               params.lowPcThreshold(), 0.0,
                               mTh, params.highSwSlope());

    }

    ///////////////////////////////////////////////////////////
    //! Getter methods to be called by regularization routines
    ///////////////////////////////////////////////////////////

    Scalar ovenDryPc() const
    { return ovenDryPc_; }

    Scalar highPcThreshold() const
    { return highPcThreshold_; }

    Scalar lowSwThreshold() const
    { return lowSwThreshold_; }

    Scalar lowSwSlope() const
    { return thresholdSlope_; }


    Scalar highSwThreshold() const
    { return highSwThreshold_; }

    Scalar highSweThreshold() const
    { return highSweThreshold_; }

    Scalar lowPcThreshold() const
    { return lowPcThreshold_; }

    Scalar highSwSlope() const
    { return highSwSlope_; }

    const Spline<Scalar>& highSwSpline() const
    { return highSwSpline_; }

private:
    //! The user only sets the oven dry pc
    void setLowSwThreshold(const Scalar sw)
    { lowSwThreshold_ = sw; }

    void setLowSwSlope(const Scalar slope)
    { thresholdSlope_ = slope; }

    void setHighPcThreshold(const Scalar pc)
    { highPcThreshold_ = pc; }



    void setHighSwThreshold(const Scalar sw)
    { highSwThreshold_ = sw; }

    void setHighSwSlope(const Scalar slope)
    { highSwSlope_ = slope; }

    void setHighSwSpline(const Scalar x0, const Scalar x1,
                         const Scalar y0, const Scalar y1,
                         const Scalar m0, const Scalar m1)
    { highSwSpline_.set(x0, x1, y0, y1, m0, m1); }

    void setLowPcThreshold(const Scalar pc)
    { lowPcThreshold_ = pc; }

    Scalar ovenDryPc_;
    Scalar highPcThreshold_;
    Scalar lowSwThreshold_;
    Scalar thresholdSlope_;

    Scalar highSwThreshold_;
    Scalar highSweThreshold_;
    Scalar lowPcThreshold_;
    Scalar highSwSlope_;
    Spline<Scalar> highSwSpline_;
};

/*!\ingroup fluidmatrixinteractionslaws
 * \brief Regularization of dry region due to Webb (2000)
 */
template <class Scalar>
class WebbRegularization
{
public:
    using Params = WebbRegularizationParams<Scalar>;

    ////////////////////////////////////////////////////////////
    //! Computing regularized capillary pressures
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low saturations
     */
    template<class ParamType>
    static bool underLowSwThreshold_pc(const ParamType& params, const Scalar sw)
    { return sw < params.lowSwThreshold(); }

    /*!
     * \brief Threshold for high saturations -> no regularization
     */
    template<class ParamType>
    static bool overHighSwThreshold_pc(const ParamType& params, const Scalar sw)
    { return sw > params.highSwThreshold(); }

    /*!
     * \brief Regularized curve for low saturations
     */
    template<class ParamType>
    static Scalar pcForLowSw(const ParamType& params, const Scalar sw)
    {
        return std::pow(10, params.lowSwSlope()*(sw - params.lowSwThreshold()))*params.highPcThreshold();
    }

    /*!
     * \brief Regularized curve for high saturations -> no regularization
     */
    template<class ParamType>
    static Scalar pcForHighSw(const ParamType& params, const Scalar sw)
    {
        const Scalar swe = ParamType::EffToAbsPolicy::swToSwe(params, sw);
        if (swe < 1.0)
            return params.highSwSpline().eval(swe);
        else
            return params.highSwSlope()*(swe - 1.0);
    }

    ////////////////////////////////////////////////////////////
    //! Computing regularized saturations
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low capilary pressures -> no regularization
     */
    template<class ParamType>
    static bool underLowPcThreshold_sw(const ParamType& params, const Scalar pc)
    { return pc < params.lowPcThreshold(); }

    /*!
     * \brief Threshold for high capilary pressures
     */
    template<class ParamType>
    static bool overHighPcThreshold_sw(const ParamType& params, const Scalar pc)
    { return pc > params.highPcThreshold(); }

    /*!
     * \brief Regularized curve for low capilary pressures -> no regularization
     */
    template<class ParamType>
    static Scalar swForLowPc(const ParamType& params, const Scalar pc)
    {
        // check if we are out of the bounds of the raw curve
        if (pc <= 0)
        {
            // for swThHigh = 1.0 the slope would get infinity
            // swThHigh > 1.0 are not sensible threshold values
            // setting swThHigh = 1.0 is a way to disable regularization
            if (params.highSweThreshold() > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return 1.0;

            // otherwise we invert the straight line
            else
                return pc/params.highSwSlope() + 1.0;
        }

        // in physical range but over regularization threshold
        // invert the spline
        else
            return params.highSwSpline().intersectInterval(params.highSweThreshold(), 1.0, 0, 0, 0, pc);
    }

    /*!
     * \brief Regularized curve for high capilary pressures
     */
    template<class ParamType>
    static Scalar swForHighPc(const ParamType& params, const Scalar pc)
    {
        return params.lowSwThreshold() - std::log(params.highPcThreshold()/pc)/(std::log(10)*params.lowSwSlope());
    }

    ////////////////////////////////////////////////////////////
    //! Computing regularized relative permeabilities
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low saturations
     */
    template<class ParamType>
    static constexpr bool underLowSwThreshold_krn(const ParamType& params, const Scalar sw)
    { return false; }

    /*!
     * \brief Threshold for high saturations -> no regularization
     */
    template<class ParamType>
    static constexpr bool overHighSwThreshold_krw(const ParamType& params, const Scalar sw)
    { return false; }

    /*!
     * \brief Regularized curve for high saturations -> no regularization
     */
    template<class ParamType>
    static Scalar krwForHighSw(const ParamType& params, const Scalar sw)
    { return ParamType::MaterialLaw::krwRaw(params, ParamType::EffToAbsPolicy::swToSwe(params, sw)); }

    /*!
     * \brief Regularized curve for low saturations -> no regularization
     */
    template<class ParamType>
    static Scalar krnForLowSw(const ParamType& params, const Scalar sw)
    { return ParamType::MaterialLaw::krnRaw(params, ParamType::EffToAbsPolicy::swToSwe(params, sw)); }

};

} // end namespace Dumux

#endif
