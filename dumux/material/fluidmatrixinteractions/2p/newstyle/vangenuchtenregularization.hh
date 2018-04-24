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
#ifndef DUMUX_VANGENUCHTEN_REGULARISATION_HH
#define DUMUX_VANGENUCHTEN_REGULARISATION_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

/*!
 * \brief Parameters that are necessary for the \em regularization of
 *        the VanGenuchten material law. For low saturations we continue linearly.
 *        For high saturation we use a spline with linear continuation.
 * \ingroup fluidmatrixinteractionsparams
 */
template<class Scalar>
class VanGenuchtenRegularizationParams
{
public:
    VanGenuchtenRegularizationParams()
    {
        setLowSweThreshold_pc(0.01);
        setHighSweThreshold_pc(0.99);
        setLowSweThreshold_krn(0.1);
        setHighSweThreshold_krw(0.9);
    }

    ////////////////////////////////////////////////////////
    //! Setter methods to be called by the user
    ////////////////////////////////////////////////////////

    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setLowSweThreshold_pc(const Scalar swe)
    { lowSweThreshold_pc_ = swe; }

    //! Set the swe threshold over which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setHighSweThreshold_pc(const Scalar swe)
    { highSweThreshold_pc_ = swe; }

    //! Set the swe threshold under which krn-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setLowSweThreshold_krn(const Scalar swe)
    { lowSweThreshold_krn_ = swe; }

    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setHighSweThreshold_krw(const Scalar swe)
    { highSweThreshold_krw_ = swe; }

    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    DUNE_DEPRECATED_MSG("Use setLowSweThreshold_pc(swe) instead!")
    void setPcLowSw(const Scalar swe)
    { lowSweThreshold_pc_ = swe; }

    //! Set the swe threshold over which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    DUNE_DEPRECATED_MSG("Use setHighSweThreshold_pc(swe) instead!")
    void setPcHighSw(const Scalar swe)
    { highSweThreshold_pc_ = swe; }

    //! Set the swe threshold under which krn-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    DUNE_DEPRECATED_MSG("Use setLowSweThreshold_krn(swe) instead!")
    void setKrnLowSw(const Scalar swe)
    { lowSweThreshold_krn_ = swe; }

    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    DUNE_DEPRECATED_MSG("Use setHighSweThreshold_krw(swe) instead!")
    void setKrwHighSw(const Scalar swe)
    { highSweThreshold_krw_ = swe; }

    //! This is a free functions, being a friend we can access private members
    template<class Params>
    friend void initRegularizationParams(Params& params)
    {
        using MaterialLaw = typename Params::MaterialLaw;
        using EffToAbsPolicy = typename Params::EffToAbsPolicy;

        ///////////////////////////////////////////////////////////////
        //! Precompute pc regularization for small wetting saturations
        ///////////////////////////////////////////////////////////////
        const Scalar lowSweTh_pc = params.lowSweThreshold_pc();
        params.setLowSwThreshold_pc(EffToAbsPolicy::sweToSw(params, lowSweTh_pc));
        params.setHighPcThreshold_sw(MaterialLaw::pcRaw(params, lowSweTh_pc));
        params.setLowSwSlope_pc(MaterialLaw::dpc_dsweRaw(params, lowSweTh_pc));

        //////////////////////////////////////////////////////////////
        //! Precompute pc regularization for high wetting saturations
        //////////////////////////////////////////////////////////////

        // convert the effective saturation to an absolute one
        // which we actually need later
        const Scalar highSweTh_pc = params.highSweThreshold_pc();
        params.setHighSwThreshold_pc(EffToAbsPolicy::sweToSw(params, highSweTh_pc));
        params.setLowPcThreshold_sw(MaterialLaw::pcRaw(params, highSweTh_pc));
        params.setHighSwSlope_pc(-2*params.lowPcThreshold_sw()/(1.0 - highSweTh_pc));
        const Scalar mTh = MaterialLaw::dpc_dsweRaw(params, highSweTh_pc);
        params.setHighSwSpline_pc(highSweTh_pc, 1.0,
                                  params.lowPcThreshold_sw(), 0.0,
                                  mTh, params.highSwSlope_pc());

        //////////////////////////////////////////////////////////////
        //! Precompute krw regularization for high wetting saturations
        //////////////////////////////////////////////////////////////
        const Scalar highSweTh_krw = params.highSweThreshold_krw();
        params.setHighSwThreshold_krw(EffToAbsPolicy::sweToSw(params, highSweTh_krw));
        params.setHighSwSpline_krw(highSweTh_krw, 1.0,
                                   MaterialLaw::krwRaw(params, highSweTh_krw), 1.0,
                                   MaterialLaw::dkrw_dsweRaw(params, highSweTh_krw), 0.0);

        //////////////////////////////////////////////////////////////
        //! Precompute krn regularization for low wetting saturations
        //////////////////////////////////////////////////////////////
        const Scalar lowSweTh_krn = params.lowSweThreshold_krn();
        params.setLowSwThreshold_krn(EffToAbsPolicy::sweToSw(params, lowSweTh_krn));
        params.setLowSwSpline_krn(0.0, lowSweTh_krn,
                                  1.0, MaterialLaw::krnRaw(params, lowSweTh_krn),
                                  0.0, MaterialLaw::dkrn_dsweRaw(params, lowSweTh_krn));
    }

    ///////////////////////////////////////////////////////////
    //! Getter methods to be called by regularization routines
    ///////////////////////////////////////////////////////////

    //! Thresholds with effective saturations

    Scalar highSweThreshold_pc() const
    { return highSweThreshold_pc_; }

    Scalar lowSweThreshold_pc() const
    { return lowSweThreshold_pc_; }

    Scalar highSweThreshold_krw() const
    { return highSweThreshold_krw_; }

    Scalar lowSweThreshold_krn() const
    { return lowSweThreshold_krn_; }

    //! Thresholds with absolute saturations

    Scalar highPcThreshold_sw() const
    { return highPcThreshold_sw_; }

    Scalar lowPcThreshold_sw() const
    { return lowPcThreshold_sw_; }

    Scalar lowSwThreshold_pc() const
    { return lowSwThreshold_pc_; }

    Scalar highSwThreshold_pc() const
    { return highSwThreshold_pc_; }

    Scalar lowSwThreshold_krn() const
    { return lowSwThreshold_krn_; }

    Scalar highSwThreshold_krw() const
    { return highSwThreshold_krw_; }

    //! Slopes and splines

    Scalar lowSwSlope_pc() const
    { return lowSwSlope_pc_; }

    Scalar highSwSlope_pc() const
    { return highSwSlope_pc_; }

    const Spline<Scalar>& highSwSpline_pc() const
    { return highSwSpline_pc_; }

    const Spline<Scalar>& highSwSpline_krw() const
    { return highSwSpline_krw_; }

    const Spline<Scalar>& lowSwSpline_krn() const
    { return lowSwSpline_krn_; }

private:
    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setLowSwThreshold_pc(const Scalar sw)
    { lowSwThreshold_pc_ = sw; }

    //! Set the swe threshold over which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setHighSwThreshold_pc(const Scalar sw)
    { highSwThreshold_pc_ = sw; }

    //! Set the swe threshold under which krn-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setLowSwThreshold_krn(const Scalar sw)
    { lowSwThreshold_krn_ = sw; }

    //! Set the swe threshold under which pc-sw will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setHighSwThreshold_krw(const Scalar sw)
    { highSwThreshold_krw_ = sw; }

    //! Set the pc threshold under which sw-pc will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setLowPcThreshold_sw(const Scalar pc)
    { lowPcThreshold_sw_ = pc; }

    //! Set the pc threshold over which sw-pc will be regularized
    //! If you change this you have to call initRegularizationParams again
    void setHighPcThreshold_sw(const Scalar pc)
    { highPcThreshold_sw_ = pc; }


    //! Regularized pc curve for high sw
    void setHighSwSlope_pc(const Scalar slope)
    { highSwSlope_pc_ = slope; }

    //! Regularized pc curve for low sw
    void setLowSwSlope_pc(const Scalar slope)
    { lowSwSlope_pc_ = slope; }

    //! Regularized pc curve for high sw
    void setHighSwSpline_pc(const Scalar x0, const Scalar x1,
                            const Scalar y0, const Scalar y1,
                            const Scalar m0, const Scalar m1)
    { highSwSpline_pc_.set(x0, x1, y0, y1, m0, m1); }

    //! Regularized krw curve for high sw
    void setHighSwSpline_krw(const Scalar x0, const Scalar x1,
                             const Scalar y0, const Scalar y1,
                             const Scalar m0, const Scalar m1)
    { highSwSpline_krw_.set(x0, x1, y0, y1, m0, m1); }

    //! Regularized krw curve for high sw
    void setLowSwSpline_krn(const Scalar x0, const Scalar x1,
                            const Scalar y0, const Scalar y1,
                            const Scalar m0, const Scalar m1)
    { lowSwSpline_krn_.set(x0, x1, y0, y1, m0, m1); }


    //! Thresholds with effective saturations
    Scalar lowSweThreshold_pc_;
    Scalar highSweThreshold_pc_;
    Scalar lowSweThreshold_krn_;
    Scalar highSweThreshold_krw_;

    //! Thresholds with absolute saturations
    Scalar lowSwThreshold_pc_;
    Scalar highSwThreshold_pc_;
    Scalar lowSwThreshold_krn_;
    Scalar highSwThreshold_krw_;
    Scalar lowPcThreshold_sw_;
    Scalar highPcThreshold_sw_;

    //! Slopes and splines
    Scalar highSwSlope_pc_;
    Scalar lowSwSlope_pc_;
    Spline<Scalar> highSwSpline_pc_;
    Spline<Scalar> highSwSpline_krw_;
    Spline<Scalar> lowSwSpline_krn_;
};

/*!\ingroup fluidmatrixinteractionslaws
 * \brief Regularization dry and wet region with linear extensions
 */
template <class Scalar>
class VanGenuchtenRegularization
{
public:
    using Params = VanGenuchtenRegularizationParams<Scalar>;

    ////////////////////////////////////////////////////////////
    //! Computing regularized capillary pressures
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low saturations
     */
    template<class ParamType>
    static bool underLowSwThreshold_pc(const ParamType& params, const Scalar sw)
    { return sw < params.lowSwThreshold_pc(); }

    /*!
     * \brief Threshold for high saturations -> no regularization
     */
    template<class ParamType>
    static bool overHighSwThreshold_pc(const ParamType& params, const Scalar sw)
    { return sw > params.highSwThreshold_pc(); }

    /*!
     * \brief Regularized curve for low saturations
     */
    template<class ParamType>
    static Scalar pcForLowSw(const ParamType& params, const Scalar sw)
    {
        const Scalar swe = ParamType::EffToAbsPolicy::swToSwe(params, sw);
        const Scalar lowSweTh_pc = params.lowSweThreshold_pc();
        return params.highPcThreshold_sw()
               + params.lowSwSlope_pc()*(swe - lowSweTh_pc);
    }

    /*!
     * \brief Regularized curve for high saturations -> no regularization
     */
    template<class ParamType>
    static Scalar pcForHighSw(const ParamType& params, const Scalar sw)
    {
        const Scalar swe = ParamType::EffToAbsPolicy::swToSwe(params, sw);
        if (swe < 1.0)
            return params.highSwSpline_pc().eval(swe);
        else
            return params.highSwSlope_pc()*(swe - 1.0);
    }

    ////////////////////////////////////////////////////////////
    //! Computing regularized saturations
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low capilary pressures -> no regularization
     */
    template<class ParamType>
    static bool underLowPcThreshold_sw(const ParamType& params, const Scalar pc)
    { return pc < params.lowPcThreshold_sw(); }

    /*!
     * \brief Threshold for high capilary pressures
     */
    template<class ParamType>
    static bool overHighPcThreshold_sw(const ParamType& params, const Scalar pc)
    { return pc > params.highPcThreshold_sw(); }

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
            if (params.highSweThreshold_pc() > 1.0 - std::numeric_limits<Scalar>::epsilon())
                return 1.0;

            // otherwise we invert the straight line
            else
                return pc/params.highSwSlope_pc() + 1.0;
        }

        // in physical range but over regularization threshold
        // invert the spline
        else
            return params.highSwSpline_pc().intersectInterval(params.highSweThreshold_pc(), 1.0, 0, 0, 0, pc);
    }

    /*!
     * \brief Regularized curve for high capilary pressures
     */
    template<class ParamType>
    static Scalar swForHighPc(const ParamType& params, const Scalar pc)
    {
        const Scalar highPcTh = params.highPcThreshold_sw();
        return (pc - highPcTh)/params.lowSwSlope_pc() + params.lowSweThreshold_pc();
    }

    ////////////////////////////////////////////////////////////
    //! Computing regularized relative permeabilities
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Threshold for low saturations
     */
    template<class ParamType>
    static bool underLowSwThreshold_krn(const ParamType& params, const Scalar sw)
    { return sw < params.lowSwThreshold_krn(); }

    /*!
     * \brief Threshold for high saturations
     */
    template<class ParamType>
    static bool overHighSwThreshold_krw(const ParamType& params, const Scalar sw)
    { return sw > params.highSwThreshold_krw(); }

    /*!
     * \brief Regularized curve for high saturations
     */
    template<class ParamType>
    static Scalar krwForHighSw(const ParamType& params, const Scalar sw)
    {
        const Scalar swe = ParamType::EffToAbsPolicy::swToSwe(params, sw);
        if (swe > 1.0 - std::numeric_limits<Scalar>::epsilon())
            return 1.0;
        // we are between swe=1 and swe=highThreshold
        else
            return params.highSwSpline_krw().eval(swe);
    }

    /*!
     * \brief Regularized curve for low saturations
     */
    template<class ParamType>
    static Scalar krnForLowSw(const ParamType& params, const Scalar sw)
    {
        const Scalar swe = ParamType::EffToAbsPolicy::swToSwe(params, sw);
        if (swe < 0.0)
            return 1.0;
        // we are between swe=0 and swe=lowThreshold
        else
            return params.lowSwSpline_krn().eval(swe);
    }

};

} // end namespace Dumux

#endif
