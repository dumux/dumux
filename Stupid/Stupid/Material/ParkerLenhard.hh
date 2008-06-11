/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file ParkerLenhard.hh
 * \brief Template implementing Parker-Lenhard capillary pressure hystersis.
 */
#ifndef PARKER_LENHARD_HH
#define PARKER_LENHARD_HH

#include <Stupid/Auxilary/Apis.hh>

#include <Stupid/Material/TwophaseSat.hh>
#include <Stupid/Material/VanGenuchten.hh>
#include <Stupid/Material/ParkerLenhardState.hh>

#include <math.h>
#include <assert.h>

#include <algorithm>

//#define PL_STRICT_ASSERTS

#ifdef PL_STRICT_ASSERTS
  // we use a macro here, since if the assertation fails at run time,
  // we want to see at which line it happened
  #define ASSERT_RANGE(VAL, LEFT, RIGHT) assert(LEFT <= VAL && VAL <= RIGHT);
#else // !PL_STRICT_ASSERTS
  template<class Scalar, class Scalar1, class Scalar2>
  void ASSERT_RANGE(Scalar &var, Scalar1 left, Scalar2 right)
  { var = std::max(Scalar(left), std::min(Scalar(right), var)); }
#endif

namespace Stupid
{
    /*!
     * \internal
     * \brief Represents a scanning curve.
     *
     * The class has pointers to the scanning curves
     * with higher and lower loop number, this saving
     * the history of the imbibitions and drainages.
     */
    template <class ScalarT>
    class PLScanningCurve
    {
    public:
        typedef ScalarT Scalar;

        /*!
         * \brief Constructs main imbibtion curve.
         *
         * Further scanning curves can be added with
         * setNext.
         */
        PLScanningCurve()
        {
            _loopNum = 0;
            _prev = new PLScanningCurve(NULL, // prev
                                        this, // next
                                        -1,   // loop number
                                        0,    // Sw
                                        1e12, // pC
                                        0,    // Sw_app
                                        0,    // SwMic
                                        0);   // SwMdc
            _next = NULL;

            _Swe = 1.0;
            _pC  = 0.0;
            _Sw_app = 1.0;
            _SwMic = 1.0;
            _SwMdc = 1.0;
        }

    protected:
        PLScanningCurve(PLScanningCurve *prev,
                        PLScanningCurve *next,
                        int    loopN,
                        Scalar Swe,
                        Scalar pC,
                        Scalar Sw_app,
                        Scalar SwMic,
                        Scalar SwMdc)
        {
            _prev = prev;
            _next = next;
            _loopNum = loopN;
            _Swe = Swe;
            _pC  = pC;
            _Sw_app = Sw_app;
            _SwMic = SwMic;
            _SwMdc = SwMdc;
        }

    public:
        /*!
         * \brief Destructor. After it was called
         *        all references to the next() curve are
         *        invalid!
         */
        ~PLScanningCurve()
        {
                if (_loopNum == 0)
                    delete _prev;
                if (_loopNum >= 0)
                    delete _next;
        }

        /*!
         * \brief Return the previous scanning curve, i.e. the curve
         *        with one less reversal than the current one.
         */
        PLScanningCurve *prev() const
        { return _prev; }

        /*!
         * \brief Return the next scanning curve, i.e. the curve
         *        with one more reversal than the current one.
         */
        PLScanningCurve *next() const
        { return _next; }

        /*!
         * \brief Set the next scanning curve.
         *
         * Next in the sense of the number of reversals
         * from imbibition to drainage or vince versa. If this
         * curve already has a list of next curves, it is
         * deleted and thus forgotten.
         */
        void setNext(Scalar Swe,
                     Scalar pC,
                     Scalar Sw_app,
                     Scalar SwMic,
                     Scalar SwMdc)
            {
                // if _next is NULL, delete does nothing, so
                // this is valid!!
                delete _next;

                _next = new PLScanningCurve(this, // prev
                                            NULL, // next
                                            loopNum() + 1,
                                            Swe,
                                            pC,
                                            Sw_app,
                                            SwMic,
                                            SwMdc);
            }

        /*!
         * \brief Erase all curves with higher loop number than the
         *        this one.
         */
        void eraseHigher()
            {
                delete _next;
                _next = NULL;
            }
        
        /*!
         * \brief Returns true iff the given effective saturation
         *        Swei is within the scope of the curve, i.e.
         *        whether Swei is part of the curve's
         *        domain and the curve thus applies to Swi.
         */
        bool isValidAt_Swe(Scalar Swei)
        {
            if (isImbib())
                // for inbibition the given saturation
                // must be between the start of the
                // current imbibition and the the start
                // of the last drainage
                return Swe() < Swei && Swei < _prev->Swe();
            else
                // for drainage the given saturation
                // must be between the start of the
                // last imbibition and the start
                // of the current drainage
                return _prev->Swe() < Swei && Swei < Swe();
        };

        /*!
         * \brief Returns true iff a given capillary pressure
         *        pC is within the scope of the curve, i.e.
         *        whether pC is part of the curve's
         *        image and the curve thus applies to pC.
         */
        bool isValidAt_pC(Scalar pCi)
        {
            if (isImbib())
                return _prev->pC() < pCi  && pCi < pC();
            else
                return pC() < pCi && pCi < _prev->pC();
        };

        /*!
         * \brief Returns true iff the scanning curve is a
         *        imbibition curve.
         */
        bool isImbib()
        { return loopNum()%2 == 1; }

        /*!
         * \brief Returns true iff the scanning curve is a
         *        drainage curve.
         */
        bool isDrain()
        { return !isImbib(); }

        /*!
         * \brief The loop number of the scanning curve.
         *
         * The MDC is 0, PISC is 1, PDSC is 2, ...
         */
        int loopNum()
        { return _loopNum; }

        /*!
         * \brief Effective saturation at the
         *        last reversal point.
         */
        Scalar Swe() const
        { return _Swe; }

        /*!
         * \brief Capillary pressure at the last reversal point.
         */
        Scalar pC() const
        { return _pC; }

        /*!
         * \brief Apparent saturation at the
         *        last reversal point.
         */
        Scalar Sw_app() const
        { return _Sw_app; }

        /*!
         * \brief Apparent saturation of the last reversal point on
         *        the pressure MIC.
         */
        Scalar SwMic()
        { return _SwMic; }

        /*!
         * \brief Apparent saturation of the last reversal point on
         *        the pressure MDC.
         */
        Scalar SwMdc()
        { return _SwMdc; }

    private:
        PLScanningCurve *_prev;
        PLScanningCurve *_next;

        int    _loopNum;

        Scalar _Swe;
        Scalar _pC;
        Scalar _Sw_app;

        Scalar _SwMdc;
        Scalar _SwMic;
    };

    /*!
     * \ingroup material
     * \brief Implements the Parker-Lenhard twophase
     *        p_c-Sw hysteresis model. This class adheres to the twophase
     *        capillary pressure API.
     */
    template <class StateT, class CapPressureT>
    class ParkerLenhard
    {
    public:
        typedef CapPressureT CapPressure;
        typedef StateT State;

        typedef typename State::Scalar Scalar;
        typedef Stupid::PLScanningCurve<Scalar> ScanningCurve;
        typedef Stupid::TwophaseSat<State> TwophaseSat;

        /*!
         * \brief Resets the hysteresis model to the
         *        initial state on the main drainage curve
         */
        static void reset(State &state)
        {
            Api::require<Api::ParkerLenhardState>(state);
            Api::require<Api::TwophaseSatParams>(state);

            delete state.mdc(); // this will work even if _mdc == NULL!
            state.setMdc(new ScanningCurve());
            state.setCsc(state.mdc());
            state.setPisc(NULL);
            state.setSnrei(0.0);
        };

        /*!
         * \brief Returns true iff the given absolute saturation
         *        is viable (i.e. larger than the residual
         *        saturation of the wetting phase but smaller
         *        than the maximum saturation of the current
         *        PISC)
         */
        static bool isValidSw(const State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

            Scalar Swe = TwophaseSat::Swe(state, Sw);

            return 0 <= Swe && Swe <= 1. - state.Snrei();
        }

        /*!
         * \brief Set the current absolute saturation for the
         *        current timestep
         */
        static void updateState(State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardState>(state);
            Api::require<Api::TwophaseSatParams>(state);

            if (Sw > 1 - 1e-5) {
                // if the absolute saturation is almost 1,
                // it means that we're back to the beginning
                reset(state);
                return;
            }

            Scalar Swe = TwophaseSat::Swe(state, Sw);
            if (state.pisc() == NULL && Swe > 0.97) {
                // for numeric reasons don't start a PISC if the
                // wetting saturation is still higher than 97%
                return;
            }
            ASSERT_RANGE(Swe, 0, 1);
            
            // find the loop number which corrosponds to the
            // given effective saturation
            ScanningCurve *curve = _findScanningCurve_Swe(state, Swe);
            
            // make sure that the distance between two consecutive
            // curves is large enough in order to avoid too steep
            // capillary pressure curves
            if (fabs(curve->Swe() - Swe) < 10e-2)
            {
                if (curve != state.mdc()) {
                    curve->eraseHigher();
                }
                return;
            }
            else if (curve != state.mdc() && 
                     fabs(curve->prev()->Swe() - Swe) < 10e-2)
            {
                return;
            }
            

            Scalar Sw_app = _Swapp(state, Swe);

            // calculate the apparent saturation on the MIC and MDC
            // which yield the same capillary pressure as the
            // Sw at the current scanning curve
            Scalar pc = pC(state, Sw);
            Scalar Sw_mic = CapPressure::Sw(state.micParams(), pc);
            Scalar Sw_mdc = CapPressure::Sw(state.mdcParams(), pc);
            curve->setNext(Swe, pc, Sw_app,
                           Sw_mic, Sw_mdc);

            state.setCsc(curve);

            // if we're back on the MDC, we also have a new PISC!
            if (state.csc() == state.mdc()) {
                state.setPisc(state.mdc()->next());
                _calcSnrei(state, Swe);
            }
        }


        /*!
         * \brief Returns the capillary pressure dependend on
         *        the _absolute_ saturation of the wetting phase.
         */
        static Scalar pC(const State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

            Scalar Swe = TwophaseSat::Swe(state, Sw);

            // if the effective saturation is smaller than zero we let
            // the underlying material law decide what to do.
            if (Swe < 0) {
                return CapPressure::pC(state.mdcParams(), Swe);
            }

            // calculate the current apparent saturation
            ScanningCurve *sc = _findScanningCurve_Swe(state, Swe);
            Scalar Sw_app = _Swapp(state, Swe);

            // if the apparent saturation exceeds the 'legal' limits,
            // we also the underlying material law decide what to do.
            if (Sw_app > 1) {
                return CapPressure::pC(state.mdcParams(), Sw_app);
            }

            // put the effective saturation into the capillary pressure model
            Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
            if (sc->isImbib()) {
                Scalar SwMic =
                    pos * (sc->prev()->SwMic() - sc->SwMic())
                    + sc->SwMic();

                return CapPressure::pC(state.micParams(), SwMic);
            }
            else { // sc->isDrain()
                Scalar SwMdc =
                    pos*(sc->prev()->SwMdc() - sc->SwMdc())
                    + sc->SwMdc();

                return CapPressure::pC(state.mdcParams(), SwMdc);
            }
        }

        /*!
         * \brief Returns the absolute saturation to a given on the
         *        capillary pressure.
         */
        static Scalar Sw(const State &state, Scalar pC)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

            assert(pC >= 0);

            // find the relevant scanning curve based on
            // the given capillary pressure
            ScanningCurve *sc = _findScanningCurve_pC(state, pC);
            // if we're on the MDC, we're already done...
            if (sc == state.mdc()) {
                Scalar Swe = CapPressure::Sw(state.mdcParams(), pC);
                return TwophaseSat::Sw(state, Swe);
            }

            // like in the forward direction, we first calculate
            // the relative position on the scanning curve sc
            Scalar pos;
            if (sc->isImbib()) {
                Scalar SwMic = CapPressure::Sw(state.micParams(), pC);
                pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
            }
            else {
                Scalar SwMdc = CapPressure::Sw(state.mdcParams(), pC);
                pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
            }

            // now we can calculate the apparent saturation
            Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();

            // which can finally converted into an absolute saturation
            return TwophaseSat::Sw(state, _SweFromSwapp(state, Sw_app));
        }

        /*!
         * \brief The derivative of the capillary pressure regarding
         *        the absolute saturation.
         */
        static Scalar dpC_dSw(const State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

            Scalar Swe = TwophaseSat::Swe(state, Sw);

            // if the effective saturation exceeds 1, we have a
            // problem, but we let the underlying capillary pressure
            // model deal with it. it might apply a regularization,
            // or it might crash ;)
            if (Swe > 1)
                return CapPressure::dpC_dSw(state.mdcParams(), Swe);

            // calculate the current apparent saturation
            ScanningCurve *sc = _findScanningCurve_Swe(state, Swe);
            Scalar Sw_app = _Swapp(state, Swe);
            ASSERT_RANGE(Sw_app, 0, 1);

            // calculate the derivative of the linear interpolation parameter
            // in regard to the apparent saturation
            Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
            Scalar dpos_dSwapp = 1/(sc->prev()->Sw_app() - sc->Sw_app());

            ASSERT_RANGE(pos, 0, 1);

            if (sc->isImbib()) {
                Scalar SwMic =
                    pos * (sc->prev()->SwMic() - sc->SwMic())
                    + sc->SwMic();
                // the factor behind the pos variable is a constant
                Scalar dSwMic_dSwapp = dpos_dSwapp *
                                       (sc->prev()->SwMic() - sc->SwMic());
                ASSERT_RANGE(SwMic, 0, 1);

                // inner times outer derivative (-> chain rule)
                return CapPressure::dpC_dSw(state.micParams(), SwMic)*
                          dSwMic_dSwapp*
                          _dSwapp_dSwe(state, sc) *
                          TwophaseSat::dSwe_dSw(state);
            }
            else { // sc->isDrain()
                Scalar SwMdc =
                    pos*(sc->prev()->SwMdc() - sc->SwMdc())
                    + sc->SwMdc();
                // the factor behind the pos variable is a constant
                Scalar dSwMdc_dSwapp = dpos_dSwapp *
                                       (sc->prev()->SwMdc() - sc->SwMdc());
                ASSERT_RANGE(SwMdc, 0, 1);

                // inner times outer derivative (-> chain rule)
                return CapPressure::dpC_dSw(state.mdcParams(), SwMdc)*
                         dSwMdc_dSwapp *
                         _dSwapp_dSwe(state, sc) *
                         TwophaseSat::dSwe_dSw(state);
            }
        }

        /*!
         * \brief The derivative of the absolute saturation regarding
         *        the capillary pressure.
         */
        static Scalar dSw_dpC (const State &state, Scalar pC)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

#ifdef PL_STRICT_ASSERTS
            assert(pC >= 0);
#else
            if (pC < 0) pC = 0;
#endif

            // find the relevant scanning curve based on
            // the given capillary pressure
            ScanningCurve *sc = _findScanningCurve_pC(state, pC);
            // if we're on the MDC, we're already done...
            if (sc == state.mdc()) {
                return CapPressure::dSw_dpC(state.mdcParams(), pC)*TwophaseSat::dSw_dSwe(state);
            }

            // like in the forward direction, we first calculate
            // the relative position on the scanning curve sc
//            Scalar pos;
            Scalar dpos_dSwMc;
            Scalar dSwMc_dpC;
            if (sc->isImbib()) {
//                Scalar SwMic = state.mic().Sw(pC);
//                pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
                dpos_dSwMc = 1/(sc->prev()->SwMic() - sc->SwMic());
                dSwMc_dpC = CapPressure::dSw_dpC(state.micParams(), pC);
            }
            else {
//                Scalar SwMdc = state.mdc().Sw(pC);
//                pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
                dpos_dSwMc = 1/(sc->prev()->SwMdc() - sc->SwMdc());
                dSwMc_dpC = CapPressure::dSw_dpC(state.mdcParams(), pC);
            }

            // now we can calculate the apparent saturation
//            Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();
            Scalar dSwapp_dpos = sc->prev()->Sw_app() - sc->Sw_app();
            Scalar dSwapp_dSwMc = dSwapp_dpos*dpos_dSwMc;

            return dSwMc_dpC*dSwapp_dSwMc*
                   _dSwe_dSwapp(state, sc)*TwophaseSat::dSw_dSwe(state);
        }

        /*!
         * \brief The relative permeability for the wetting phase of
         *        the medium.
         */
        static Scalar krw(const State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);

            ASSERT_RANGE(Sw, 0, 1);

            // for the effective permeability we use play-type
            // hystersis. (That's because it's basically impossible
            // to invert the krw-Sw realtion effectivly.)
            Scalar Swe = TwophaseSat::Swe(state, Sw);

            return CapPressure::krw(state.mdcParams(), Swe);
            
#if 0 // TODO: saturation-permebility hysteresis
            Scalar Sw_app = _Swapp(state, Swe);
            ASSERT_RANGE(Sw_app, 0, 1);

            if (state.pisc() && Swe > state.csc()->Swe()) {
                return CapPressure::krw(state.micParams(), Sw_app);
            }
            else { // sc->isDrain()
                return CapPressure::krw(state.mdcParams(), Sw_app);
            }
#endif
        }

        /*!
         * \brief The relative permeability for the non-wetting phase
         *        of the params.
         */
        static Scalar krn(const State &state, Scalar Sw)
        {
            Api::require<Api::ParkerLenhardParams>(state);
            Api::require<Api::TwophaseSatParams>(state);
            
            ASSERT_RANGE(Sw, 0, 1);

            // for the effective permeability we use play-type
            // hystersis. (That's because it's basically impossible
            // to invert the krw-Sw realtion effectivly.)
            Scalar Swe = TwophaseSat::Swe(state, Sw);

            return CapPressure::krn(state.mdcParams(), Swe);

#if 0 // TODO: saturation-permebility hysteresis
            Scalar Sw_app = _Swapp(state, Swe);
            ASSERT_RANGE(Sw_app, 0, 1);

            if (state.pisc() && Swe > state.csc()->Swe()) {
                return CapPressure::krn(state.micParams(), Sw_app);
            }
            else { // sc->isDrain()
                return CapPressure::krn(state.mdcParams(), Sw_app);
            }
#endif
        }

        static Scalar SwToSwapp(const State &state, Scalar Sw)
        {
            return _Swapp(state, TwophaseSat::Swe(state, Sw));
        }
    private:
        // find the loop on which the an effective
        // saturation has to be
        static ScanningCurve *_findScanningCurve_Swe(const State &state, Scalar Swe)
        {
            if (state.pisc() == NULL || Swe <= state.pisc()->Swe()) {
                // we don't have a PISC yet, or the effective
                // saturation is smaller than the saturation where the
                // PISC begins. In this case are on the MDC
                return state.mdc();
            }

            // if we have a primary imbibition curve, and our current
            // effective saturation is higher than the beginning of
            // the secondary drainage curve. this means we are on the
            // PISC again.
            if (state.pisc()->next() == NULL ||
                state.pisc()->next()->Swe() < Swe)
            {
                return state.pisc();
            }

            ScanningCurve *curve = state.csc()->next();
            while (true) {
                assert(curve != state.mdc()->prev());
                if (curve->isValidAt_Swe(Swe)) {
                    return curve;
                }
                curve = curve->prev();
            }
        }

        // find the loop on which an capillary pressure belongs to
        static ScanningCurve *_findScanningCurve_pC(const State &state, Scalar pC)
        {
            if (state.mdc()->next() == NULL) {
                // we don't have a PISC yet,
                // so we must be on the MDC
                // (i.e. _csc == _mdc)
                return state.mdc();
            }

            ScanningCurve *curve = state.csc()->next();
            while (true) {
                assert(curve != state.mdc()->prev());
                if (curve->isValidAt_pC(pC))
                    return curve;
                curve = curve->prev();
            }
        }

        // calculate and save _Snrei
        static void _calcSnrei(State &state, Scalar Swei)
        {
            if (state.Snre() == 0.0) {
                state.setSnrei(0.0);
            }
            else {
                // use Land's law
                Scalar R = 1.0/state.Snre() - 1;
                state.setSnrei((1 - Swei)/(1 + R*(1 - Swei)));
            }

            // if the effective saturation is smaller or equal 100%,
            // the current trapped saturation must be smaller than the
            // residual saturation
            assert(Swei > 1 || state.Snrei() <= state.Snre());
        };

        // returns the trapped effective saturation at j
        static Scalar _Snrij(const State &state, Scalar Swej)
        {
            return state.Snrei()*(Swej - state.pisc()->Swe()) /
                                 (1 - state.Snrei() - state.pisc()->Swe());
        };

        // returns the apparent saturation of the
        // wetting phase depending on the effective saturation
        static Scalar _Swapp(const State &state, Scalar Swe)
        {
            if (state.pisc() == NULL || Swe <= state.pisc()->Swe()) {
                // we are on the main drainage curve, i.e.
                // no non-wetting fluid is trapped
                // -> apparent saturation == effective saturation
                return Swe;
            }


            // we are on a imbibition or drainage curve
            // which is not the main drainage curve
            // -> apparent saturation ==
            //    effective saturation + trapped effective saturation
            return Swe + _Snrij(state, Swe);
        };

        // Returns the effective saturation to a given apparent one
        static Scalar _SweFromSwapp(const State &state, Scalar Sw_app)
        {
            if (state.pisc() == NULL || Sw_app <= state.pisc()->Swe()) {
                // we are on the main drainage curve, i.e.
                // no non-wetting fluid is trapped
                // -> apparent saturation == effective saturation
                return Sw_app;
            }

            return (Sw_app*(1 - state.Snrei()
                              - state.pisc()->Swe())
                              + state.pisc()->Swe()*state.Snrei())
                   /(1 - state.pisc()->Swe());
        };

        // returns the derivative of the apparent saturation in
        // regard to the effective one
        static Scalar _dSwapp_dSwe(const State &state, ScanningCurve *sc)
            {
                if (sc == state.mdc())
                    return 1.0;

                return 1 + state.Snrei()/(1 - state.Snrei()
                                            - state.pisc()->Swe());
            }

        // returns the derivative of the apparent saturation in
        // regard to the effective one
        static Scalar _dSwe_dSwapp(const State &state, ScanningCurve *sc)
            {
                if (sc == state.mdc())
                    return 1.0;

                return (1 - state.Snrei() - state.pisc()->Swe())
                     / (1 - state.pisc()->Swe());
            }
    };

#undef ASSERT_RANGE
}; // namespace Stupid

#endif // PARKER_LENHARD_HH
