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
 * \file
 * \brief A simple plotting application for f(Sw) material laws.
 */
#include "StupidConfig.hh"

#include <Stupid/Material/ParkerLenhard.hh>
#include <Stupid/Material/PlayType.hh>

#include <Stupid/Material/RegularizedVanGenuchten.hh>
#include <Stupid/Material/BrooksCorey.hh>

#include <Stupid/Material/TwophaseSat.hh>

#include <Stupid/Auxilary/StateHelperMacros.hh>

#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

using namespace Stupid;

////////////////////
// Some global parameters
////////////////////
typedef long double Scalar;
const   int numSamples = 1500;
const   Scalar Swr = 0.17;
//const   Scalar Swr = 0.0;
const   Scalar Snr = 0.15;
//const   Scalar Snr = 0.0;

template <class CapPressureParamsT>
class MyGlobalState
{
public:
    typedef CapPressureParamsT                  CapPressureParams;
    typedef typename CapPressureParams::Scalar  Scalar;

    MyGlobalState(Scalar Swr,
                    Scalar Snr,
                    const CapPressureParams &imbibParams,
                    const CapPressureParams &drainParams)
        : _micParams(imbibParams), _mdcParams(drainParams)
        {
            Api::require<Api::TwophaseSatState>(*this);
            _globalStateObj = this;

            _Swr = Swr;
            setSnr(Snr);
        }

    static MyGlobalState *get()
        { return _globalStateObj; }

    PROPERTY(CapPressureParams, micParams, setMicParams);
    PROPERTY(CapPressureParams, mdcParams, setMdcParams);

    // we also need to update Snre the Swr is changed, so that
    // we have to define the setter for Swr ourselfs
    PARAMETER(Scalar, Swr);
    void setSwr(Scalar Swr){ _Swr = Swr; _Snre = _Snr/(1 - _Swr); }

    // we also need to update Snre the Snr is changed, so that
    // we have to define the setter for Snr ourselfs
    PARAMETER(Scalar, Snr);
    void setSnr(Scalar Snr){ _Snr = Snr; _Snre = _Snr/(1 - _Swr); }
    PARAMETER(Scalar, Snre);

private:
    static MyGlobalState *_globalStateObj;
};

// Declaration of MyGlobalState::_globalStateObj. This must
// be done outside the class. (looks pretty ugly, eh?)
template <class PCParamsT>
MyGlobalState<PCParamsT> *MyGlobalState<PCParamsT>::_globalStateObj;

///////////// Parker-Lenhard State
template <class GlobalState>
class MyPLState
{
public:
    typedef typename GlobalState::Scalar            Scalar;
    typedef PLScanningCurve<Scalar>                 ScanningCurve;
    typedef typename GlobalState::CapPressureParams CapPressureParams;

    MyPLState()
        {
            Api::require<Api::TwophaseSatParams>(*this);
            Api::require<Api::ParkerLenhardState>(*this);

            _Snrei = 0;
            _mdc = new ScanningCurve();
            _pisc = _csc = NULL;
        }

    ~MyPLState()
        {
            delete _mdc;
        }

    PROXY_PARAMETER(GlobalState::get(), Scalar, Swr);
    PROXY_PARAMETER(GlobalState::get(), Scalar, Snr);
    PROXY_PARAMETER(GlobalState::get(), Scalar, Snre);
    PROXY_PARAMETER(GlobalState::get(), CapPressureParams, micParams);
    PROXY_PARAMETER(GlobalState::get(), CapPressureParams, mdcParams);

    PROPERTY(Scalar, Snrei, setSnrei);
    MUTABLE_PTR_PROPERTY(ScanningCurve, mdc, setMdc);
    MUTABLE_PTR_PROPERTY(ScanningCurve, pisc, setPisc);
    MUTABLE_PTR_PROPERTY(ScanningCurve, csc, setCsc);
};

///////////// Play-Type State
template <class GlobalState>
class MyPTState
{
public:
    typedef typename GlobalState::Scalar            Scalar;
    typedef typename GlobalState::CapPressureParams CapPressureParams;

    MyPTState()
        {
            Api::require<Api::PlayTypeState>(*this);

            _SweRef = 0;
        }

    PROXY_PARAMETER(GlobalState::get(), Scalar, Swr);
    PROXY_PARAMETER(GlobalState::get(), Scalar, Snr);
    PROXY_PARAMETER(GlobalState::get(), Scalar, Snre);
    PROXY_PARAMETER(GlobalState::get(), CapPressureParams, micParams);
    PROXY_PARAMETER(GlobalState::get(), CapPressureParams, mdcParams);

    PROPERTY(Scalar, SweRef,  setSweRef);
    PROPERTY(bool,   isImbib, setImbib);

    CONST_PARAMETER(Scalar, deltaSwe, 0.01);
};


template <class HystModel, class State>
void printCurve(const std::string &title,
                State &state,
                typename HystModel::Scalar a0,
                typename HystModel::Scalar a1,
                int n)
{
    typedef typename HystModel::Scalar Scalar;

    Scalar inc = (a1 - a0)/n;
    int i = 0;
    for (; i <= n; ++i) {
        Scalar absSat = a0 + inc*i;
        if (HystModel::isValidSw(state, absSat))
            break;
    }
    assert(i <= n);
    HystModel::updateState(state, a0 + inc*i);
    cerr << "Snrei: " << state.Snrei() << "\n";
    

    cout << "snip: " << title << "\n";
    for (; i <= n; ++i) {
        Scalar absSat = a0 + inc*i;

        if (!HystModel::isValidSw(state, absSat))
            continue;

//        cout << absSat << "  " << HystModel::krn(state, absSat) << endl;
//        cout << absSat << "  " << HystModel::krw(state, absSat) << endl;
        cout << absSat << "  " << HystModel::pC(state, absSat) << endl;
//        cout << absSat << "  " << HystModel::Sw(state, HystModel::pC(state, absSat)) << endl;
//        cout << absSat << "  " << HystModel::dpC_dSw(state, absSat) << endl;
//        cout << absSat << "  " << HystModel::dSw_dpC(state, HystModel::pC(state, absSat)) << endl;
//        cout << absSat << "  " << HystModel::dSw_dpC(state, HystModel::pC(state, absSat)) - 1/HystModel::dpC_dSw(state, absSat) << endl;

//        cout << HystModel::SwToSwapp(state, absSat) << "  " << HystModel::krn(state, absSat) << endl;
    }
    cout << "snap\n";
}

template <class HystModel, class State>
void printCycle(State &state, int n)
{
    printCurve<HystModel, State>("mic", state, 0, 1, n);
    HystModel::reset(state);
    printCurve<HystModel, State>("mdc", state, 1, 0, n);

#if 0
    typedef typename HystModel::Scalar Scalar;
    for (Scalar bla = 1.0; bla > 0.0; bla -= 0.01) {
        printCurve<HystModel, State>((boost::format("curve-%02.0f")%(bla*100)).str(), state, bla, 1.0, 100);
    }
#else
    printCurve<HystModel, State>("imbib-50->100", state, 0.5, 1.0, n);
    printCurve<HystModel, State>("drain-80->0", state,   0.80, 0.0, n);
    printCurve<HystModel, State>("imbib-60->100", state, 0.6, 1.0, n);
    printCurve<HystModel, State>("drain-75->0", state,   0.75, 0.0, n);
    printCurve<HystModel, State>("imbib-30->100", state, 0.30, 1.0, n);
#endif
}

void printVanGenuchtenMain()
{
    ////////////////////
    // Some type definitions
    ////////////////////
     typedef RegularizedVanGenuchtenState<Scalar> VanGenuchtenState;
     typedef RegularizedVanGenuchten<VanGenuchtenState>  VanGenuchten;
//    typedef VanGenuchtenState<Scalar> VanGenuchtenState;
//    typedef VanGenuchten<VanGenuchtenState>  VanGenuchten;

    typedef MyGlobalState<VanGenuchtenState> PLVGGlobalState;
    typedef MyPLState<PLVGGlobalState>         PLVGState;
    typedef ParkerLenhard<PLVGState, VanGenuchten> ParkerLenhardVG;

    ////////////////////
    // State instantiation
    ////////////////////
    // The two states for the MIC and the MDC required for hystersis
    VanGenuchtenState vgDrain(0.00045, 6.3);
    VanGenuchtenState vgImbib(0.00045*2, 6.5);
    vgImbib.setVgMaxPC(vgDrain.vgMaxPC());

    // For the global state we assume, that the parameters for the MIC
    // and the MDC are constant over the whole domain, as well as
    // the residual saturations of the phases.
    PLVGGlobalState  plvgGlobalState(Swr, Snr, vgImbib, vgDrain);

    // The local state is responsible for all parameters which vary
    // spatially (the current hysteresis state, ...). Here we only
    // look at one individual point, so this is not cruicial, but if
    // large domains (think > 10 million vertices) are used RAM usage
    // and cache coherence should _drastically_ improove.
    PLVGState        plvgState;

    ////////////////////
    // Start the actual "calculation"
    ////////////////////
    fprintf(stderr,
            "sizeof(Scalar): %d, sizeof(plvgGlobalState): %d, sizeof(plvgState): %d\n",
            (int) sizeof(Scalar),
            (int) sizeof(plvgGlobalState),
            (int) sizeof(plvgState));
    printCycle<ParkerLenhardVG, PLVGState>(plvgState, numSamples);
}
/*
void printBrooksCoreyMain()
{
    ////////////////////
    // Some type definitions
    ////////////////////
    typedef BrooksCoreyState<Scalar>        BrooksCoreyState;
    typedef BrooksCorey<BrooksCoreyState>   BrooksCorey;

    typedef MyGlobalState<BrooksCoreyState>     PLBCGlobalState;
    typedef MyPLState<PLBCGlobalState>            PLBCState;
    typedef ParkerLenhard<PLBCState, BrooksCorey> ParkerLenhardBC;

    ////////////////////
    // State instantiation
    ////////////////////
    // The two states for the MIC and the MDC required for hystersis
    BrooksCoreyState bcImbib(850., 4.1);
    BrooksCoreyState bcDrain(1550., 4.2);

    // For the global state we assume, that the parameters for the MIC
    // and the MDC are constant over the whole domain, as well as
    // the residual saturations of the phases.
    PLBCGlobalState  plbcGlobalState(Swr, Snr, bcImbib, bcDrain);

    // The local state is responsible for all parameters which vary
    // spatially (the current hysteresis state, ...). Here we only
    // look at one individual point, so this is not cruicial, but if
    // large domains (think > 10 million vertices) are used RAM usage
    // and cache coherence should _drastically_ improove.
    PLBCState        plbcState;

    ////////////////////
    // Start the actual "calculation"
    ////////////////////
    fprintf(stderr,
            "sizeof(Scalar): %d, sizeof(plbcGlobalState): %d, sizeof(plbcState): %d\n",
            (int) sizeof(Scalar),
            (int) sizeof(plbcGlobalState),
            (int) sizeof(plbcState));
    printCycle<ParkerLenhardBC, PLBCState>(plbcState, numSamples);
}

void printPlayBrooksCoreyMain()
{
    ////////////////////
    // Some type definitions
    ////////////////////
    typedef BrooksCoreyState<Scalar>        BrooksCoreyState;
    typedef BrooksCorey<BrooksCoreyState>   BrooksCorey;

    typedef MyGlobalState<BrooksCoreyState>       PTBCGlobalState;
    typedef MyPTState<PTBCGlobalState>            PTBCState;
    typedef PlayType<PTBCState, BrooksCorey>      PlayTypeBC;

    ////////////////////
    // State instantiation
    ////////////////////
    // The two states for the MIC and the MDC required for hystersis
    BrooksCoreyState bcImbib(850., 4.1);
    BrooksCoreyState bcDrain(1550., 4.2);

    // For the global state we assume, that the parameters for the MIC
    // and the MDC are constant over the whole domain, as well as
    // the residual saturations of the phases.
    PTBCGlobalState  ptbcGlobalState(Swr, Snr, bcImbib, bcDrain);

    // The local state is responsible for all parameters which vary
    // spatially (the current hysteresis state, ...). Here we only
    // look at one individual point, so this is not cruicial, but if
    // large domains (think > 10 million vertices) are used RAM usage
    // and cache coherence should _drastically_ improove.
    PTBCState        ptbcState;

    ////////////////////
    // Start the actual "calculation"
    ////////////////////
    fprintf(stderr,
            "sizeof(Scalar): %d, sizeof(plbcGlobalState): %d, sizeof(plbcState): %d\n",
            (int) sizeof(Scalar),
            (int) sizeof(ptbcGlobalState),
            (int) sizeof(ptbcState));
    printCycle<PlayTypeBC, PTBCState>(ptbcState, numSamples);
}
*/

int main()
{
    try
    {
        printVanGenuchtenMain();
//        printBrooksCoreyMain();
//        printPlayBrooksCoreyMain();
    }
    catch(std::exception &e) {
        cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}
