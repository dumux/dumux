// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EOS
 * \brief Implements the Peng-Robinson equation of state for liquids
 *        and gases.
 *
 * See:
 *
 * D.-Y. Peng, D.B. Robinson (1976, pp. 59–64) \cite peng1976 <BR>
 * https://doi.org/10.1021/i160057a011
 *
 * R. Reid, et al. "The properties of gases and liquids" (1987, pp. 42-44, 82)
 * \cite reid1987
 */
#ifndef DUMUX_PENG_ROBINSON_HH
#define DUMUX_PENG_ROBINSON_HH

#include <dune/common/math.hh>
#include <dumux/material/fluidstates/temperatureoverlay.hh>
#include <dumux/material/idealgas.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/tabulated2dfunction.hh>

#include <iostream>

namespace Dumux {

/*!
 * \ingroup EOS
 * \brief Implements the Peng-Robinson equation of state for liquids
 *        and gases.
 *
 * See:
 *
 * D.-Y. Peng, D.B. Robinson (1976, pp. 59–64) \cite peng1976 <BR>
 *
 * R. Reid, et al. (1987, pp. 42-44, 82) \cite reid1987
 */
template <class Scalar>
class PengRobinson
{
    PengRobinson()
    { }

public:
    static void  init(Scalar aMin, Scalar aMax, int na,
                      Scalar bMin, Scalar bMax, int nb)
    {
        // resize the tabulation for the critical points
        criticalTemperature_.resize(aMin, aMax, na, bMin, bMax, nb);
        criticalPressure_.resize(aMin, aMax, na, bMin, bMax, nb);
        criticalMolarVolume_.resize(aMin, aMax, na, bMin, bMax, nb);

        Scalar VmCrit = 18e-6;
        Scalar pCrit = 1e7;
        Scalar TCrit = 273;
        for (int i = 0; i < na; ++i) {
            Scalar a = criticalTemperature_.iToX(i);
            for (int j = 0; j < nb; ++j) {
                Scalar b = criticalTemperature_.jToY(j);

                findCriticalPoint_(TCrit, pCrit, VmCrit,  a, b);

                criticalTemperature_.setSamplePoint(i, j, TCrit);
                criticalPressure_.setSamplePoint(i, j, pCrit);
                criticalMolarVolume_.setSamplePoint(i, j, VmCrit);
            }
        }
    };

    /*!
     * \brief Predicts the vapor pressure \f$\mathrm{[Pa]}\f$ for the temperature given in
     *        setTP().
     * \param T temperature in \f$\mathrm{[K]}\f$
     * \param params Parameters
     *
     * Initially, the vapor pressure is roughly estimated by using the
     * Ambrose-Walton method, then the Newton method is used to make
     * difference between the gas and liquid phase fugacity zero.
     */
    template <class Params>
    static Scalar computeVaporPressure(const Params &params, Scalar T)
    {
        using Component = typename Params::Component;
        if (T >= Component::criticalTemperature())
            return Component::criticalPressure();

        // initial guess of the vapor pressure
        Scalar Vm[3];
        const Scalar eps = Component::criticalPressure()*1e-10;

        // use the Ambrose-Walton method to get an initial guess of
        // the vapor pressure
        Scalar pVap = ambroseWalton_(params, T);

        // Newton-Raphson method
        for (int i = 0; i < 5; ++i) {
            // calculate the molar densities
            assert(molarVolumes(Vm, params, T, pVap) == 3);

            Scalar f = fugacityDifference_(params, T, pVap, Vm[0], Vm[2]);
            Scalar df_dp =
                fugacityDifference_(params, T, pVap  + eps, Vm[0], Vm[2])
                -
                fugacityDifference_(params, T, pVap - eps, Vm[0], Vm[2]);
            df_dp /= 2*eps;

            Scalar delta = f/df_dp;
            pVap = pVap - delta;
            using std::abs;
            if (abs(delta/pVap) < 1e-10)
                break;
        }

        return pVap;
    }

    /*!
     * \brief Computes molar volumes \f$\mathrm{[m^3 / mol]}\f$ where the Peng-Robinson EOS is
     *        true.
     * \param fs Thermodynamic state of the fluids
     * \param params Parameters
     * \param phaseIdx The phase index
     * \param isGasPhase Specifies the phase state
     * \param handleUnphysicalPhase Special handling of the case when the EOS has only one
              intersection with the pressure, and the intersection does not correspond to
              the given phase (the phase is thus considered unphysical). If it happens in
              the case of critical fluid, the critical molar volume is returned for the
              unphysical phase. If the fluid is not critical, a proper extremum of the
              EOS is returned for the unphysical phase. If the parameter is false and the
              EOS has only one intersection with the pressure, the molar volume is computed
              from that single intersection, not depending of the given phase (gas or
              fluid). If the EOS has three intersections with the pressure, this parameter
              is ignored.
     */
    template <class FluidState, class Params>
    static Scalar computeMolarVolume(const FluidState &fs,
                                     Params &params,
                                     int phaseIdx,
                                     bool isGasPhase,
                                     bool handleUnphysicalPhase = true)
    {
        Scalar Vm = 0;

        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);
        
        Scalar a = params.a(phaseIdx); // "attractive factor"
        Scalar b = params.b(phaseIdx); // "co-volume"

        Scalar RT = Constants<Scalar>::R*T;
        Scalar Astar = a*p/(RT*RT);
        Scalar Bstar = b*p/RT;

        Scalar a1 = 1.0;
        Scalar a2 = - (1 - Bstar);
        Scalar a3 = Astar - Bstar*(3*Bstar + 2);
        Scalar a4 = Bstar*(- Astar + Bstar*(1 + Bstar));

        Scalar Z[3] = {}; // zero-initializer
        int numSol = invertCubicPolynomial(Z, a1, a2, a3, a4);

        // ignore the first two results if the smallest
        // compressibility factor is <= 0.0. (this means that if we
        // would get negative molar volumes for the liquid phase, we
        // consider the liquid phase non-existent.)
        // Note that invertCubicPolynomial sorts the roots in ascending order
        if (numSol > 1 && Z[0] <= 0.0) {
            Z[0] = Z[numSol - 1];
            numSol = 1;
        }

        if (numSol > 0 && Z[0] <= 0)
            DUNE_THROW(NumericalProblem, "No positive compressibility factors found");
        
        if (numSol == 3) {
            // the EOS has three intersections with the pressure,
            // i.e. the molar volume of gas is the largest one and the
            // molar volume of liquid is the smallest one
            if (isGasPhase)
                Vm = Z[2]*RT/p;
            else
                Vm = Z[0]*RT/p;
        }
        else if (numSol == 1) {
            // the EOS has only one intersection with the pressure
            
            Vm = Z[0]*RT/p;

            // Handle single root case specially unless told otherwise
            if (handleUnphysicalPhase)
                Vm = handleSingleRoot_(Vm, fs, params, phaseIdx, isGasPhase);
        }
        else
            DUNE_THROW(NumericalProblem, "Unexpected number of roots in the equation for compressibility factor: " << numSol);

        using std::isfinite;
        
        if (!isfinite(Vm) || Vm <= 0)
            DUNE_THROW(NumericalProblem, "Bad molar volume");

        return Vm;
    }

    /*!
     * \brief Returns the critical temperature for a given mix
     *
     * \param params EOS (equation of state) parameters of a single-component fluid \
     *        (usually PengRobinsonParms) or a mixture (usually PengRobinsonMixtureParams)
     */
    template <class Params>
    static Scalar criticalTemperature(const Params &params)
    {
        // Here, a() is the attractive force parameter and b()
        // is the repulsive force parameter of the EOS
        return criticalTemperature_(params.a(), params.b());
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ for a given pressure
     *        and molar volume.
     *
     * This is the same value as computeFugacity() because the mole
     * fraction of a component in a pure fluid is obviously always
     * 100%.
     *
     * \param params Parameters
     */
    template <class Params>
    static Scalar computeFugacityCoeffient(const Params &params)
    {
        Scalar T = params.temperature();
        Scalar p = params.pressure();
        Scalar Vm = params.molarVolume();

        Scalar RT = Constants<Scalar>::R*T;
        Scalar Z = p*Vm/RT;
        Scalar Bstar = p*params.b() / RT;
        using std::sqrt;
        Scalar tmp =
            (Vm + params.b()*(1 + sqrt(2))) /
            (Vm + params.b()*(1 - sqrt(2)));
        Scalar expo = - params.a()/(RT * 2 * params.b() * sqrt(2));
        using std::exp;
        using std::pow;
        Scalar fugCoeff =
            exp(Z - 1) / (Z - Bstar) *
            pow(tmp, expo);

        return fugCoeff;
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ for a given pressure
     *        and molar volume.
     *
     * This is the fugacity coefficient times the pressure. The mole
     * fraction of a component in a pure fluid is obviously always
     * 100%, so it is not required.
     *
     * \param params Parameters
     */
    template <class Params>
    static Scalar computeFugacity(const Params &params)
    { return params.pressure()*computeFugacityCoeff(params); }

protected:

    /*
     * Handle single-root case in the equation of state. The problem here is
     * that only one phase exists in this case, while we should also figure out
     * something to return for the other phase. For reference, see Andreas
     * Lauser's thesis: http://dx.doi.org/10.18419/opus-516 (page 32).
     */
    template <class FluidState, class Params>
    static Scalar handleSingleRoot_(Scalar Vm,
                                  const FluidState &fs,
                                  const Params &params,
                                  int phaseIdx,
                                  bool isGasPhase)
    {
        Scalar T = fs.temperature(phaseIdx);
        
        // for the other phase, we take the extremum of the EOS
        // with the largest distance from the intersection.
        if (T > criticalTemperature_(params.a(phaseIdx), params.b(phaseIdx))) {
            // if the EOS does not exhibit any extrema, the fluid
            // is critical...
            return handleCriticalFluid_(Vm, fs, params, phaseIdx, isGasPhase);
        }
        else {
            // find the extrema (if they are present)
            Scalar Vmin, Vmax, pmin, pmax;
            using std::min;
            using std::max;
            if (findExtrema_(Vmin, Vmax, pmin, pmax,
                             params.a(phaseIdx), params.b(phaseIdx), T))
                return isGasPhase ? max(Vmax, Vm) : min(Vmin, Vm);
        }
        return Vm;
    }
    
    template <class FluidState, class Params>
    static Scalar handleCriticalFluid_(Scalar Vm,
                                     const FluidState &fs,
                                     const Params &params,
                                     int phaseIdx,
                                     bool isGasPhase)
    {
        Scalar Vcrit = criticalMolarVolume_(params.a(phaseIdx), params.b(phaseIdx));

        using std::max;
        using std::min;
        return isGasPhase ? max(Vm, Vcrit) : min(Vm, Vcrit);
    }

    static void findCriticalPoint_(Scalar &Tcrit,
                                   Scalar &pcrit,
                                   Scalar &Vcrit,
                                   Scalar a,
                                   Scalar b)
    {
        Scalar minVm;
        Scalar maxVm;

        Scalar minP(0);
        Scalar maxP(1e100);

        // first, we need to find an isotherm where the EOS exhibits
        // a maximum and a minimum
        Scalar Tmin = 250; // [K]
        for (int i = 0; i < 30; ++i) {
            bool hasExtrema = findExtrema_(minVm, maxVm, minP, maxP, a, b, Tmin);
            if (hasExtrema)
                break;
            Tmin /= 2;
        }

        Scalar T = Tmin;

        // Newton's method: Start at minimum temperature and minimize
        // the "gap" between the extrema of the EOS
        for (int i = 0; i < 25; ++i) {
            // calculate function we would like to minimize
            Scalar f = maxVm - minVm;

            // check if we're converged
            if (f < 1e-10) {
                Tcrit = T;
                pcrit = (minP + maxP)/2;
                Vcrit = (maxVm + minVm)/2;
                return;
            }

            // backward differences. Forward differences are not
            // robust, since the isotherm could be critical if an
            // epsilon was added to the temperature. (this is case
            // rarely happens, though)
            const Scalar eps = - 1e-8;
            [[maybe_unused]] bool hasExtrema = findExtrema_(minVm, maxVm, minP, maxP, a, b, T + eps);
            assert(hasExtrema);
            Scalar fStar = maxVm - minVm;

            // derivative of the difference between the maximum's
            // molar volume and the minimum's molar volume regarding
            // temperature
            Scalar fPrime = (fStar - f)/eps;

            // update value for the current iteration
            Scalar delta = f/fPrime;
            if (delta > 0)
                delta = -10;

            // line search (we have to make sure that both extrema
            // still exist after the update)
            for (int j = 0; ; ++j) {
                if (j >= 20) {
                    DUNE_THROW(NumericalProblem,
                               "Could not determine the critical point for a=" << a << ", b=" << b);
                }

                if (findExtrema_(minVm, maxVm, minP, maxP, a, b, T - delta)) {
                    // if the isotherm for T - delta exhibits two
                    // extrema the update is finished
                    T -= delta;
                    break;
                }
                else
                    delta /= 2;
            }
        }
    }

    // find the two molar volumes where the EOS exhibits extrema and
    // which are larger than the covolume of the phase
    static bool findExtrema_(Scalar &Vmin,
                             Scalar &Vmax,
                             Scalar &pMin,
                             Scalar &pMax,
                             Scalar a,
                             Scalar b,
                             Scalar T)
    {
        Scalar u = 2;
        Scalar w = -1;

        Scalar RT = Constants<Scalar>::R*T;

        // calculate coefficients of the 4th order polynominal in
        // monomial basis
        Scalar a1 = RT;
        Scalar a2 = 2*RT*u*b - 2*a;
        Scalar a3 = 2*RT*w*b*b + RT*u*u*b*b  + 4*a*b - u*a*b;
        Scalar a4 = 2*RT*u*w*b*b*b + 2*u*a*b*b - 2*a*b*b;
        Scalar a5 = RT*w*w*b*b*b*b - u*a*b*b*b;

        // Newton method to find first root

        // if the values which we got on Vmin and Vmax are useful, we
        // will reuse them as initial value, else we will start 10%
        // above the covolume
        Scalar V = b*1.1;
        Scalar delta = 1.0;
        using std::abs;
        for (int i = 0; abs(delta) > 1e-9; ++i) {
            Scalar f = a5 + V*(a4 + V*(a3 + V*(a2 + V*a1)));
            Scalar fPrime = a4 + V*(2*a3 + V*(3*a2 + V*4*a1));

            if (abs(fPrime) < 1e-20) {
                // give up if the derivative is zero
                Vmin = 0;
                Vmax = 0;
                return false;
            }


            delta = f/fPrime;
            V -= delta;

            if (i > 20) {
                // give up after 20 iterations...
                Vmin = 0;
                Vmax = 0;
                return false;
            }
        }

        // polynomial division
        Scalar b1 = a1;
        Scalar b2 = a2 + V*b1;
        Scalar b3 = a3 + V*b2;
        Scalar b4 = a4 + V*b3;

        // invert resulting cubic polynomial analytically
        Scalar allV[4];
        allV[0] = V;
        int numSol = 1 + invertCubicPolynomial(allV + 1, b1, b2, b3, b4);

        // sort all roots of the derivative
        std::sort(allV + 0, allV + numSol);

        // check whether the result is physical
        if (allV[numSol - 2] < b) {
            // the second largest extremum is smaller than the phase's
            // covolume which is physically impossible
            Vmin = Vmax = 0;
            return false;
        }

        // it seems that everything is okay...
        Vmin = allV[numSol - 2];
        Vmax = allV[numSol - 1];
        return true;
    }

    /*!
     * \brief The Ambrose-Walton method to estimate the vapor
     *        pressure.
     *
     * \return Vapor pressure estimate in bar
     *
     * See:
     *
     * D. Ambrose, J. Walton (1989, pp. 1395-1403) \cite ambrose1989
     */
    template <class Params>
    static Scalar ambroseWalton_(const Params &params, Scalar T)
    {
        using Component = typename Params::Component;

        Scalar Tr = T / Component::criticalTemperature();
        Scalar tau = 1 - Tr;
        Scalar omega = Component::acentricFactor();
        using std::sqrt;
        using Dune::power;
        Scalar f0 = (tau*(-5.97616 + sqrt(tau)*(1.29874 - tau*0.60394)) - 1.06841*power(tau, 5))/Tr;
        Scalar f1 = (tau*(-5.03365 + sqrt(tau)*(1.11505 - tau*5.41217)) - 7.46628*power(tau, 5))/Tr;
        Scalar f2 = (tau*(-0.64771 + sqrt(tau)*(2.41539 - tau*4.26979)) + 3.25259*power(tau, 5))/Tr;
        using std::exp;
        return Component::criticalPressure()*exp(f0 + omega * (f1 + omega*f2));
    }

    /*!
     * \brief Returns the difference between the liquid and the gas phase
     *        fugacities in [bar]
     *
     * \param params Parameters
     * \param T Temperature [K]
     * \param p Pressure [bar]
     * \param VmLiquid Molar volume of the liquid phase [cm^3/mol]
     * \param VmGas Molar volume of the gas phase [cm^3/mol]
     */
    template <class Params>
    static Scalar fugacityDifference_(const Params &params,
                                      Scalar T,
                                      Scalar p,
                                      Scalar VmLiquid,
                                      Scalar VmGas)
    { return fugacity(params, T, p, VmLiquid) - fugacity(params, T, p, VmGas); }

    static Tabulated2DFunction<Scalar> criticalTemperature_;
    static Tabulated2DFunction<Scalar> criticalPressure_;
    static Tabulated2DFunction<Scalar> criticalMolarVolume_;
};

template <class Scalar>
Tabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalTemperature_;

template <class Scalar>
Tabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalPressure_;

template <class Scalar>
Tabulated2DFunction<Scalar> PengRobinson<Scalar>::criticalMolarVolume_;

} // end namespace

#endif
