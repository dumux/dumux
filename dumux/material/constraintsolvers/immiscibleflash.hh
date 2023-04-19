// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ConstraintSolvers
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 */
#ifndef DUMUX_IMMISCIBLE_FLASH_HH
#define DUMUX_IMMISCIBLE_FLASH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>

namespace Dumux {

/*!
 * \ingroup ConstraintSolvers
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 *
 * In a N-phase, N-component context, we have the following
 * unknowns if assuming immiscibility:
 *
 * - N pressures
 * - N saturations
 *
 * This sums up to 2*N unknowns. On the equations side of things, we
 * have:
 *
 * - N total component molarities
 * - 1 The sum of all saturations is 1
 * - N-1 Relations from capillary pressure
 *
 * this also sums up to 2*N. We include the capillary pressures and
 * the sum of the saturations explicitly. This means that we only
 * solve for the first pressure and N-1 saturations.
 *
 * If a fluid phase is incompressible, the pressure cannot determined
 * by this, though. In this case the original pressure is kept, and
 * the saturation of the phase is calculated by dividing the global
 * molarity of the component by the phase density.
 */
template <class Scalar, class FluidSystem>
class ImmiscibleFlash
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static_assert(numPhases == numComponents,
                  "Immiscibility assumes that the number of phases"
                  " is equal to the number of components");

    using ParameterCache = typename FluidSystem::ParameterCache;

    static constexpr int numEq = numPhases;

    using Matrix = Dune::FieldMatrix<Scalar, numEq, numEq>;
    using Vector = Dune::FieldVector<Scalar, numEq>;

public:
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    /*!
     * \brief Construct a new flash
     * \param wettingPhaseIdx the phase index of the wetting phase
     */
    explicit ImmiscibleFlash(int wettingPhaseIdx = 0)
    : wettingPhaseIdx_(wettingPhaseIdx)
    {}

    /*!
     * \brief Guess initial values for all quantities.
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param globalMolarities
     */
    template <class FluidState>
    void guessInitial(FluidState &fluidState,
                      ParameterCache &paramCache,
                      const ComponentVector &globalMolarities)
    {
        // the sum of all molarities
        Scalar sumMoles = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoles += globalMolarities[compIdx];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 2e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param globalMolarities
     * \param material The material law object
     *
     * The phase's fugacities must already be set.
     */
    template <class MaterialLaw, class FluidState>
    void solve(FluidState& fluidState,
               ParameterCache& paramCache,
               const MaterialLaw& material,
               const ComponentVector& globalMolarities)
    {
        /////////////////////////
        // Check if all fluid phases are incompressible
        /////////////////////////
        bool allIncompressible = true;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::isCompressible(phaseIdx)) {
                allIncompressible = false;
                break;
            }
        }

        if (allIncompressible) {
            paramCache.updateAll(fluidState);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);
                fluidState.setMolarDensity(phaseIdx, rhoMolar);

                Scalar saturation =
                    globalMolarities[/*compIdx=*/phaseIdx]
                    / fluidState.molarDensity(phaseIdx);
                fluidState.setSaturation(phaseIdx, saturation);
            }
        }

//         /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Matrix J;
        // solution, i.e. phase composition
        Vector deltaX;
        // right hand side
        Vector b;

        completeFluidState_(fluidState, paramCache, material);

        const int nMax = 50; // <- maximum number of newton iterations
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_(J, b, fluidState, paramCache, material, globalMolarities);

            // Solve J*x = b
            deltaX = 0;

            try { J.solve(deltaX, b); }
            catch (Dune::FMatrixError& e)
            {
                /*
                std::cout << "error: " << e << "\n";
                std::cout << "b: " << b << "\n";
                */

                throw NumericalProblem(e.what());
            }

            /*
            printFluidState_(fluidState);
            std::cout << "J:\n";
            for (int i = 0; i < numEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    std::ostringstream os;
                    os << J[i][j];

                    std::string s(os.str());
                    do {
                        s += " ";
                    } while (s.size() < 20);
                    std::cout << s;
                }
                std::cout << "\n";
            }
            std::cout << "residual: " << b << "\n";
            std::cout << "deltaX: " << deltaX << "\n";
            std::cout << "---------------\n";
            */

            // update the fluid quantities.
            Scalar relError = update_(fluidState, paramCache, material, deltaX);

            if (relError < 1e-9)
                return;
        }

        /*
        printFluidState_(fluidState);
        std::cout << "globalMolarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */

        DUNE_THROW(NumericalProblem,
                   "Flash calculation failed."
                   " {c_alpha^kappa} = {" << globalMolarities << "}, T = "
                   << fluidState.temperature(/*phaseIdx=*/0));
    }


protected:
    template <class FluidState>
    void printFluidState_(const FluidState &fs)
    {
        std::cout << "saturations: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.density(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "global component molarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar sum = 0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fs.saturation(phaseIdx)*fs.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw, class FluidState>
    void linearize_(Matrix &J,
                    Vector &b,
                    FluidState &fluidState,
                    ParameterCache &paramCache,
                    const MaterialLaw& material,
                    const ComponentVector &globalMolarities)
    {
        FluidState origFluidState(fluidState);
        ParameterCache origParamCache(paramCache);

        Vector tmp;

        // reset jacobian
        J = 0;

        calculateDefect_(b, fluidState, fluidState, globalMolarities);

        // assemble jacobian matrix
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar x_i = getQuantity_(fluidState, pvIdx);
            const Scalar eps = 1e-10/quantityWeight_(fluidState, pvIdx);
            setQuantity_<MaterialLaw>(fluidState, paramCache, material, pvIdx, x_i + eps);

            // compute derivative of the defect
            calculateDefect_(tmp, origFluidState, fluidState, globalMolarities);
            tmp -= b;
            tmp /= eps;
            // store derivative in jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                J[eqIdx][pvIdx] = tmp[eqIdx];

            // fluid state and parameter cache to their original values
            fluidState = origFluidState;
            paramCache = origParamCache;

            // end forward differences
            ////////
        }
    }

    template <class FluidState>
    void calculateDefect_(Vector &b,
                          const FluidState &fluidStateEval,
                          const FluidState &fluidState,
                          const ComponentVector &globalMolarities)
    {
        // global molarities are given
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            b[phaseIdx] =
                fluidState.saturation(phaseIdx)
                * fluidState.molarity(phaseIdx, /*compIdx=*/phaseIdx);

            b[phaseIdx] -= globalMolarities[/*compIdx=*/phaseIdx];
        }
    }

    template <class MaterialLaw, class FluidState>
    Scalar update_(FluidState &fluidState,
                   ParameterCache &paramCache,
                   const MaterialLaw& material,
                   const Vector &deltaX)
    {
        Scalar relError = 0;
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            Scalar tmp = getQuantity_(fluidState, pvIdx);
            Scalar delta = deltaX[pvIdx];
            using std::max;
            using std::abs;
            using std::min;
            relError = max(relError, abs(delta)*quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 20% change in saturation per
                // iteration
                delta = min(0.2, max(-0.2, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 30% change in pressure per
                // iteration
                delta = min(0.30*fluidState.pressure(0), max(-0.30*fluidState.pressure(0), delta));
            }

            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }

        completeFluidState_<MaterialLaw>(fluidState, paramCache, material);

        return relError;
    }

    template <class MaterialLaw, class FluidState>
    void completeFluidState_(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const MaterialLaw& material)
    {
        // calculate the saturation of the last phase as a function of
        // the other saturations
        Scalar sumSat = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += fluidState.saturation(phaseIdx);

        if (sumSat > 1.0) {
            // make sure that the last saturation does not become
            // negative
            for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            {
                Scalar S = fluidState.saturation(phaseIdx);
                fluidState.setSaturation(phaseIdx, S/sumSat);
            }
            sumSat = 1;
        }

        // set the last saturation
        fluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        const auto pc = material.capillaryPressures(fluidState, wettingPhaseIdx_);
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx,
                                   fluidState.pressure(0)
                                   + (pc[phaseIdx] - pc[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParameterCache::Temperature|ParameterCache::Composition);

        // update all densities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
            fluidState.setMolarDensity(phaseIdx, rhoMolar);
        }
    }

    bool isPressureIdx_(int pvIdx)
    { return pvIdx == 0; }

    bool isSaturationIdx_(int pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    Scalar getQuantity_(const FluidState &fs, int pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            return fs.pressure(phaseIdx);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            int phaseIdx = pvIdx - 1;
            return fs.saturation(phaseIdx);
        }
    }

    // set a quantity in the fluid state
    template <class MaterialLaw, class FluidState>
    void setQuantity_(FluidState &fs,
                      ParameterCache &paramCache,
                      const MaterialLaw& material,
                      int pvIdx,
                      Scalar value)
    {
        assert(0 <= pvIdx && pvIdx < numEq);

        if (pvIdx < 1) {
            // -> first pressure
            Scalar delta = value - fs.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, fs.pressure(phaseIdx) + delta);
            paramCache.updateAllPressures(fs);

            // update all densities
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                Scalar rhoMolar = FluidSystem::molarDensity(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
                fs.setMolarDensity(phaseIdx, rhoMolar);
            }
        }
        else {
            assert(pvIdx < numPhases);

            // -> first M-1 saturations
            Scalar delta = value - fs.saturation(/*phaseIdx=*/pvIdx - 1);
            fs.setSaturation(/*phaseIdx=*/pvIdx - 1, value);

            // set last saturation (-> minus the change of the
            // saturation of the other phases)
            fs.setSaturation(/*phaseIdx=*/numPhases - 1,
                             fs.saturation(numPhases - 1) - delta);

            // update all fluid pressures using the capillary pressure
            // law
            const auto pc = material.capillaryPressures(fs, wettingPhaseIdx_);
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx,
                               fs.pressure(0)
                               + (pc[phaseIdx] - pc[0]));
            paramCache.updateAllPressures(fs);

            // update all densities
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                Scalar rhoMolar = FluidSystem::molarDensity(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
                fs.setMolarDensity(phaseIdx, rhoMolar);

            }
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    void setQuantityRaw_(FluidState &fs, int pvIdx, Scalar value)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            fs.setPressure(phaseIdx, value);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            int phaseIdx = pvIdx - 1;

            // make sure that the first M-1 saturations does not get
            // negative
            using std::max;
            value = max(0.0, value);
            fs.setSaturation(phaseIdx, value);
        }
    }

    template <class FluidState>
    Scalar quantityWeight_(const FluidState &fs, int pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1.0/1e5;
        // first M - 1 saturations
        else {
            assert(pvIdx < numPhases);
            return 1.0;
        }
    }

private:
    int wettingPhaseIdx_ = 0; //!< the phase index of the wetting phase

};

} // end namespace Dumux

#endif
