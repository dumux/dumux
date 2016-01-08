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
 * \brief This local operator extends the 2cstokes2p2clocaloperator
 *        by non-isothermal conditions.
 */
#ifndef DUMUX_TWOCNISTOKES2P2CNILOCALOPERATOR_HH
#define DUMUX_TWOCNISTOKES2P2CNILOCALOPERATOR_HH

#include <dumux/multidomain/2cstokes2p2c/2cstokes2p2clocaloperator.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief The extension of the local operator for the coupling of a two-component Stokes model
 *        and a two-phase two-component Darcy model for non-isothermal conditions.
 *
 * This model implements the coupling between a free-flow model
 * and a porous-medium flow model under non-isothermal conditions.
 * Here the coupling conditions for the individual balance are presented:
 *
 * The total mass balance equation:
 * \f[
 *  \left[
 *    \left( \varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g} \right) \cdot \boldsymbol{n}
 *  \right]^\textrm{ff}
 *  = -\left[
 *      \left( \varrho_\textrm{g} \boldsymbol{v}_\textrm{g}
 *             + \varrho_\textrm{l} \boldsymbol{v}_\textrm{l} \right) \cdot \boldsymbol{n}
 *    \right]^\textrm{pm}
 * \f]
 * in which \f$n\f$ represents a vector normal to the interface pointing outside of
 * the specified subdomain.
 *
 * The momentum balance (tangential), which corresponds to the Beavers-Jospeh Saffman condition:
 * \f[
 *  \left[
 *   \left( {\boldsymbol{v}}_\textrm{g}
 *          + \frac{\sqrt{\left(\boldsymbol{K} \boldsymbol{t}_i \right) \cdot \boldsymbol{t}_i}}
 *                 {\alpha_\textrm{BJ} \mu_\textrm{g}} \boldsymbol{{\tau}}_\textrm{t} \boldsymbol{n}
 *          \right) \cdot \boldsymbol{t}_i
 *  \right]^\textrm{ff}
 *  = 0
 * \f]
 * with
 * \f$
 * \boldsymbol{{\tau}_\textrm{t}} = \left[ \mu_\textrm{g} + \mu_\textrm{g,t} \right]
 *                                  \nabla \left( \boldsymbol{v}_\textrm{g}
 *                                                + \boldsymbol{v}_\textrm{g}^\intercal \right)
 * \f$
 * in which the eddy viscosity \f$ \mu_\textrm{g,t} = 0 \f$ for the Stokes equation.
 *
 * The momentum balance (normal):
 * \f[
 *  \left[
 *    \left(
 *      \left\lbrace
 *        \varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g} {\boldsymbol{v}}_\textrm{g}^\intercal
 *        - \boldsymbol{{\tau}}_\textrm{t}
 *        + {p}_\textrm{g} \boldsymbol{I}
 *      \right\rbrace \boldsymbol{n}
 *    \right) \cdot \boldsymbol{n}
 *  \right]^\textrm{ff}
 *  = p_\textrm{g}^\textrm{pm}
 * \f]
 *
 * The component mass balance equation (continuity of fluxes):
 * \f[
 *  \left[
 *    \left(
 *      \varrho_\textrm{g} {X}^\kappa_\textrm{g} {\boldsymbol{v}}_\textrm{g}
 *      - {\boldsymbol{j}}^\kappa_\textrm{g,ff,t,diff}
 *    \right) \cdot \boldsymbol{n}
 *  \right]^\textrm{ff}
 *  = -\left[
 *    \left(
 *      \varrho_\textrm{g} X^\kappa_\textrm{g} \boldsymbol{v}_\textrm{g}
 *      - \boldsymbol{j}^\kappa_\textrm{g,pm,diff}
 *      + \varrho_\textrm{l} \boldsymbol{v}_\textrm{l} X^\kappa_\textrm{l}
 *      - \boldsymbol{j}^\kappa_\textrm{l,pm,diff}
 *    \right) \cdot \boldsymbol{n}
 *  \right]^\textrm{pm}
 *  = 0
 * \f]
 * in which the diffusive fluxes \f$ j_\textrm{diff} \f$ are the diffusive fluxes as
 * they are implemented in the individual subdomain models.
 *
 * The component mass balance equation (continuity of mass/ mole fractions):
 * \f[
 *  \left[ {X}^{\kappa}_\textrm{g} \right]^\textrm{ff}
 *  = \left[ X^{\kappa}_\textrm{g} \right]^\textrm{pm}
 * \f]
 *
 * The energy balance equation (continuity of fluxes):
 * \f[
 *  \left[
 *    \left(
 *      \varrho_\textrm{g} {h}_\textrm{g} {\boldsymbol{v}}_\textrm{g}
 *      - {h}^\textrm{a}_\textrm{g} {\boldsymbol{j}}^\textrm{a}_\textrm{g,ff,t,diff}
 *      - {h}^\textrm{w}_\textrm{g} {\boldsymbol{j}}^\textrm{w}_\textrm{g,ff,t,diff}
 *      - \left( \lambda_\textrm{g} + \lambda_\textrm{g,t} \right) \nabla {T}
 *    \right) \cdot \boldsymbol{n}
 *  \right]^\textrm{ff}
 *  = -\left[
 *       \left(
 *         \varrho_\textrm{g} h_\textrm{g} \boldsymbol{v}_\textrm{g}
 *         + \varrho_\textrm{l} h_\textrm{l} \boldsymbol{v}_\textrm{l}
 *         - \lambda_\textrm{pm} \nabla T
 *      \right) \cdot \boldsymbol{n}
 *    \right]^\textrm{pm}
 * \f]
 *
 * The energy balance equation (continuity of temperature):
 * \f[
 *  \left[ {T} \right]^\textrm{ff}
 *  = \left[ T \right]^\textrm{pm}
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class TwoCNIStokesTwoPTwoCNILocalOperator :
        public TwoCStokesTwoPTwoCLocalOperator<TypeTag>
{
public:
    typedef TwoCStokesTwoPTwoCLocalOperator<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) GlobalProblem;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) Stokes2cniTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TwoPTwoCNITypeTag;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, FluxVariables) BoundaryVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, FluxVariables) BoundaryVariables2;

    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, GridView) Stokes2cniGridView;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, GridView) TwoPTwoCNIGridView;
    typedef typename Stokes2cniGridView::template Codim<0>::Entity SDElement1;
    typedef typename TwoPTwoCNIGridView::template Codim<0>::Entity SDElement2;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, Indices) Stokes2cniIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, Indices) TwoPTwoCNIIndices;

    enum {
        dimWorld = MDGrid::dimensionworld
    };

    // Stokes
    enum { numComponents1 = Stokes2cniIndices::numComponents };
    enum { // equation indices
        energyEqIdx1 = Stokes2cniIndices::energyEqIdx             //!< Index of the energy balance equation
    };
    enum { // component indices
        transportCompIdx1 = Stokes2cniIndices::transportCompIdx,  //!< Index of transported component
        phaseCompIdx1 = Stokes2cniIndices::phaseCompIdx           //!< Index of main component of the phase
    };

    // Darcy
    enum { numPhases2 = GET_PROP_VALUE(TwoPTwoCNITypeTag, NumPhases) };
    enum { // equation indices
        energyEqIdx2 = TwoPTwoCNIIndices::energyEqIdx      //!< Index of the energy balance equation
    };
    enum { // phase indices
        wPhaseIdx2 = TwoPTwoCNIIndices::wPhaseIdx,          //!< Index for the liquid phase
        nPhaseIdx2 = TwoPTwoCNIIndices::nPhaseIdx           //!< Index for the gas phase
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename MDGrid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    // multidomain flags
    static const bool doAlphaCoupling = true;
    static const bool doPatternCoupling = true;

    TwoCNIStokesTwoPTwoCNILocalOperator(GlobalProblem& globalProblem)
        : ParentType(globalProblem)
    { }

public:
    //! \copydoc Dumux::TwoCStokesTwoPTwoCLocalOperator::evalCoupling()
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    void evalCoupling(const LFSU1& lfsu1, const LFSU2& lfsu2,
                      const int vertInElem1, const int vertInElem2,
                      const SDElement1& sdElement1, const SDElement2& sdElement2,
                      const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                      const CParams &cParams,
                      RES1& couplingRes1, RES2& couplingRes2) const
    {
        // evaluate coupling of mass and momentum balances
        ParentType::evalCoupling(lfsu1, lfsu2,
                                 vertInElem1, vertInElem2,
                                 sdElement1, sdElement2,
                                 boundaryVars1, boundaryVars2,
                                 cParams,
                                 couplingRes1, couplingRes2);

        const GlobalPosition& globalPos1 = cParams.fvGeometry1.subContVol[vertInElem1].global;
        const GlobalPosition& globalPos2 = cParams.fvGeometry2.subContVol[vertInElem2].global;

        // ENERGY Balance
        // Neumann-like conditions
        if (cParams.boundaryTypes1.isCouplingNeumann(energyEqIdx1))
        {
            if (this->globalProblem().sdProblem2().isCornerPoint(globalPos2))
            {
                // convective energy flux (enthalpy is mass based, mass flux also needed for useMoles)
                Scalar convectiveFlux = 0.0;
                for (int phaseIdx=0; phaseIdx<numPhases2; ++phaseIdx)
                {
                    convectiveFlux -= boundaryVars2.volumeFlux(phaseIdx)
                                      * cParams.elemVolVarsCur2[vertInElem2].density(phaseIdx)
                                      * cParams.elemVolVarsCur2[vertInElem2].enthalpy(phaseIdx);
                }

                // conductive energy flux
                Scalar conductiveFlux = boundaryVars2.normalMatrixHeatFlux();

                couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                        -(convectiveFlux - conductiveFlux));
            }
            else
            {
                couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                        this->globalProblem().localResidual2().residual(vertInElem2)[energyEqIdx2]);
            }
        }

        // TODO: unify the behavior for cParams.boundaryTypes2.isCouplingNeumann()
        //       with the different one in the isothermal LOP
        if (cParams.boundaryTypes2.isCouplingNeumann(energyEqIdx2))
        {
            const GlobalPosition& bfNormal1 = boundaryVars1.face().normal;
            // only enter here, if a boundary layer model is used for the computation of the diffusive fluxes
            if (ParentType::blModel_)
            {
                // convective energy flux (enthalpy is mass based, mass flux also needed for useMoles)
                Scalar convectiveFlux = boundaryVars1.normalVelocity()
                                        * cParams.elemVolVarsCur1[vertInElem1].density()
                                        * cParams.elemVolVarsCur1[vertInElem1].enthalpy();

                // conductive energy flux
                Scalar conductiveFlux = bfNormal1.two_norm()
                                        * evalBoundaryLayerTemperatureGradient(cParams, vertInElem1)
                                        * (boundaryVars1.thermalConductivity()
                                           + boundaryVars1.thermalEddyConductivity());

                // enthalpy transported by diffusive fluxes
                Scalar sumDiffusiveFluxes = 0.0;
                Scalar sumDiffusiveEnergyFlux = 0.0;
                for (int compIdx=0; compIdx < numComponents1; compIdx++)
                {
                    if (compIdx != phaseCompIdx1)
                    {
                        Scalar diffusiveFlux = bfNormal1.two_norm()
                                               * ParentType::evalBoundaryLayerConcentrationGradient(cParams, vertInElem1)
                                               * (boundaryVars1.diffusionCoeff(compIdx)
                                                  + boundaryVars1.eddyDiffusivity())
                                               * boundaryVars1.molarDensity()
                                               * ParentType::evalMassTransferCoefficient(cParams, vertInElem1, vertInElem2);
                        sumDiffusiveFluxes += diffusiveFlux;
                        sumDiffusiveEnergyFlux += diffusiveFlux
                                                  * boundaryVars1.componentEnthalpy(compIdx)
                                                  * FluidSystem::molarMass(compIdx); // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s]
                    }
                }
                sumDiffusiveEnergyFlux -= sumDiffusiveFluxes
                                          * boundaryVars1.componentEnthalpy(phaseCompIdx1)
                                          * FluidSystem::molarMass(phaseCompIdx1);

                // TODO: use mass transfer coefficient here?
                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        -(convectiveFlux - sumDiffusiveEnergyFlux - conductiveFlux));
            }
            else if (this->globalProblem().sdProblem1().isCornerPoint(globalPos1))
            {
                // convective energy flux (enthalpy is mass based, mass flux also needed for useMoles)
                Scalar convectiveFlux = boundaryVars1.normalVelocity()
                                        * cParams.elemVolVarsCur1[vertInElem1].density()
                                        * cParams.elemVolVarsCur1[vertInElem1].enthalpy();

                // conductive energy flux
                Scalar conductiveFlux = bfNormal1
                                        * boundaryVars1.temperatureGrad()
                                        * (boundaryVars1.thermalConductivity()
                                           + boundaryVars1.thermalEddyConductivity());

                // enthalpy transported by diffusive fluxes
                Scalar sumDiffusiveFluxes = 0.0;
                Scalar sumDiffusiveEnergyFlux = 0.0;
                for (int compIdx=0; compIdx < numComponents1; compIdx++)
                {
                    if (compIdx != phaseCompIdx1)
                    {
                        Scalar diffusiveFlux = bfNormal1
                                               * boundaryVars1.moleFractionGrad(compIdx)
                                               * (boundaryVars1.diffusionCoeff(compIdx)
                                                  + boundaryVars1.eddyDiffusivity())
                                               * boundaryVars1.molarDensity();
                        sumDiffusiveFluxes += diffusiveFlux;
                        sumDiffusiveEnergyFlux += diffusiveFlux
                                                  * boundaryVars1.componentEnthalpy(compIdx)
                                                  * FluidSystem::molarMass(compIdx); // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s]
                    }
                }
                sumDiffusiveEnergyFlux -= sumDiffusiveFluxes
                                          * boundaryVars1.componentEnthalpy(phaseCompIdx1)
                                          * FluidSystem::molarMass(phaseCompIdx1);

                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        -(convectiveFlux - sumDiffusiveEnergyFlux - conductiveFlux));
            }
            else
            {
                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        this->globalProblem().localResidual1().residual(vertInElem1)[energyEqIdx1]);
            }
        }

        // Dirichlet-like conditions
        if (cParams.boundaryTypes1.isCouplingDirichlet(energyEqIdx1))
        {
            // set residualStokes[energyIdx1] = T in stokesncnicouplinglocalresidual.hh
            couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                    -cParams.elemVolVarsCur2[vertInElem2].temperature());
        }

        if (cParams.boundaryTypes2.isCouplingDirichlet(energyEqIdx2))
        {
            // set residualDarcy[energyEqIdx2] = T in 2p2cnicouplinglocalresidual.hh
            couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                    -cParams.elemVolVarsCur1[vertInElem1].temperature());
        }
    }

    //! \copydoc Dumux::TwoCStokesTwoPTwoCLocalOperator::evalCoupling12()
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    DUNE_DEPRECATED_MSG("evalCoupling12() is deprecated. Use evalCoupling() instead.")
    void evalCoupling12(const LFSU1& lfsu1, const LFSU2& lfsu2,
                        const int vertInElem1, const int vertInElem2,
                        const SDElement1& sdElement1, const SDElement2& sdElement2,
                        const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                        const CParams &cParams,
                        RES1& couplingRes1, RES2& couplingRes2) const
    {
        const GlobalPosition& globalPos1 = cParams.fvGeometry1.subContVol[vertInElem1].global;
        const GlobalPosition& bfNormal1 = boundaryVars1.face().normal;
        const Scalar normalMassFlux1 = boundaryVars1.normalVelocity() *
            cParams.elemVolVarsCur1[vertInElem1].density();
        GlobalProblem& globalProblem = this->globalProblem();

        // evaluate coupling of mass and momentum balances
        ParentType::evalCoupling12(lfsu1, lfsu2,
                                   vertInElem1, vertInElem2,
                                   sdElement1, sdElement2,
                                   boundaryVars1, boundaryVars2,
                                   cParams,
                                   couplingRes1, couplingRes2);

        if (cParams.boundaryTypes2.isCouplingNeumann(energyEqIdx2))
        {
            // only enter here, if a boundary layer model is used for the computation of the diffusive fluxes
            if (ParentType::blModel_)
            {
                // convective energy flux
                const Scalar convectiveFlux = normalMassFlux1
                                              * cParams.elemVolVarsCur1[vertInElem1].enthalpy();

                // enthalpy transported by diffusive fluxes
                // multiply the diffusive flux with the mass transfer coefficient
                static_assert(numComponents1 == 2,
                              "This coupling condition is only implemented for two components.");
                Scalar diffusiveEnergyFlux = 0.0;
                Scalar diffusiveFlux = bfNormal1.two_norm()
                                       * ParentType::evalBoundaryLayerConcentrationGradient(cParams, vertInElem1)
                                       * (boundaryVars1.diffusionCoeff(transportCompIdx1)
                                          + boundaryVars1.eddyDiffusivity())
                                       * boundaryVars1.molarDensity()
                                       * ParentType::evalMassTransferCoefficient(cParams, vertInElem1, vertInElem2);

                diffusiveEnergyFlux += diffusiveFlux * FluidSystem::molarMass(transportCompIdx1)
                                       * boundaryVars1.componentEnthalpy(transportCompIdx1);
                diffusiveEnergyFlux -= diffusiveFlux * FluidSystem::molarMass(phaseCompIdx1)
                                       * boundaryVars1.componentEnthalpy(phaseCompIdx1);

                // conductive transported energy
                const Scalar conductiveFlux = bfNormal1.two_norm()
                                              * evalBoundaryLayerTemperatureGradient(cParams, vertInElem1)
                                              * (boundaryVars1.thermalConductivity()
                                                 + boundaryVars1.thermalEddyConductivity());

                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        -(convectiveFlux - diffusiveEnergyFlux - conductiveFlux));
            }
            else if (globalProblem.sdProblem1().isCornerPoint(globalPos1))
            {
                const Scalar convectiveFlux =
                    normalMassFlux1 *
                    cParams.elemVolVarsCur1[vertInElem1].enthalpy();
                const Scalar conductiveFlux =
                    bfNormal1 *
                    boundaryVars1.temperatureGrad() *
                    (boundaryVars1.thermalConductivity() + boundaryVars1.thermalEddyConductivity());
                Scalar sumDiffusiveFluxes = 0.0;
                Scalar sumDiffusiveEnergyFlux = 0.0;
                for (int compIdx=0; compIdx < numComponents1; compIdx++)
                {
                    if (compIdx != phaseCompIdx1)
                    {
                        Scalar diffusiveFlux = boundaryVars1.moleFractionGrad(compIdx)
                                               * boundaryVars1.face().normal
                                               *(boundaryVars1.diffusionCoeff(compIdx) + boundaryVars1.eddyDiffusivity())
                                               * boundaryVars1.molarDensity();
                        sumDiffusiveFluxes += diffusiveFlux;
                        sumDiffusiveEnergyFlux += diffusiveFlux
                                                  * boundaryVars1.componentEnthalpy(compIdx)
                                                  * FluidSystem::molarMass(compIdx); // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s]
                    }
                }
                sumDiffusiveEnergyFlux -= sumDiffusiveFluxes * boundaryVars1.componentEnthalpy(phaseCompIdx1)
                                          * FluidSystem::molarMass(phaseCompIdx1);
                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        -(convectiveFlux - sumDiffusiveEnergyFlux - conductiveFlux));
            }
            else
            {
                // the energy flux from the Stokes domain
                couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                        globalProblem.localResidual1().residual(vertInElem1)[energyEqIdx1]);
            }
        }
        if (cParams.boundaryTypes2.isCouplingDirichlet(energyEqIdx2))
        {
            // set residualDarcy[energyEqIdx2] = T in 2p2cnilocalresidual.hh
            couplingRes2.accumulate(lfsu2.child(energyEqIdx2), vertInElem2,
                                    -cParams.elemVolVarsCur1[vertInElem1].temperature());
        }
    }

    //! \copydoc Dumux::TwoCStokesTwoPTwoCLocalOperator::evalCoupling21()
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    DUNE_DEPRECATED_MSG("evalCoupling21() is deprecated. Use evalCoupling() instead.")
    void evalCoupling21(const LFSU1& lfsu1, const LFSU2& lfsu2,
                        const int vertInElem1, const int vertInElem2,
                        const SDElement1& sdElement1, const SDElement2& sdElement2,
                        const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                        const CParams &cParams,
                        RES1& couplingRes1, RES2& couplingRes2) const
    {
        GlobalProblem& globalProblem = this->globalProblem();

        // evaluate coupling of mass and momentum balances
        ParentType::evalCoupling21(lfsu1, lfsu2,
                                   vertInElem1, vertInElem2,
                                   sdElement1, sdElement2,
                                   boundaryVars1, boundaryVars2,
                                   cParams,
                                   couplingRes1, couplingRes2);

        const GlobalPosition& globalPos2 = cParams.fvGeometry2.subContVol[vertInElem2].global;
        GlobalPosition normalMassFlux2(0.);

        // velocity*normal*area*rho
        // mass flux is needed for both (mass/mole) formulation, as the enthalpy is mass based
        for (int phaseIdx=0; phaseIdx<numPhases2; ++phaseIdx)
            normalMassFlux2[phaseIdx] = -boundaryVars2.volumeFlux(phaseIdx)*
                cParams.elemVolVarsCur2[vertInElem2].density(phaseIdx);

        if (cParams.boundaryTypes1.isCouplingDirichlet(energyEqIdx1))
        {
            // set residualStokes[energyIdx1] = T in stokes2cnilocalresidual.hh
            couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                    -cParams.elemVolVarsCur2[vertInElem2].temperature());
        }
        if (cParams.boundaryTypes1.isCouplingNeumann(energyEqIdx1))
        {
            if (globalProblem.sdProblem2().isCornerPoint(globalPos2))
            {
                const Scalar convectiveFlux = normalMassFlux2[nPhaseIdx2] *
                    cParams.elemVolVarsCur2[vertInElem2].enthalpy(nPhaseIdx2)
                    +
                    normalMassFlux2[wPhaseIdx2] *
                    cParams.elemVolVarsCur2[vertInElem2].enthalpy(wPhaseIdx2);
                const Scalar conductiveFlux = boundaryVars2.normalMatrixHeatFlux();

                couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                        -(convectiveFlux - conductiveFlux));
            }
            else
            {
                couplingRes1.accumulate(lfsu1.child(energyEqIdx1), vertInElem1,
                                        globalProblem.localResidual2().residual(vertInElem2)[energyEqIdx2]);
            }
        }
    }

    /*!
     * \brief Returns the temperature gradient through the boundary layer
     *
     * \todo This function could be moved to a more model specific place, because
     *       of its runtime parameters.
     *
     * \param cParams a parameter container
     * \param scvIdx1 The local index of the sub-control volume of the Stokes domain
     */
    template<typename CParams>
    Scalar evalBoundaryLayerTemperatureGradient(CParams cParams, const int scvIdx) const
    {
        const Scalar temperatureOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
        Scalar normalTemperatureGrad = cParams.elemVolVarsCur1[scvIdx].temperature()
                                       - temperatureOut;
        return normalTemperatureGrad
               / ParentType::evalBoundaryLayerModel(cParams, scvIdx).thermalBoundaryLayerThickness();
    }
};
} // end namespace Dumux

#endif // DUMUX_TWOCNISTOKES2P2CNILOCALOPERATOR_HH
