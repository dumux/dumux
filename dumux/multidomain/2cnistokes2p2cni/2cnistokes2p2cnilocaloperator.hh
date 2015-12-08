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
    typedef typename GET_PROP_TYPE(TypeTag, Problem) GlobalProblem;

    // Get the TypeTags of the subproblems
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) Stokes2cniTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TwoPTwoCNITypeTag;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, FluxVariables) BoundaryVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, FluxVariables) BoundaryVariables2;

    // Multidomain Grid and Subgrid types
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, GridView) Stokes2cniGridView;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, GridView) TwoPTwoCNIGridView;

    typedef typename Stokes2cniGridView::template Codim<0>::Entity SDElement1;
    typedef typename TwoPTwoCNIGridView::template Codim<0>::Entity SDElement2;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, Indices) Stokes2cniIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, Indices) TwoPTwoCNIIndices;

    enum { dim = MDGrid::dimension };
    enum {
        energyEqIdx1 = Stokes2cniIndices::energyEqIdx          //!< Index of the energy balance equation
    };
    enum {  numComponents = Stokes2cniIndices::numComponents };
    enum { phaseCompIdx = Stokes2cniIndices::phaseCompIdx};
    enum { // indices in the Darcy domain
        numPhases2 = GET_PROP_VALUE(TwoPTwoCNITypeTag, NumPhases),

        // equation index
        energyEqIdx2 = TwoPTwoCNIIndices::energyEqIdx,      //!< Index of the energy balance equation

        wPhaseIdx2 = TwoPTwoCNIIndices::wPhaseIdx,          //!< Index for the liquid phase
        nPhaseIdx2 = TwoPTwoCNIIndices::nPhaseIdx           //!< Index for the gas phase
    };
    enum { nPhaseIdx1 = Stokes2cniIndices::phaseIdx };            //!< Index of the free-flow phase of the fluidsystem
    enum { // indices of the components
        transportCompIdx1 = Stokes2cniIndices::transportCompIdx,  //!< Index of transported component
        phaseCompIdx1 = Stokes2cniIndices::phaseCompIdx           //!< Index of main component of the phase
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;       //!< A field vector with dim entries
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef TwoCStokesTwoPTwoCLocalOperator<TypeTag> ParentType;

    TwoCNIStokesTwoPTwoCNILocalOperator(GlobalProblem& globalProblem)
        : ParentType(globalProblem)
    { }

    static const bool doAlphaCoupling = true;
    static const bool doPatternCoupling = true;

    //! \copydoc Dumux::TwoCStokesTwoPTwoCLocalOperator::evalCoupling12()
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    void evalCoupling12(const LFSU1& lfsu_s, const LFSU2& lfsu_n,
                        const int vertInElem1, const int vertInElem2,
                        const SDElement1& sdElement1, const SDElement2& sdElement2,
                        const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                        const CParams &cParams,
                        RES1& couplingRes1, RES2& couplingRes2) const
    {
        const DimVector& globalPos1 = cParams.fvGeometry1.subContVol[vertInElem1].global;
        const DimVector& bfNormal1 = boundaryVars1.face().normal;
        const Scalar normalMassFlux1 = boundaryVars1.normalVelocity() *
            cParams.elemVolVarsCur1[vertInElem1].density();
        GlobalProblem& globalProblem = this->globalProblem();

        // evaluate coupling of mass and momentum balances
        ParentType::evalCoupling12(lfsu_s, lfsu_n,
                                   vertInElem1, vertInElem2,
                                   sdElement1, sdElement2,
                                   boundaryVars1, boundaryVars2,
                                   cParams,
                                   couplingRes1, couplingRes2);

        if (cParams.boundaryTypes2.isCouplingInflow(energyEqIdx2))
        {
            unsigned int blModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, BoundaryLayer, Model);
            // only enter here, if a boundary layer model  is used for the computation of the diffusive fluxes
            if (blModel)
            {
                // convective energy flux
                const Scalar convectiveFlux =
                    normalMassFlux1 *
                    cParams.elemVolVarsCur1[vertInElem1].enthalpy();

                // diffusive transported energy
                const Scalar massFractionOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
                const Scalar M1 = FluidSystem::molarMass(transportCompIdx1);
                const Scalar M2 = FluidSystem::molarMass(phaseCompIdx1);
                const Scalar X2 = 1.0 - massFractionOut;
                const Scalar massToMoleDenominator = M2 + X2*(M1 - M2);
                const Scalar moleFractionOut = massFractionOut * M2 /massToMoleDenominator;
                Scalar normalMoleFracGrad =
                    cParams.elemVolVarsCur1[vertInElem1].moleFraction(transportCompIdx1) -
                    moleFractionOut;

                const Scalar velocity = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);
                // current position + additional virtual runup distance
                const Scalar distance = globalPos1[0] + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, Offset);
                const Scalar kinematicViscosity = cParams.elemVolVarsCur1[vertInElem2].kinematicViscosity();
                BoundaryLayerModel<TypeTag> boundaryLayerModel(velocity, distance, kinematicViscosity, blModel);
                if (blModel == 1)
                    boundaryLayerModel.setConstThickness(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, ConstThickness));
                if (blModel >= 4)
                    boundaryLayerModel.setYPlus(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, YPlus));
                if (blModel >= 5)
                    boundaryLayerModel.setRoughnessLength(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, RoughnessLength));
                if (blModel == 7)
                    boundaryLayerModel.setHydraulicDiameter(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, HydraulicDiameter));

                normalMoleFracGrad /= boundaryLayerModel.massBoundaryLayerThickness();

                Scalar diffusiveFlux =
                    bfNormal1.two_norm() *
                    normalMoleFracGrad *
                    (boundaryVars1.diffusionCoeff(transportCompIdx1) + boundaryVars1.eddyDiffusivity()) *
                    boundaryVars1.molarDensity();

                // multiply the diffusive flux with the mass transfer coefficient
                unsigned int mtModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, MassTransfer, Model);
                if (mtModel != 0)
                {
                    MassTransferModel<TypeTag> massTransferModel(cParams.elemVolVarsCur2[vertInElem2].saturation(wPhaseIdx2),
                                                                 cParams.elemVolVarsCur2[vertInElem2].porosity(),
                                                                 boundaryLayerModel.massBoundaryLayerThickness(),
                                                                 mtModel);
                    if (mtModel == 1)
                        massTransferModel.setMassTransferCoeff(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, Coefficient));
                    if (mtModel == 2 || mtModel == 4)
                        massTransferModel.setCharPoreRadius(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, CharPoreRadius));
                    if (mtModel == 3)
                        massTransferModel.setCapillaryPressure(cParams.elemVolVarsCur2[vertInElem2].capillaryPressure());
                    diffusiveFlux *= massTransferModel.massTransferCoefficient();
                }

                // the boundary layer models only work for two components
                assert(numComponents == 2);
                Scalar diffusiveEnergyFlux = 0.0;
                diffusiveEnergyFlux += diffusiveFlux * FluidSystem::molarMass(transportCompIdx1)
                                       * boundaryVars1.componentEnthalpy(transportCompIdx1);
                diffusiveEnergyFlux -= diffusiveFlux * FluidSystem::molarMass(phaseCompIdx1)
                                       * boundaryVars1.componentEnthalpy(phaseCompIdx1);

                // conductive transported energy
                const Scalar temperatureOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
                Scalar normalTemperatureGrad =
                    cParams.elemVolVarsCur1[vertInElem1].temperature() -
                    temperatureOut;

                normalTemperatureGrad /= boundaryLayerModel.thermalBoundaryLayerThickness();

                const Scalar conductiveFlux =
                    bfNormal1.two_norm() *
                    normalTemperatureGrad *
                    (boundaryVars1.thermalConductivity() + boundaryVars1.thermalEddyConductivity());

                couplingRes2.accumulate(lfsu_n.child(energyEqIdx2), vertInElem2,
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
                for (int compIdx=0; compIdx < numComponents; compIdx++)
                {
                    if (compIdx != phaseCompIdx)
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
                sumDiffusiveEnergyFlux -= sumDiffusiveFluxes * boundaryVars1.componentEnthalpy(phaseCompIdx)
                                          * FluidSystem::molarMass(phaseCompIdx);
                couplingRes2.accumulate(lfsu_n.child(energyEqIdx2), vertInElem2,
                                        -(convectiveFlux - sumDiffusiveEnergyFlux - conductiveFlux));
            }
            else
            {
                // the energy flux from the stokes domain
                couplingRes2.accumulate(lfsu_n.child(energyEqIdx2), vertInElem2,
                                        globalProblem.localResidual1().residual(vertInElem1)[energyEqIdx1]);
            }
        }
        if (cParams.boundaryTypes2.isCouplingOutflow(energyEqIdx2))
        {
            // set residualDarcy[energyEqIdx2] = T in 2p2cnilocalresidual.hh
            couplingRes2.accumulate(lfsu_n.child(energyEqIdx2), vertInElem2,
                                    -cParams.elemVolVarsCur1[vertInElem1].temperature());
        }
    }

    //! \copydoc Dumux::TwoCStokesTwoPTwoCLocalOperator::evalCoupling21()
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    void evalCoupling21(const LFSU1& lfsu_s, const LFSU2& lfsu_n,
                        const int vertInElem1, const int vertInElem2,
                        const SDElement1& sdElement1, const SDElement2& sdElement2,
                        const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                        const CParams &cParams,
                        RES1& couplingRes1, RES2& couplingRes2) const
    {
        GlobalProblem& globalProblem = this->globalProblem();

        // evaluate coupling of mass and momentum balances
        ParentType::evalCoupling21(lfsu_s, lfsu_n,
                                   vertInElem1, vertInElem2,
                                   sdElement1, sdElement2,
                                   boundaryVars1, boundaryVars2,
                                   cParams,
                                   couplingRes1, couplingRes2);

        const DimVector& globalPos2 = cParams.fvGeometry2.subContVol[vertInElem2].global;
        DimVector normalMassFlux2(0.);

        // velocity*normal*area*rho
        // mass flux is needed for both (mass/mole) formulation, as the enthalpy is mass based
        for (int phaseIdx=0; phaseIdx<numPhases2; ++phaseIdx)
            normalMassFlux2[phaseIdx] = -boundaryVars2.volumeFlux(phaseIdx)*
                cParams.elemVolVarsCur2[vertInElem2].density(phaseIdx);

        if (cParams.boundaryTypes1.isCouplingOutflow(energyEqIdx1))
        {
            // set residualStokes[energyIdx1] = T in stokes2cnilocalresidual.hh
            couplingRes1.accumulate(lfsu_s.child(energyEqIdx1), vertInElem1,
                                    -cParams.elemVolVarsCur2[vertInElem2].temperature());
        }
        if (cParams.boundaryTypes1.isCouplingInflow(energyEqIdx1))
        {
            if (globalProblem.sdProblem2().isCornerPoint(globalPos2))
            {
                const Scalar convectiveFlux = normalMassFlux2[nPhaseIdx2] *
                    cParams.elemVolVarsCur2[vertInElem2].enthalpy(nPhaseIdx2)
                    +
                    normalMassFlux2[wPhaseIdx2] *
                    cParams.elemVolVarsCur2[vertInElem2].enthalpy(wPhaseIdx2);
                const Scalar conductiveFlux = boundaryVars2.normalMatrixHeatFlux();

                couplingRes1.accumulate(lfsu_s.child(energyEqIdx1), vertInElem1,
                                        -(convectiveFlux - conductiveFlux));
            }
            else
            {
                couplingRes1.accumulate(lfsu_s.child(energyEqIdx1), vertInElem1,
                                        globalProblem.localResidual2().residual(vertInElem2)[energyEqIdx2]);
            }
        }
    }
};
} // end namespace Dumux

#endif // DUMUX_TWOCNISTOKES2P2CNILOCALOPERATOR_HH
