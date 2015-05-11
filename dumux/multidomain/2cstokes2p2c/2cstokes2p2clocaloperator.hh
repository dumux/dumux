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
 * \brief The local operator for the coupling of a two-component Stokes model
 *        and a two-phase two-component porous-medium model under isothermal conditions.
 */
#ifndef DUMUX_2CSTOKES_2P2C_LOCALOPERATOR_HH
#define DUMUX_2CSTOKES_2P2C_LOCALOPERATOR_HH

#include <iostream>

#include <dune/common/version.hh>

#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dumux/multidomain/common/multidomainproperties.hh>
#include <dumux/multidomain/2cstokes2p2c/2cstokes2p2cpropertydefaults.hh>
#include <dumux/freeflow/boundarylayermodel.hh>
#include <dumux/freeflow/masstransfermodel.hh>
#include <dumux/freeflow/stokesnc/stokesncmodel.hh>
#include <dumux/implicit/2p2c/2p2cmodel.hh>


namespace Dumux {

/*!
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief The local operator for the coupling of a two-component Stokes model
 *        and a two-phase two-component porous-medium model under isothermal conditions.
 */
template<class TypeTag>
class TwoCStokesTwoPTwoCLocalOperator :
        public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags,
        public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<TwoCStokesTwoPTwoCLocalOperator<TypeTag>>,
        public Dune::PDELab::MultiDomain::FullCouplingPattern,
        public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
 public:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) GlobalProblem;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCouplingLocalOperator) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) Stokes2cTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TwoPTwoCTypeTag;

    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, ElementVolumeVariables) ElementVolumeVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, ElementVolumeVariables) ElementVolumeVariables2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, FluxVariables) BoundaryVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, FluxVariables) BoundaryVariables2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, BoundaryTypes) BoundaryTypes1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, BoundaryTypes) BoundaryTypes2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, FVElementGeometry) FVElementGeometry1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, FVElementGeometry) FVElementGeometry2;

    // Multidomain Grid and Subgrid types
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename MDGrid::Traits::template Codim<0>::Entity MDElement;
    typedef typename MDGrid::SubDomainGrid SDGrid;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, GridView) Stokes2cGridView;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, GridView) TwoPTwoCGridView;
    typedef typename Stokes2cGridView::template Codim<0>::Entity SDElement1;
    typedef typename TwoPTwoCGridView::template Codim<0>::Entity SDElement2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, Indices) Stokes2cIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, Indices) TwoPTwoCIndices;

    enum { dim = MDGrid::dimension };

    // FREE FLOW
    enum { numEq1 = GET_PROP_VALUE(Stokes2cTypeTag, NumEq) };
    enum { nPhaseIdx1 = Stokes2cIndices::phaseIdx };           	//!< Index of the free-flow phase of the fluidsystem
    enum { // equation indices in the Stokes domain
        momentumXIdx1 = Stokes2cIndices::momentumXIdx, 			//!< Index of the x-component of the momentum balance
        momentumYIdx1 = Stokes2cIndices::momentumYIdx, 			//!< Index of the y-component of the momentum balance
        momentumZIdx1 = Stokes2cIndices::momentumZIdx, 			//!< Index of the z-component of the momentum balance
        lastMomentumIdx1 = Stokes2cIndices::lastMomentumIdx, 	//!< Index of the last component of the momentum balance
        massBalanceIdx1 = Stokes2cIndices::massBalanceIdx, 		//!< Index of the mass balance
        transportEqIdx1 = Stokes2cIndices::transportEqIdx 		//!< Index of the transport equation
    };
    enum { // indices of the components
        transportCompIdx1 = Stokes2cIndices::transportCompIdx, 	//!< Index of transported component
        phaseCompIdx1 = Stokes2cIndices::phaseCompIdx         	//!< Index of main component of the phase
    };

    // POROUS MEDIUM
    enum { numEq2 = GET_PROP_VALUE(TwoPTwoCTypeTag, NumEq) };
    enum { numPhases2 = GET_PROP_VALUE(TwoPTwoCTypeTag, NumPhases) };
    enum { // equation indices in the Darcy domain
        contiWEqIdx2 = TwoPTwoCIndices::contiWEqIdx,    //!< Index of the continuity equation for water component
        massBalanceIdx2 = TwoPTwoCIndices::contiNEqIdx  //!< Index of the total mass balance (if one comopnent balance is replaced)
    };
    enum { // component indices
    	wCompIdx2 = TwoPTwoCIndices::wCompIdx,          //!< Index of the liquids main component
        nCompIdx2 = TwoPTwoCIndices::nCompIdx           //!< Index of the main component of the gas
    };
    enum { // phase indices
        wPhaseIdx2 = TwoPTwoCIndices::wPhaseIdx,        //!< Index for the liquid phase
        nPhaseIdx2 = TwoPTwoCIndices::nPhaseIdx          //!< Index for the gas phase
    };

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename Stokes2cGridView::template Codim<dim>::EntityPointer VertexPointer1;
    typedef typename TwoPTwoCGridView::template Codim<dim>::EntityPointer VertexPointer2;
    typedef typename MDGrid::Traits::template Codim<0>::EntityPointer MDElementPointer;

    TwoCStokesTwoPTwoCLocalOperator(GlobalProblem& globalProblem)
        : globalProblem_(globalProblem)
    {
        blModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, BoundaryLayer, Model);
        massTransferModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, MassTransfer, Model);

        if (blModel_ != 0)
            std::cout << "Using boundary layer model " << blModel_ << std::endl;
        if (massTransferModel_ != 0)
            std::cout << "Using mass transfer model " << massTransferModel_ << std::endl;
    }

    // multidomain flags
    static const bool doAlphaCoupling = true;
    static const bool doPatternCoupling = true;

    /*!
     * \brief Do the coupling. The unknowns are transferred from dune-multidomain.
     *        Based on them, a coupling residual is calculated and added at the
     *        respective positions in the matrix.
     *
     * \param intersectionGeometry the geometry of the intersection
     * \param lfsu_s local basis for the trial space of the Stokes domain
     * \param unknowns1 the unknowns vector of the Stokes element (formatted according to PDELab)
     * \param lfsv_s local basis for the test space of the Stokes domain
     * \param lfsu_n local basis for the trail space of the Darcy domain
     * \param unknowns2 the unknowns vector of the Darcy element (formatted according to PDELab)
     * \param lfsv_n local basis for the test space of the Darcy domain
     * \param couplingRes1 the coupling residual from the Stokes domain
     * \param couplingRes2 the coupling residual from the Darcy domain
     *
     */
    template<typename IntersectionGeom, typename LFSU1, typename LFSU2,
             typename X, typename LFSV1, typename LFSV2,typename RES>
    void alpha_coupling (const IntersectionGeom& intersectionGeometry,
                         const LFSU1& lfsu_s, const X& unknowns1, const LFSV1& lfsv_s,
                         const LFSU2& lfsu_n, const X& unknowns2, const LFSV2& lfsv_n,
                         RES& couplingRes1, RES& couplingRes2) const
    {
        const MDElementPointer mdElementPointer1 = intersectionGeometry.inside();
        const MDElementPointer mdElementPointer2 = intersectionGeometry.outside();

        const MDElement& mdElement1 = *mdElementPointer1;
        const MDElement& mdElement2 = *mdElementPointer2;

        // the subodmain elements
        const SDElement1& sdElement1 = *(globalProblem_.sdElementPointer1(mdElement1));
        const SDElement2& sdElement2 = *(globalProblem_.sdElementPointer2(mdElement2));

        // a container for the parameters on each side of the coupling interface (see below)
        CParams cParams;

        // update fvElementGeometry and the element volume variables
        updateElemVolVars(lfsu_s, lfsu_n,
                          unknowns1, unknowns2,
                          sdElement1, sdElement2,
                          cParams);

        // first element
        const int faceIdx1 = intersectionGeometry.indexInInside();
        const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement1 =
            Dune::ReferenceElements<typename MDGrid::ctype,dim>::general(mdElement1.type());
        const int numVerticesOfFace = referenceElement1.size(faceIdx1, 1, dim);

        // second element
        const int faceIdx2 = intersectionGeometry.indexInOutside();
        const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement2 =
            Dune::ReferenceElements<typename MDGrid::ctype,dim>::general(mdElement2.type());

        for (int vertexInFace = 0; vertexInFace < numVerticesOfFace; ++vertexInFace)
        {
            const int vertInElem1 = referenceElement1.subEntity(faceIdx1, 1, vertexInFace, dim);
            const int vertInElem2 = referenceElement2.subEntity(faceIdx2, 1, vertexInFace, dim);

            const int boundaryFaceIdx1 = cParams.fvGeometry1.boundaryFaceIndex(faceIdx1, vertexInFace);
            const int boundaryFaceIdx2 = cParams.fvGeometry2.boundaryFaceIndex(faceIdx2, vertexInFace);

            // obtain the boundary types
            const VertexPointer1 vPtr1 = sdElement1.template subEntity<dim>(vertInElem1);
            const VertexPointer2 vPtr2 = sdElement2.template subEntity<dim>(vertInElem2);

            globalProblem_.sdProblem1().boundaryTypes(cParams.boundaryTypes1, *vPtr1);
            globalProblem_.sdProblem2().boundaryTypes(cParams.boundaryTypes2, *vPtr2);

            const BoundaryVariables1 boundaryVars1(globalProblem_.sdProblem1(),
                                                   sdElement1,
                                                   cParams.fvGeometry1,
                                                   boundaryFaceIdx1,
                                                   cParams.elemVolVarsCur1,
                                                   /*onBoundary=*/true);
            const BoundaryVariables2 boundaryVars2(globalProblem_.sdProblem2(),
                                                   sdElement2,
                                                   cParams.fvGeometry2,
                                                   boundaryFaceIdx2,
                                                   cParams.elemVolVarsCur2,
                                                   /*onBoundary=*/true);

            asImp_()->evalCoupling12(lfsu_s, lfsu_n, // local function spaces
                                     vertInElem1, vertInElem2,
                                     sdElement1, sdElement2,
                                     boundaryVars1, boundaryVars2,
                                     cParams,
                                     couplingRes1, couplingRes2);
            asImp_()->evalCoupling21(lfsu_s, lfsu_n, // local function spaces
                                     vertInElem1, vertInElem2,
                                     sdElement1, sdElement2,
                                     boundaryVars1, boundaryVars2,
                                     cParams,
                                     couplingRes1, couplingRes2);
        }
    }

    /*!
     * \brief Update the volume variables of the element and extract the unknowns from dune-pdelab vectors
     *        and bring them into a form which fits to dumux.
     *
     * \param lfsu_s local basis for the trial space of the Stokes domain
     * \param lfsu_n local basis for the trial space of the Darcy domain
     * \param unknowns1 the unknowns vector of the Stokes element (formatted according to PDELab)
     * \param unknowns2 the unknowns vector of the Darcy element (formatted according to PDELab)
     * \param sdElement1 the element in the Stokes domain
     * \param sdElement2 the element in the Darcy domain
     * \param cParams a parameter container
     *
     */
    template<typename LFSU1, typename LFSU2, typename X, typename CParams>
    void updateElemVolVars (const LFSU1& lfsu_s, const LFSU2& lfsu_n,
                            const X& unknowns1, const X& unknowns2,
                            const SDElement1& sdElement1, const SDElement2& sdElement2,
                            CParams &cParams) const
    {
        cParams.fvGeometry1.update(globalProblem_.sdGridView1(), sdElement1);
        cParams.fvGeometry2.update(globalProblem_.sdGridView2(), sdElement2);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        const int numVertsOfElem1 = sdElement1.subEntities(dim);
        const int numVertsOfElem2 = sdElement2.subEntities(dim);
#else
        const int numVertsOfElem1 = sdElement1.template count<dim>();
        const int numVertsOfElem2 = sdElement2.template count<dim>();
#endif

        //bring the local unknowns x_s into a form that can be passed to elemVolVarsCur.update()
        Dune::BlockVector<Dune::FieldVector<Scalar,1>> elementSol1(0.);
        Dune::BlockVector<Dune::FieldVector<Scalar,1>> elementSol2(0.);
        elementSol1.resize(unknowns1.size());
        elementSol2.resize(unknowns2.size());

        for (int idx=0; idx<numVertsOfElem1; ++idx)
        {
            for (int eqIdx1=0; eqIdx1<numEq1; ++eqIdx1)
                elementSol1[eqIdx1*numVertsOfElem1+idx] = unknowns1(lfsu_s.child(eqIdx1),idx);
            for (int eqIdx2=0; eqIdx2<numEq2; ++eqIdx2)
                elementSol2[eqIdx2*numVertsOfElem2+idx] = unknowns2(lfsu_n.child(eqIdx2),idx);
        }
#if HAVE_VALGRIND
        for (unsigned int i = 0; i < elementSol1.size(); i++)
            Valgrind::CheckDefined(elementSol1[i]);
        for (unsigned int i = 0; i < elementSol2.size(); i++)
            Valgrind::CheckDefined(elementSol2[i]);
#endif // HAVE_VALGRIND


        // evaluate the local residual with the PDELab solution
        globalProblem_.localResidual1().evalPDELab(sdElement1, cParams.fvGeometry1, elementSol1,
                                                   cParams.elemVolVarsPrev1, cParams.elemVolVarsCur1);
        globalProblem_.localResidual2().evalPDELab(sdElement2, cParams.fvGeometry2, elementSol2,
                                                   cParams.elemVolVarsPrev2, cParams.elemVolVarsCur2);

    }

    /*!
     * \brief Evaluation of the coupling from Stokes (1 or s) to Darcy (2 or n).
     *
     * \param lfsu_s local basis for the trial space of the Stokes domain
     * \param lfsu_n local basis for the trial space of the Darcy domain
     * \param vertInElem1 local vertex index in element1
     * \param vertInElem2 local vertex index in element2
     * \param sdElement1 the element in the Stokes domain
     * \param sdElement2 the element in the Darcy domain
     * \param boundaryVars1 the boundary variables at the interface of the Stokes domain
     * \param boundaryVars2 the boundary variables at the interface of the Darcy domain
     * \param cParams a parameter container
     * \param couplingRes1 the coupling residual from the Stokes domain
     * \param couplingRes2 the coupling residual from the Darcy domain
     */
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
        const Scalar normalMassFlux = boundaryVars1.normalVelocity() *
            cParams.elemVolVarsCur1[vertInElem1].density();

        //rho*v*n as NEUMANN condition for porous medium (set, if B&J defined as NEUMANN condition)
        if (cParams.boundaryTypes2.isCouplingInflow(massBalanceIdx2))
        {
            if (globalProblem_.sdProblem1().isCornerPoint(globalPos1))
            {
                couplingRes2.accumulate(lfsu_n.child(massBalanceIdx2), vertInElem2,
                                        -normalMassFlux);
            }
            else
            {
                couplingRes2.accumulate(lfsu_n.child(massBalanceIdx2), vertInElem2,
                                        globalProblem_.localResidual1().residual(vertInElem1)[massBalanceIdx1]);
            }
        }
        if (cParams.boundaryTypes2.isCouplingOutflow(massBalanceIdx2))
        {
            couplingRes2.accumulate(lfsu_n.child(massBalanceIdx2), vertInElem2,
                                    globalProblem_.localResidual1().residual(vertInElem1)[momentumYIdx1]
                                    -cParams.elemVolVarsCur1[vertInElem1].pressure());
        }

        if (cParams.boundaryTypes2.isCouplingInflow(contiWEqIdx2))
        {
            // only enter here, if a BOUNDARY LAYER MODEL is used for the computation of the diffusive fluxes
            if (blModel_)
            {
                const Scalar massFractionOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
                const Scalar M1 = FluidSystem::molarMass(transportCompIdx1);
                const Scalar M2 = FluidSystem::molarMass(phaseCompIdx1);
                const Scalar X2 = 1.0 - massFractionOut;
                const Scalar massToMoleDenominator = M2 + X2*(M1 - M2);
                const Scalar moleFractionOut = massFractionOut * M2 /massToMoleDenominator;

                const Scalar advectiveFlux =
                    normalMassFlux *
                    cParams.elemVolVarsCur1[vertInElem1].massFraction(transportCompIdx1);
                Scalar normalMoleFracGrad =
                    cParams.elemVolVarsCur1[vertInElem1].moleFraction(transportCompIdx1) -
                    moleFractionOut;

                const Scalar velocity = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);
                // current position + additional virtual runup distance
                const Scalar distance = globalPos1[0] + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, Offset);
                const Scalar kinematicViscosity = cParams.elemVolVarsCur1[vertInElem2].kinematicViscosity();
                BoundaryLayerModel<TypeTag> boundaryLayerModel(velocity, distance, kinematicViscosity, blModel_);
                if (blModel_ == 1)
                    boundaryLayerModel.setConstThickness(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, ConstThickness));
                if (blModel_ >= 4)
                    boundaryLayerModel.setYPlus(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, YPlus));
                if (blModel_ >= 5)
                    boundaryLayerModel.setRoughnessLength(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, RoughnessLength));

                normalMoleFracGrad /= boundaryLayerModel.massBoundaryLayerThickness();

                const Scalar diffusiveFlux =
                    bfNormal1.two_norm() *
                    normalMoleFracGrad *
                    (boundaryVars1.diffusionCoeff(transportCompIdx1) + boundaryVars1.eddyDiffusivity()) *
                    boundaryVars1.molarDensity() *
                    FluidSystem::molarMass(transportCompIdx1);

                if (massTransferModel_ == 0)
                {
                    couplingRes2.accumulate(lfsu_n.child(contiWEqIdx2), vertInElem2,
                                            -(advectiveFlux - diffusiveFlux));
                }
                // transition from the mass transfer coefficient concept to the coupling via the local residual,
                // when saturations become small; only diffusive fluxes are scaled!
                else
                {
                    MassTransferModel<TypeTag> massTransferModel(cParams.elemVolVarsCur2[vertInElem2].saturation(wPhaseIdx2),
                                                    cParams.elemVolVarsCur2[vertInElem2].porosity(),
                                                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, CharPoreDiameter),
                                                    boundaryLayerModel.massBoundaryLayerThickness(),
                                                    massTransferModel_);
                    if (massTransferModel_ == 1)
                        massTransferModel.setMassTransferCoeff(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, Coefficient));
                    if (massTransferModel_ == 3)
                        massTransferModel.setCapillaryPressure(cParams.elemVolVarsCur2[vertInElem2].capillaryPressure());
                    const Scalar massTransferCoeff = massTransferModel.massTransferCoefficient();

                    if (massTransferCoeff > 1.0 || massTransferCoeff < 0.0)
                        std::cout << "MTC out of bounds, should be in between 0.0 and 1.0! >>> " << massTransferCoeff << std::endl;

                    if (globalProblem_.sdProblem1().isCornerPoint(globalPos1))
                    {
                        const Scalar diffusiveFluxAtCorner =
                            bfNormal1 *
                            boundaryVars1.moleFractionGrad(transportCompIdx1) *
                            (boundaryVars1.diffusionCoeff(transportCompIdx1) + boundaryVars1.eddyDiffusivity()) *
                            boundaryVars1.molarDensity() *
                            FluidSystem::molarMass(transportCompIdx1);

                        couplingRes2.accumulate(lfsu_n.child(contiWEqIdx2), vertInElem2,
                                                -massTransferCoeff*(advectiveFlux - diffusiveFlux) -
                                                (1.-massTransferCoeff)*(advectiveFlux - diffusiveFluxAtCorner));
                    }
                    else
                    {
                        couplingRes2.accumulate(lfsu_n.child(contiWEqIdx2), vertInElem2,
                                                -massTransferCoeff*(advectiveFlux - diffusiveFlux) +
                                                (1.-massTransferCoeff)*globalProblem_.localResidual1().residual(vertInElem1)[transportEqIdx1]);
                    }
                }
            }
            else
            {
                // compute fluxes explicitly at corner points - only quarter control volume
                if (globalProblem_.sdProblem1().isCornerPoint(globalPos1))
                {
                    const Scalar advectiveFlux =
                        normalMassFlux *
                        cParams.elemVolVarsCur1[vertInElem1].massFraction(transportCompIdx1);
                    const Scalar diffusiveFlux =
                        bfNormal1 *
                        boundaryVars1.moleFractionGrad(transportCompIdx1) *
                        (boundaryVars1.diffusionCoeff(transportCompIdx1) + boundaryVars1.eddyDiffusivity()) *
                        boundaryVars1.molarDensity() *
                        FluidSystem::molarMass(transportCompIdx1);

                    couplingRes2.accumulate(lfsu_n.child(contiWEqIdx2), vertInElem2,
                                            -(advectiveFlux - diffusiveFlux));
                }
                // coupling via the defect
                else
                {
                    // the component mass flux from the stokes domain
                    couplingRes2.accumulate(lfsu_n.child(contiWEqIdx2), vertInElem2,
                                            globalProblem_.localResidual1().residual(vertInElem1)[transportEqIdx1]);
                }
            }
        }
        if (cParams.boundaryTypes2.isCouplingOutflow(contiWEqIdx2))
        	std::cerr << "Upwind PM -> FF does not work for the transport equation for a 2-phase system!" << std::endl;
    }

    /*!
     * \brief Evaluation of the coupling from Darcy (2 or n) to Stokes (1 or s).
     *
     * \param lfsu_s local basis for the trial space of the Stokes domain
     * \param lfsu_n local basis for the trial space of the Darcy domain
     * \param vertInElem1 local vertex index in element1
     * \param vertInElem2 local vertex index in element2
     * \param sdElement1 the element in the Stokes domain
     * \param sdElement2 the element in the Darcy domain
     * \param boundaryVars1 the boundary variables at the interface of the Stokes domain
     * \param boundaryVars2 the boundary variables at the interface of the Darcy domain
     * \param cParams a parameter container
     * \param couplingRes1 the coupling residual from the Stokes domain
     * \param couplingRes2 the coupling residual from the Darcy domain
     */
    template<typename LFSU1, typename LFSU2, typename RES1, typename RES2, typename CParams>
    void evalCoupling21(const LFSU1& lfsu_s, const LFSU2& lfsu_n,
                        const int vertInElem1, const int vertInElem2,
                        const SDElement1& sdElement1, const SDElement2& sdElement2,
                        const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                        const CParams &cParams,
                        RES1& couplingRes1, RES2& couplingRes2) const
    {
        const DimVector& globalPos2 = cParams.fvGeometry2.subContVol[vertInElem2].global;
        DimVector normalFlux2(0.);

        // velocity*normal*area*rho
        for (int phaseIdx=0; phaseIdx<numPhases2; ++phaseIdx)
            normalFlux2[phaseIdx] = -boundaryVars2.volumeFlux(phaseIdx)*
                cParams.elemVolVarsCur2[vertInElem2].density(phaseIdx);

        //p*n as NEUMANN condition for free flow (set, if B&J defined as NEUMANN condition)
        if (cParams.boundaryTypes1.isCouplingOutflow(momentumYIdx1))
        {
            //p*A*n as NEUMANN condition for free flow (set, if B&J defined as NEUMANN condition)
            //pressure correction in stokeslocalresidual.hh
            couplingRes1.accumulate(lfsu_s.child(momentumYIdx1), vertInElem1,
                                    cParams.elemVolVarsCur2[vertInElem2].pressure(nPhaseIdx2) *
                                    boundaryVars2.face().area);
        }
        if (cParams.boundaryTypes1.isCouplingInflow(momentumYIdx1))
        {

            // v.n as Dirichlet condition for the Stokes domain
            // set residualStokes[momentumYIdx1] = vy in stokeslocalresidual.hh
            if (globalProblem_.sdProblem2().isCornerPoint(globalPos2))
            {
                couplingRes1.accumulate(lfsu_s.child(momentumYIdx1), vertInElem1,
                                        -((normalFlux2[nPhaseIdx2] + normalFlux2[wPhaseIdx2])
                                          / cParams.elemVolVarsCur1[vertInElem1].density()));
            }
            else
            {
                // v.n as DIRICHLET condition for the Stokes domain (negative sign!)
                couplingRes1.accumulate(lfsu_s.child(momentumYIdx1), vertInElem1,
                                        globalProblem_.localResidual2().residual(vertInElem2)[massBalanceIdx2]
                                        / cParams.elemVolVarsCur1[vertInElem1].density());
            }
        }

        //coupling residual is added to "real" residual
        //here each node is passed twice, hence only half of the dirichlet condition has to be set
        if (cParams.boundaryTypes1.isCouplingOutflow(transportEqIdx1))
        {
            // set residualStokes[transportEqIdx1] = x in stokes2clocalresidual.hh
            couplingRes1.accumulate(lfsu_s.child(transportEqIdx1), vertInElem1,
                                    -cParams.elemVolVarsCur2[vertInElem2].massFraction(nPhaseIdx2, wCompIdx2));
        }
        if (cParams.boundaryTypes1.isCouplingInflow(transportEqIdx1))
        	std::cerr << "Upwind PM -> FF does not work for the transport equation for a 2-phase system!" << std::endl;
    }

 protected:
    GlobalProblem& globalProblem() const
    { return globalProblem_; }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

 private:
    /*!
     * \brief A struct that contains data of the FF and PM including boundary types,
     *        volume variables in both subdomains and geometric information
     */
    struct CParams
    {
        BoundaryTypes1 boundaryTypes1;
        BoundaryTypes2 boundaryTypes2;
        ElementVolumeVariables1 elemVolVarsPrev1;
        ElementVolumeVariables1 elemVolVarsCur1;
        ElementVolumeVariables2 elemVolVarsPrev2;
        ElementVolumeVariables2 elemVolVarsCur2;
        FVElementGeometry1 fvGeometry1;
        FVElementGeometry2 fvGeometry2;
    };

    unsigned int blModel_;
    unsigned int massTransferModel_;

    GlobalProblem& globalProblem_;
};

} // end namespace Dumux

#endif // DUMUX_2CSTOKES_2P2C_LOCALOPERATOR_HH
