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

#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dumux/freeflow/boundarylayermodel.hh>
#include <dumux/freeflow/masstransfermodel.hh>
#include <dumux/freeflow/stokesnc/model.hh>
#include <dumux/multidomain/properties.hh>
#include <dumux/porousmediumflow/2p2c/implicit/model.hh>

#include "propertydefaults.hh"

namespace Dumux {

/*!
 * \ingroup TwoPTwoCStokesTwoCModel
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief The local operator for the coupling of a two-component Stokes model
 *        and a two-phase two-component porous-medium model under isothermal conditions.
 *
 * This model implements the coupling between a free-flow model
 * and a porous-medium flow model under isothermal conditions.
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
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
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

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, SpatialParams) SpatialParams;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, ElementVolumeVariables) ElementVolumeVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, ElementVolumeVariables) ElementVolumeVariables2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, FluxVariables) BoundaryVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, FluxVariables) BoundaryVariables2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, ElementBoundaryTypes) ElementBoundaryTypes1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, ElementBoundaryTypes) ElementBoundaryTypes2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, BoundaryTypes) BoundaryTypes1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, BoundaryTypes) BoundaryTypes2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, FVElementGeometry) FVElementGeometry1;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, FVElementGeometry) FVElementGeometry2;

    // Multidomain Grid types
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename MDGrid::Traits::template Codim<0>::Entity MDElement;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, GridView) Stokes2cGridView;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, GridView) TwoPTwoCGridView;
    typedef typename Stokes2cGridView::template Codim<0>::Entity SDElement1;
    typedef typename TwoPTwoCGridView::template Codim<0>::Entity SDElement2;

    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, Indices) Stokes2cIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, Indices) TwoPTwoCIndices;

    enum {
        dim = MDGrid::dimension,
        dimWorld = MDGrid::dimensionworld
    };

    // Stokes
    enum { numEq1 = GET_PROP_VALUE(Stokes2cTypeTag, NumEq) };
    enum { numComponents1 = Stokes2cIndices::numComponents };
    enum { // equation indices
        momentumXIdx1 = Stokes2cIndices::momentumXIdx,         //!< Index of the x-component of the momentum balance
        momentumYIdx1 = Stokes2cIndices::momentumYIdx,         //!< Index of the y-component of the momentum balance
        momentumZIdx1 = Stokes2cIndices::momentumZIdx,         //!< Index of the z-component of the momentum balance
        lastMomentumIdx1 = Stokes2cIndices::lastMomentumIdx,   //!< Index of the last component of the momentum balance
        massBalanceIdx1 = Stokes2cIndices::massBalanceIdx,     //!< Index of the mass balance
        transportEqIdx1 = Stokes2cIndices::transportEqIdx      //!< Index of the transport equation
    };
    enum { // component indices
        transportCompIdx1 = Stokes2cIndices::transportCompIdx, //!< Index of transported component
        phaseCompIdx1 = Stokes2cIndices::phaseCompIdx    //!< Index of main component of the phase
    };

    // Darcy
    enum { numEq2 = GET_PROP_VALUE(TwoPTwoCTypeTag, NumEq) };
    enum { numPhases2 = GET_PROP_VALUE(TwoPTwoCTypeTag, NumPhases) };
    enum { // equation indices
        contiWEqIdx2 = TwoPTwoCIndices::contiWEqIdx,     //!< Index of the continuity equation for water component
        massBalanceIdx2 = TwoPTwoCIndices::contiNEqIdx   //!< Index of the total mass balance (if one component balance is replaced)
    };
    enum { // component indices
        wCompIdx2 = TwoPTwoCIndices::wCompIdx,           //!< Index of the liquids main component
        nCompIdx2 = TwoPTwoCIndices::nCompIdx            //!< Index of the main component of the gas
    };
    enum { // phase indices
        wPhaseIdx2 = TwoPTwoCIndices::wPhaseIdx,         //!< Index for the liquid phase
        nPhaseIdx2 = TwoPTwoCIndices::nPhaseIdx          //!< Index for the gas phase
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename MDGrid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename Stokes2cGridView::template Codim<dim>::EntityPointer VertexPointer1;
    typedef typename TwoPTwoCGridView::template Codim<dim>::EntityPointer VertexPointer2;

    // multidomain flags
    static const bool doAlphaCoupling = true;
    static const bool doPatternCoupling = true;

public:
    //! \brief The constructor
    TwoCStokesTwoPTwoCLocalOperator(GlobalProblem& globalProblem)
        : globalProblem_(globalProblem)
    {
        static_assert(GET_PROP_VALUE(Stokes2cTypeTag, UseMoles) == GET_PROP_VALUE(TwoPTwoCTypeTag, UseMoles),
                      "The coupling conditions is only implemented for same formulations (mass or mole) in both subdomains.");

        blModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, BoundaryLayer, Model);
        massTransferModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, MassTransfer, Model);

        if (blModel_ != 0)
            std::cout << "Using boundary layer model " << blModel_ << std::endl;
        if (massTransferModel_ != 0)
            std::cout << "Using mass transfer model " << massTransferModel_ << std::endl;
    }

    /*!
     * \brief Do the coupling. The unknowns are transferred from dune-multidomain.
     *        Based on them, a coupling residual is calculated and added at the
     *        respective positions in the matrix.
     *
     * \param intersectionGeometry the geometry of the intersection
     * \param lfsu1 local basis for the trial space of the Stokes domain
     * \param unknowns1 the unknowns vector of the Stokes element (formatted according to PDELab)
     * \param lfsv1 local basis for the test space of the Stokes domain
     * \param lfsu2 local basis for the trail space of the Darcy domain
     * \param unknowns2 the unknowns vector of the Darcy element (formatted according to PDELab)
     * \param lfsv2 local basis for the test space of the Darcy domain
     * \param couplingRes1 the coupling residual from the Stokes domain
     * \param couplingRes2 the coupling residual from the Darcy domain
     */
    template<typename IntersectionGeom, typename LFSU1, typename LFSU2,
             typename X, typename LFSV1, typename LFSV2,typename RES>
    void alpha_coupling(const IntersectionGeom& intersectionGeometry,
                        const LFSU1& lfsu1, const X& unknowns1, const LFSV1& lfsv1,
                        const LFSU2& lfsu2, const X& unknowns2, const LFSV2& lfsv2,
                        RES& couplingRes1, RES& couplingRes2) const
    {
        const std::shared_ptr<MDElement> mdElement1
            = std::make_shared<MDElement>(intersectionGeometry.inside());
        const std::shared_ptr<MDElement> mdElement2
            = std::make_shared<MDElement>(intersectionGeometry.outside());

        // the subdomain elements
        const std::shared_ptr<SDElement1> sdElement1
            = std::make_shared<SDElement1>(globalProblem_.sdElementPointer1(*mdElement1));
        const std::shared_ptr<SDElement2> sdElement2
            = std::make_shared<SDElement2>(globalProblem_.sdElementPointer2(*mdElement2));

        // a container for the parameters on each side of the coupling interface (see below)
        CParams cParams;

        // update fvElementGeometry and the element volume variables
        updateElemVolVars(lfsu1, lfsu2,
                          unknowns1, unknowns2,
                          *sdElement1, *sdElement2,
                          cParams);

        // first element
        const int faceIdx1 = intersectionGeometry.indexInInside();
        const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement1 =
            Dune::ReferenceElements<typename MDGrid::ctype,dim>::general((*mdElement1).type());
        const int numVerticesOfFace = referenceElement1.size(faceIdx1, 1, dim);

        // second element
        const int faceIdx2 = intersectionGeometry.indexInOutside();
        const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement2 =
            Dune::ReferenceElements<typename MDGrid::ctype,dim>::general((*mdElement2).type());

        for (int vertexInFace = 0; vertexInFace < numVerticesOfFace; ++vertexInFace)
        {
            const int vertInElem1 = referenceElement1.subEntity(faceIdx1, 1, vertexInFace, dim);
            const int vertInElem2 = referenceElement2.subEntity(faceIdx2, 1, vertexInFace, dim);

            const int boundaryFaceIdx1 = cParams.fvGeometry1.boundaryFaceIndex(faceIdx1, vertexInFace);
            const int boundaryFaceIdx2 = cParams.fvGeometry2.boundaryFaceIndex(faceIdx2, vertexInFace);

            // obtain the boundary types
            const VertexPointer1 vPtr1 = (*sdElement1).template subEntity<dim>(vertInElem1);
            const VertexPointer2 vPtr2 = (*sdElement2).template subEntity<dim>(vertInElem2);

            globalProblem_.sdProblem1().boundaryTypes(cParams.boundaryTypes1, vPtr1);
            globalProblem_.sdProblem2().boundaryTypes(cParams.boundaryTypes2, vPtr2);

            BoundaryVariables1 boundaryVars1;
            boundaryVars1.update(globalProblem_.sdProblem1(),
                                 *sdElement1,
                                 cParams.fvGeometry1,
                                 boundaryFaceIdx1,
                                 cParams.elemVolVarsCur1,
                                 /*onBoundary=*/true);
            BoundaryVariables2 boundaryVars2;
            boundaryVars2.update(globalProblem_.sdProblem2(),
                                 *sdElement2,
                                 cParams.fvGeometry2,
                                 boundaryFaceIdx2,
                                 cParams.elemVolVarsCur2,
                                 /*onBoundary=*/true);

            asImp_()->evalCoupling(lfsu1, lfsu2,
                                   vertInElem1, vertInElem2,
                                   *sdElement1, *sdElement2,
                                   boundaryVars1, boundaryVars2,
                                   cParams,
                                   couplingRes1, couplingRes2);
        }
    }

    /*!
     * \brief Update the volume variables of the element and extract the unknowns from dune-pdelab vectors
     *        and bring them into a form which fits to dumux.
     *
     * \param lfsu1 local basis for the trial space of the Stokes domain
     * \param lfsu2 local basis for the trial space of the Darcy domain
     * \param unknowns1 the unknowns vector of the Stokes element (formatted according to PDELab)
     * \param unknowns2 the unknowns vector of the Darcy element (formatted according to PDELab)
     * \param sdElement1 the element in the Stokes domain
     * \param sdElement2 the element in the Darcy domain
     * \param cParams a parameter container
     */
    template<typename LFSU1, typename LFSU2, typename X, typename CParams>
    void updateElemVolVars(const LFSU1& lfsu1, const LFSU2& lfsu2,
                           const X& unknowns1, const X& unknowns2,
                           const SDElement1& sdElement1, const SDElement2& sdElement2,
                           CParams &cParams) const
    {
        cParams.fvGeometry1.update(globalProblem_.sdGridView1(), sdElement1);
        cParams.fvGeometry2.update(globalProblem_.sdGridView2(), sdElement2);

        const int numVertsOfElem1 = sdElement1.subEntities(dim);
        const int numVertsOfElem2 = sdElement2.subEntities(dim);

        // bring the local unknowns x1 into a form that can be passed to elemVolVarsCur.update()
        Dune::BlockVector<Dune::FieldVector<Scalar,1>> elementSol1(0.);
        Dune::BlockVector<Dune::FieldVector<Scalar,1>> elementSol2(0.);
        elementSol1.resize(unknowns1.size());
        elementSol2.resize(unknowns2.size());

        for (int idx=0; idx<numVertsOfElem1; ++idx)
        {
            for (int eqIdx1=0; eqIdx1<numEq1; ++eqIdx1)
                elementSol1[eqIdx1*numVertsOfElem1+idx] = unknowns1(lfsu1.child(eqIdx1),idx);
            for (int eqIdx2=0; eqIdx2<numEq2; ++eqIdx2)
                elementSol2[eqIdx2*numVertsOfElem2+idx] = unknowns2(lfsu2.child(eqIdx2),idx);
        }
#if HAVE_VALGRIND
        for (unsigned int i = 0; i < elementSol1.size(); i++)
            Valgrind::CheckDefined(elementSol1[i]);
        for (unsigned int i = 0; i < elementSol2.size(); i++)
            Valgrind::CheckDefined(elementSol2[i]);
#endif // HAVE_VALGRIND

        cParams.elemVolVarsPrev1.update(globalProblem_.sdProblem1(),
                                        sdElement1,
                                        cParams.fvGeometry1,
                                        true /* oldSol? */);
        cParams.elemVolVarsCur1.updatePDELab(globalProblem_.sdProblem1(),
                                             sdElement1,
                                             cParams.fvGeometry1,
                                             elementSol1);
        cParams.elemVolVarsPrev2.update(globalProblem_.sdProblem2(),
                                        sdElement2,
                                        cParams.fvGeometry2,
                                        true /* oldSol? */);
        cParams.elemVolVarsCur2.updatePDELab(globalProblem_.sdProblem2(),
                                             sdElement2,
                                             cParams.fvGeometry2,
                                             elementSol2);

        ElementBoundaryTypes1 bcTypes1;
        ElementBoundaryTypes2 bcTypes2;
        bcTypes1.update(globalProblem_.sdProblem1(), sdElement1, cParams.fvGeometry1);
        bcTypes2.update(globalProblem_.sdProblem2(), sdElement2, cParams.fvGeometry2);

        globalProblem_.localResidual1().evalPDELab(sdElement1, cParams.fvGeometry1,
                                                   cParams.elemVolVarsPrev1, cParams.elemVolVarsCur1,
                                                   bcTypes1);
        globalProblem_.localResidual2().evalPDELab(sdElement2, cParams.fvGeometry2,
                                                   cParams.elemVolVarsPrev2, cParams.elemVolVarsCur2,
                                                   bcTypes2);
    }

    /*!
     * \brief Evaluation of the coupling between the Stokes (1) and Darcy (2).
     *
     * Dirichlet-like and Neumann-like conditions for the respective domain are evaluated.
     *
     * \param lfsu1 local basis for the trial space of the Stokes domain
     * \param lfsu2 local basis for the trial space of the Darcy domain
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
    void evalCoupling(const LFSU1& lfsu1, const LFSU2& lfsu2,
                      const int vertInElem1, const int vertInElem2,
                      const SDElement1& sdElement1, const SDElement2& sdElement2,
                      const BoundaryVariables1& boundaryVars1, const BoundaryVariables2& boundaryVars2,
                      const CParams &cParams,
                      RES1& couplingRes1, RES2& couplingRes2) const
    {
        const GlobalPosition& globalPos1 = cParams.fvGeometry1.subContVol[vertInElem1].global;
        const GlobalPosition& globalPos2 = cParams.fvGeometry2.subContVol[vertInElem2].global;

        const GlobalPosition& bfNormal1 = boundaryVars1.face().normal;
        const Scalar normalMassFlux1 = boundaryVars1.normalVelocity()
                                       * cParams.elemVolVarsCur1[vertInElem1].density();

        // MASS Balance
        // Neumann-like conditions
        if (cParams.boundaryTypes1.isCouplingNeumann(massBalanceIdx1))
        {
            DUNE_THROW(Dune::NotImplemented, "The boundary condition isCouplingNeumann(massBalanceIdx1) for the Stokes side is not implemented.");
        }
        if (cParams.boundaryTypes2.isCouplingNeumann(massBalanceIdx2))
        {
            static_assert(!GET_PROP_VALUE(TwoPTwoCTypeTag, UseMoles),
                          "This coupling condition is only implemented for mass fraction formulation.");

            if (globalProblem_.sdProblem1().isCornerPoint(globalPos1))
            {
                couplingRes2.accumulate(lfsu2.child(massBalanceIdx2), vertInElem2,
                                        -normalMassFlux1);
            }
            else
            {
                couplingRes2.accumulate(lfsu2.child(massBalanceIdx2), vertInElem2,
                                        globalProblem_.localResidual1().residual(vertInElem1)[massBalanceIdx1]);
            }
        }

        // Dirichlet-like
        if (cParams.boundaryTypes1.isCouplingDirichlet(massBalanceIdx1))
        {
            DUNE_THROW(Dune::NotImplemented, "The boundary condition isCouplingDirichlet(massBalanceIdx1) for the Stokes side is not implemented.");
        }
        if (cParams.boundaryTypes2.isCouplingDirichlet(massBalanceIdx2))
        {
            couplingRes2.accumulate(lfsu2.child(massBalanceIdx2), vertInElem2,
                                    globalProblem_.localResidual1().residual(vertInElem1)[momentumYIdx1]
                                    -cParams.elemVolVarsCur1[vertInElem1].pressure());
        }


        // MOMENTUM_X Balance
        SpatialParams spatialParams = globalProblem_.sdProblem2().spatialParams();
        Scalar beaversJosephCoeff = spatialParams.beaversJosephCoeffAtPos(globalPos1);
        assert(beaversJosephCoeff > 0);
        beaversJosephCoeff /= std::sqrt(spatialParams.intrinsicPermeability(sdElement2, cParams.fvGeometry2, vertInElem2));

        // Neumann-like conditions
        if (cParams.boundaryTypes1.isCouplingNeumann(momentumXIdx1))
        {
            // v_tau = v - (v.n)n
            const Scalar normalComp = boundaryVars1.velocity()*bfNormal1;
            GlobalPosition normalV = bfNormal1;
            normalV *= normalComp;
            const GlobalPosition tangentialV = boundaryVars1.velocity() - normalV;

            // Implementation as Neumann-like condition: (v.n)n
            for (int dimIdx=0; dimIdx < dim; ++dimIdx)
            {
                couplingRes1.accumulate(lfsu1.child(momentumXIdx1), vertInElem1,
                                        beaversJosephCoeff
                                        * boundaryVars1.face().area
                                        * tangentialV[dimIdx]
                                        * (boundaryVars1.dynamicViscosity()
                                          + boundaryVars1.dynamicEddyViscosity()));
            }
        }

        // Dirichlet-like conditions
        if (cParams.boundaryTypes1.isCouplingDirichlet(momentumXIdx1))
        {
            // NOTE: This boundary condition is not implemented anymore because curPrimaryVars_ is protected
            DUNE_THROW(Dune::NotImplemented, "The boundary condition isCouplingNeumann(momentumXIdx1) on the Stokes side is not implemented anymore.");

            // tangential component: vx = sqrt K /alpha * (grad v n(unity))t
            // GlobalPosition tangentialVelGrad(0);
            // boundaryVars1.velocityGrad().umv(elementUnitNormal, tangentialVelGrad);
            // tangentialVelGrad /= -beaversJosephCoeff; // was - before
            // this->residual_[vertInElem1][momentumXIdx1] =
            //        tangentialVelGrad[momentumXIdx1] - globalProblem_.localResidual1().curPriVars_(vertInElem1)[momentumXIdx1]);
        }


        // MOMENTUM_Y Balance
        // Neumann-like conditions
        if (cParams.boundaryTypes1.isCouplingNeumann(momentumYIdx1))
        {
            // p*A as condition for free flow
            // pressure correction is done in stokeslocalresidual.hh
            couplingRes1.accumulate(lfsu1.child(momentumYIdx1), vertInElem1,
                                    cParams.elemVolVarsCur2[vertInElem2].pressure(nPhaseIdx2) *
                                    boundaryVars2.face().area);
        }

        // Dirichlet-like conditions
        if (cParams.boundaryTypes1.isCouplingDirichlet(momentumYIdx1))
        {
            // v.n as Dirichlet-like condition for the Stokes domain
            if (globalProblem_.sdProblem2().isCornerPoint(globalPos2))
            {
                Scalar sumNormalPhaseFluxes = 0.0;
                for (int phaseIdx=0; phaseIdx<numPhases2; ++phaseIdx)
                {
                    sumNormalPhaseFluxes -= boundaryVars2.volumeFlux(phaseIdx)
                                            * cParams.elemVolVarsCur2[vertInElem2].density(phaseIdx);
                }
                couplingRes1.accumulate(lfsu1.child(momentumYIdx1), vertInElem1,
                                        -sumNormalPhaseFluxes
                                        / cParams.elemVolVarsCur1[vertInElem1].density());
            }
            else
            {
                // set residualStokes[momentumYIdx1] = v_y in stokesnccouplinglocalresidual.hh
                couplingRes1.accumulate(lfsu1.child(momentumYIdx1), vertInElem1,
                                        globalProblem_.localResidual2().residual(vertInElem2)[massBalanceIdx2]
                                        / cParams.elemVolVarsCur1[vertInElem1].density());
            }
        }


        // COMPONENT Balance
        // Neumann-like conditions
        if (cParams.boundaryTypes1.isCouplingNeumann(transportEqIdx1))
        {
            DUNE_THROW(Dune::NotImplemented, "The boundary condition isCouplingNeumann(transportEqIdx1) is not implemented \
                                              for the Stokes side for multicomponent systems.");
        }
        if (cParams.boundaryTypes2.isCouplingNeumann(contiWEqIdx2))
        {
            // only enter here, if a boundary layer model is used for the computation of the diffusive fluxes
            if (blModel_)
            {
                Scalar advectiveFlux = normalMassFlux1
                                       * cParams.elemVolVarsCur1[vertInElem1].massFraction(transportCompIdx1);

                Scalar diffusiveFlux = bfNormal1.two_norm()
                                       * globalProblem_.evalBoundaryLayerConcentrationGradient(cParams, vertInElem1)
                                       * (boundaryVars1.diffusionCoeff(transportCompIdx1)
                                         + boundaryVars1.eddyDiffusivity())
                                       * boundaryVars1.molarDensity()
                                       * FluidSystem::molarMass(transportCompIdx1);

                const Scalar massTransferCoeff = globalProblem_.evalMassTransferCoefficient(cParams, vertInElem1, vertInElem2);

                if (massTransferModel_ && globalProblem_.sdProblem1().isCornerPoint(globalPos1))
                {
                    Scalar diffusiveFluxAtCorner = bfNormal1
                                                   * boundaryVars1.moleFractionGrad(transportCompIdx1)
                                                   * (boundaryVars1.diffusionCoeff(transportCompIdx1)
                                                      + boundaryVars1.eddyDiffusivity())
                                                   * boundaryVars1.molarDensity()
                                                   * FluidSystem::molarMass(transportCompIdx1);

                    couplingRes2.accumulate(lfsu2.child(contiWEqIdx2), vertInElem2,
                                            -massTransferCoeff*(advectiveFlux - diffusiveFlux) -
                                            (1.-massTransferCoeff)*(advectiveFlux - diffusiveFluxAtCorner));
                }
                else
                {
                    couplingRes2.accumulate(lfsu2.child(contiWEqIdx2), vertInElem2,
                                            -massTransferCoeff*(advectiveFlux - diffusiveFlux) +
                                            (1.-massTransferCoeff)*globalProblem_.localResidual1().residual(vertInElem1)[transportEqIdx1]);
                }
            }
            else if (globalProblem_.sdProblem1().isCornerPoint(globalPos1))
            {
                static_assert(!GET_PROP_VALUE(TwoPTwoCTypeTag, UseMoles),
                              "This coupling condition is only implemented for mass fraction formulation.");

                Scalar advectiveFlux = normalMassFlux1
                                       * cParams.elemVolVarsCur1[vertInElem1].massFraction(transportCompIdx1);

                Scalar diffusiveFlux = bfNormal1
                                       * boundaryVars1.moleFractionGrad(transportCompIdx1)
                                       * (boundaryVars1.diffusionCoeff(transportCompIdx1)
                                          + boundaryVars1.eddyDiffusivity())
                                       * boundaryVars1.molarDensity()
                                       * FluidSystem::molarMass(transportCompIdx1);

                couplingRes2.accumulate(lfsu2.child(contiWEqIdx2), vertInElem2,
                                        -(advectiveFlux - diffusiveFlux));
            }
            else
            {
                static_assert(GET_PROP_VALUE(Stokes2cTypeTag, UseMoles) == GET_PROP_VALUE(TwoPTwoCTypeTag, UseMoles),
                              "This coupling condition is not implemented for different formulations (mass/mole) in the subdomains.");

                // the component mass flux from the stokes domain
                couplingRes2.accumulate(lfsu2.child(contiWEqIdx2), vertInElem2,
                                        globalProblem_.localResidual1().residual(vertInElem1)[transportEqIdx1]);
            }
        }

        // Dirichlet-like conditions
        if (cParams.boundaryTypes1.isCouplingDirichlet(transportEqIdx1))
        {
            static_assert(!GET_PROP_VALUE(Stokes2cTypeTag, UseMoles),
                          "This coupling condition is only implemented for mass fraction formulation.");

            // set residualStokes[transportEqIdx1] = x in stokesnccouplinglocalresidual.hh
            // coupling residual is added to "real" residual
            couplingRes1.accumulate(lfsu1.child(transportEqIdx1), vertInElem1,
                                    -cParams.elemVolVarsCur2[vertInElem2].massFraction(nPhaseIdx2, wCompIdx2));
        }
        if (cParams.boundaryTypes2.isCouplingDirichlet(contiWEqIdx2))
        {
            DUNE_THROW(Dune::NotImplemented, "The boundary condition isCouplingDirichlet(contiWEqIdx2) is not implemented \
                                              for the Darcy side for multicomponent systems.");
        }
    }

 protected:
    GlobalProblem& globalProblem() const
    { return globalProblem_; }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

    unsigned int blModel_;
    unsigned int massTransferModel_;

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

    GlobalProblem& globalProblem_;
};

} // end namespace Dumux

#endif // DUMUX_2CSTOKES_2P2C_LOCALOPERATOR_HH
