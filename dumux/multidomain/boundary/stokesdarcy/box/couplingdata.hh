// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingData
 */

#ifndef DUMUX_STOKES_DARCY_BOX_COUPLINGDATA_HH
#define DUMUX_STOKES_DARCY_BOX_COUPLINGDATA_HH

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/math.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include <dumux/freeflow/navierstokes/staggered/velocitygradients.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <optional>

namespace Dumux {
/*!
 * \ingroup StokesDarcyCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataBoxBase : public StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;

    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FluidSystem = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id> using ModelTraits = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;

    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    using VelocityVector = typename Element<freeFlowIdx>::Geometry::GlobalCoordinate;
    template<std::size_t id> using BoundaryTypes = typename ProblemTraits<Problem<id>>::BoundaryTypes;
    using StokesVelocityGradients = StaggeredVelocityGradients<Scalar, GridGeometry<freeFlowIdx>, BoundaryTypes<freeFlowIdx>, Indices<freeFlowIdx>>;

    using AdvectionType = GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::AdvectionType>;
    using DarcysLaw = DarcysLawImplementation<SubDomainTypeTag<porousMediumIdx>, GridGeometry<porousMediumIdx>::discMethod>;
    using ForchheimersLaw = ForchheimersLawImplementation<SubDomainTypeTag<porousMediumIdx>, GridGeometry<porousMediumIdx>::discMethod>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    using ProjectionMethod = typename CouplingManager::ProjectionMethod;
    static constexpr auto projectionMethod = CouplingManager::projectionMethod;
public:
    StokesDarcyCouplingDataBoxBase(const CouplingManager& couplingmanager): ParentType(couplingmanager) {}

    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the momentum flux across the coupling boundary.
     *
     * Calculates the classical or new (normal)momentumCouplingCondition depending on the value of the parameter "Problem.NewIc"
     * Defaults to classical momentumCouplingCondition.
     *
     * For the normal momentum coupling, the porous medium side of the coupling condition
     * is evaluated, i.e. -[p n]^pm. For the new normal momentum coupling, the stokes term -viscosity*Nsbl*tau.T*grad(v_ff)*n is added.
     *
     */
    template<class ElementFaceVariables>
    Scalar momentumCouplingCondition(const Element<freeFlowIdx>& element,
                                     const FVElementGeometry<freeFlowIdx>& fvGeometry,
                                     const ElementVolumeVariables<freeFlowIdx>& stokesElemVolVars,
                                     const ElementFaceVariables& stokesElemFaceVars,
                                     const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        //################# pm contribution #################
        const auto darcyPhaseIdx = couplingPhaseIdx(porousMediumIdx);
        auto pressure = [darcyPhaseIdx](const auto& elemVolVars, const auto& scv)
                            { return elemVolVars[scv].pressure(darcyPhaseIdx); };

        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);
        Scalar momentumFlux = this->calculateProjection(scvf, stokesContext, pressure);

        // normalize pressure
        if(getPropValue<SubDomainTypeTag<freeFlowIdx>, Properties::NormalizePressure>())
            momentumFlux -= this->couplingManager().problem(freeFlowIdx).initial(scvf)[Indices<freeFlowIdx>::pressureIdx];

        //################# ff contribution #################
        static const bool newIc = getParamFromGroup<bool>("Problem", "NewIc", false);
        if (newIc)
        {
            static const bool unsymmetrizedGradientForBeaversJoseph = [&]()
            {
                const bool tmp = getParamFromGroup<bool>(this->couplingManager().problem(freeFlowIdx).paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);
                if (tmp)
                {
                    std::cerr << "Warning: You are using the deprecated parameter 'EnableUnsymmetrizedVelocityGradientForBeaversJoseph'. Use 'EnableUnsymmetrizedVelocityGradientForIC' instead."  << std::endl;
                }
                return tmp;
            }();

            // TODO: Replace unsymmetrizedGradientForBeaversJoseph below by false, when deprecation period expired
            static const bool unsymmetrizedGradientForIC = getParamFromGroup<bool>(this->couplingManager().problem(freeFlowIdx).paramGroup(),
                                                           "FreeFlow.EnableUnsymmetrizedVelocityGradientForIC", unsymmetrizedGradientForBeaversJoseph);
            // NewIc not verified for symmetrized gradient
            if(!unsymmetrizedGradientForIC)
                DUNE_THROW(Dune::NotImplemented, "Interface Conditions not verified for symmetrized stress tensors");

            const std::size_t numSubFaces = scvf.pairData().size();

            // Account for all sub faces
            for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
            {
                const auto eIdx = scvf.insideScvIdx();
                const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

                // Create a boundaryTypes object (will be empty if not at a boundary)
                std::optional<BoundaryTypes<freeFlowIdx>> currentScvfBoundaryTypes;
                if (scvf.boundary())
                {
                    currentScvfBoundaryTypes.emplace(this->couplingManager().problem(freeFlowIdx).boundaryTypes(element, scvf));
                }

                std::optional<BoundaryTypes<freeFlowIdx>> lateralFaceBoundaryTypes;
                if (lateralScvf.boundary())
                {
                    lateralFaceBoundaryTypes.emplace(this->couplingManager().problem(freeFlowIdx).boundaryTypes(element, lateralScvf));
                }

                // Get velocity gradients
                const Scalar velocityGrad_ji = StokesVelocityGradients::velocityGradJI(
                    this->couplingManager().problem(freeFlowIdx), element, fvGeometry, scvf , stokesElemFaceVars[scvf],
                    currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
                const Scalar velocityGrad_ij = unsymmetrizedGradientForIC ? 0.0 : StokesVelocityGradients::velocityGradIJ(
                    this->couplingManager().problem(freeFlowIdx), element, fvGeometry, scvf , stokesElemFaceVars[scvf],
                    currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);

                // Calculate stokes contribution to momentum flux: N_s^{bl} \tau T n
                const Scalar Nsbl = this->couplingManager().problem(freeFlowIdx).factorNMomentum(scvf);
                const Scalar viscosity = stokesElemVolVars[scvf.insideScvIdx()].viscosity();
                // Averaging the gradients over the subfaces to get evaluation at the center
                momentumFlux -= 1.0/numSubFaces * viscosity * Nsbl * (velocityGrad_ji + velocityGrad_ij);
            }
        }
        momentumFlux *= scvf.directionSign();
        return momentumFlux;
    }

    /*!
    * \brief Returns the averaged velocity vector at the interface of the porous medium according to darcys law
    *
    * The tangential porous medium velocity needs to be evaluated for the classical and new tangential coupling (slipCondition) at the
    * stokes-darcy interface, "Effective coupling conditions for arbitrary flows in Stokes-Darcy systems" by Elissa Eggenweiler.
    * We use darcys law and perform an integral average over all coupling segments.
    *
    * Depending on the parameter "Problem.NewIc" the standard permeability tensor K
    * or an altered permeability tensor M is used to evaluate the velocity. The method is using K per default.
    *
    */
    VelocityVector porousMediumVelocity(const Element<freeFlowIdx>& element, const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        static const bool newIc = getParamFromGroup<bool>("Problem", "NewIc", false);

        static constexpr int darcyDim = GridGeometry<porousMediumIdx>::GridView::dimension;
        using JacobianType = Dune::FieldMatrix<Scalar, 1, darcyDim>;
        std::vector<JacobianType> shapeDerivatives;
        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;

        VelocityVector velocity(0.0);       // velocity darcy                // density darcy
        Scalar intersectionLength = 0.0;    // (total)intersection length, could differ from scvf length

        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);
        static const bool enableGravity = getParamFromGroup<bool>(this->couplingManager().problem(porousMediumIdx).paramGroup(), "Problem.EnableGravity");

        // iteration over the different coupling segments
        for (const auto& data : stokesContext)
        {
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto darcyPhaseIdx = couplingPhaseIdx(porousMediumIdx);
                const auto& elemVolVars = *(data.elementVolVars);
                const auto& darcyFvGeometry = data.fvGeometry;
                const auto& darcyScvf = darcyFvGeometry.scvf(data.darcyScvfIdx);
                const auto& localBasis = darcyFvGeometry.feLocalBasis();

                // darcy Permeability
                const auto& K = data.permeability();

                // do second order integration as box provides linear functions
                const auto& rule = Dune::QuadratureRules<Scalar, darcyDim-1>::rule(data.segmentGeometry.type(), 2);
                for (const auto& qp : rule)
                {
                    const auto& ipLocal = qp.position();
                    const auto& ipGlobal = data.segmentGeometry.global(ipLocal);
                    const auto& ipElementLocal = data.element.geometry().local(ipGlobal);

                    VelocityVector gradP(0.0);
                    Scalar rho(0.0);

                    //calculate the shape and derivative values at the qp
                    localBasis.evaluateFunction(ipElementLocal, shapeValues);
                    localBasis.evaluateJacobian(ipElementLocal, shapeDerivatives);

                    //calc pressure gradient and rho at qp, every scv belongs to one node
                    for (const auto& scv : scvs(data.fvGeometry))
                    {
                        //gradP += p_i* (J^-T * L'_i)
                        data.element.geometry().jacobianInverseTransposed(ipElementLocal).usmv(elemVolVars[scv].pressure(darcyPhaseIdx), shapeDerivatives[scv.indexInElement()][0], gradP);
                        if (enableGravity)
                        {
                            rho += elemVolVars[scv].density(darcyPhaseIdx)*shapeValues[scv.indexInElement()][0];
                        }
                    }
                    //account for gravity
                    if (enableGravity)
                    {
                        gradP.axpy(-rho, this->couplingManager().problem(porousMediumIdx).spatialParams().gravity(ipGlobal));
                    }

                    if(newIc)
                    {
                        // darcy spatial dependent parameters
                        const auto& epsInterface = this->couplingManager().problem(freeFlowIdx).epsInterface(scvf);
                        const auto& M = this->couplingManager().problem(freeFlowIdx).matrixNTangential(scvf);
                        //Add the integrated segment velocity to the sum: v+= -w_k * sqrt(det(A^T*A))*eps**2*M/mu*gradP
                        velocity += mv(M, mv(qp.weight()*data.segmentGeometry.integrationElement(ipLocal)/elemVolVars[darcyScvf.insideScvIdx()].viscosity(darcyPhaseIdx)*epsInterface*epsInterface, gradP));
                    }
                    else
                    {
                        //add the integrated segment velocity to the sum: v+= -weight_k * sqrt(det(A^T*A))*K/mu*gradP
                        velocity += mv(K, mv(-qp.weight()*data.segmentGeometry.integrationElement(ipLocal)/elemVolVars[darcyScvf.insideScvIdx()].viscosity(darcyPhaseIdx), gradP));
                    }
                }
                intersectionLength += data.segmentGeometry.volume();
            }
        }
      velocity /= intersectionLength; //averaging
      return velocity;
    }

protected:
    // calculate projection needed for darcy residual
    template<class StokesContext, class Function>
    Scalar calculateProjection(const SubControlVolumeFace<freeFlowIdx>& stokesScvf,
                               const StokesContext& stokesContext,
                               const Element<porousMediumIdx>& darcyElement,
                               const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                               Function evalPriVar) const
    {
        Scalar projection = 0.0;

        // integrate darcy pressure over each coupling segment and average
        for (const auto& data : stokesContext)
        {
            //ToDo Is this if really necessary?
            if (stokesScvf.index() == data.stokesScvfIdx)
            {
                const auto darcyEIdxI = this->couplingManager().problem(porousMediumIdx).gridGeometry().elementMapper().index(darcyElement);
                const auto darcyEIdxJ = this->couplingManager().problem(porousMediumIdx).gridGeometry().elementMapper().index(data.element);
                const auto& elemVolVars = (darcyEIdxI == darcyEIdxJ) ? darcyElemVolVars : *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                projection += calculateSegmentIntegral(data.element, data.fvGeometry, darcyScvf, elemVolVars, data.segmentGeometry, evalPriVar);
            }
        }

        projection /= stokesScvf.area();

        return projection;
    }

    // calculate projection needed for stokes residual
    template<class StokesContext, class Function>
    Scalar calculateProjection(const SubControlVolumeFace<freeFlowIdx>& scvf,
                               const StokesContext& stokesContext,
                               Function evalPriVar) const
    {
        Scalar projection = 0.0;

        // integrate darcy pressure over each coupling segment and average
        for (const auto& data : stokesContext)
        {
            //ToDo Is this if really necessary?
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto& elemVolVars = *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                projection += calculateSegmentIntegral(data.element, data.fvGeometry, darcyScvf, elemVolVars, data.segmentGeometry, evalPriVar);
            }
        }

        projection /= scvf.area();

        return projection;
    }

    template<class SegementGeometry, class Function>
    Scalar calculateSegmentIntegral(const Element<porousMediumIdx>& element,
                                    const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                    const SubControlVolumeFace<porousMediumIdx>& scvf,
                                    const ElementVolumeVariables<porousMediumIdx>& elemVolVars,
                                    const SegementGeometry& segmentGeometry,
                                    Function evalPriVar) const
    {
        Scalar segmentProjection = 0.0;
        if constexpr (projectionMethod == ProjectionMethod::L2Projection)
        {
            const auto& localBasis = fvGeometry.feLocalBasis();

            // do second order integration as box provides linear functions
            static constexpr int darcyDim = GridGeometry<porousMediumIdx>::GridView::dimension;
            const auto& rule = Dune::QuadratureRules<Scalar, darcyDim-1>::rule(segmentGeometry.type(), 2);
            for (const auto& qp : rule)
            {
                const auto& ipLocal = qp.position();
                const auto& ipGlobal = segmentGeometry.global(ipLocal);
                const auto& ipElementLocal = element.geometry().local(ipGlobal);

                std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
                localBasis.evaluateFunction(ipElementLocal, shapeValues);

                Scalar value = 0.0;
                for (const auto& scv : scvs(fvGeometry))
                    value += evalPriVar(elemVolVars, scv)*shapeValues[scv.indexInElement()][0];

                segmentProjection += value*segmentGeometry.integrationElement(qp.position())*qp.weight();
            }
        }
        else if constexpr (projectionMethod == ProjectionMethod::AreaWeightedDofEvaluation)
        {
            segmentProjection = segmentGeometry.volume()*evalPriVar(elemVolVars, fvGeometry.scv(scvf.insideScvIdx()));
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Unkown projection method!");
        }

        return segmentProjection;
    }
};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, false, DiscretizationMethod::box>
: public StokesDarcyCouplingDataBoxBase<MDTraits, CouplingManager, enableEnergyBalance>
{
    using ParentType = StokesDarcyCouplingDataBoxBase<MDTraits, CouplingManager, enableEnergyBalance>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

    static_assert(GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    using ProjectionMethod = typename CouplingManager::ProjectionMethod;
    static constexpr auto projectionMethod = CouplingManager::projectionMethod;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    Scalar massCouplingCondition(const Element<porousMediumIdx>& element,
                                 const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                 const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                                 const ElementFluxVariablesCache<porousMediumIdx>& elementFluxVarsCache,
                                 const SubControlVolumeFace<porousMediumIdx>& scvf) const
    {
        const auto darcyPhaseIdx = couplingPhaseIdx(porousMediumIdx);
        const auto& darcyContext = this->couplingManager().darcyCouplingContextVector(element, scvf);

        Scalar flux = 0.0;
        for(const auto& data : darcyContext)
        {
            if(scvf.index() == data.darcyScvfIdx)
            {
                const Scalar velocity = data.velocity * scvf.unitOuterNormal();

                Scalar darcyDensity = darcyElemVolVars[scvf.insideScvIdx()].density(darcyPhaseIdx);

                const Scalar stokesDensity = data.volVars.density();
                const bool insideIsUpstream = velocity > 0.0;
                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += massFlux_(velocity, darcyDensity, stokesDensity, insideIsUpstream)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    Scalar massCouplingCondition(const Element<freeFlowIdx>& element,
                                 const FVElementGeometry<freeFlowIdx>& fvGeometry,
                                 const ElementVolumeVariables<freeFlowIdx>& stokesElemVolVars,
                                 const ElementFaceVariables<freeFlowIdx>& stokesElemFaceVars,
                                 const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);

        Scalar flux = 0.0;
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf() * scvf.directionSign();
        const Scalar stokesDensity = stokesElemVolVars[scvf.insideScvIdx()].density();
        for (const auto& data : stokesContext)
        {
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto darcyPhaseIdx = couplingPhaseIdx(porousMediumIdx);
                const auto& elemVolVars = *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                const auto darcyDensity = elemVolVars[darcyScvf.insideScvIdx()].density(darcyPhaseIdx);

                const bool insideIsUpstream = velocity > 0.0;
                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += massFlux_(velocity, stokesDensity, darcyDensity, insideIsUpstream)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<porousMediumIdx>& element,
                                   const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                   const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                                   const ElementFluxVariablesCache<porousMediumIdx>& elementFluxVarsCache,
                                   const SubControlVolumeFace<porousMediumIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContextVector(element, scvf);

        Scalar flux = 0.0;
        for(const auto& data : darcyContext)
        {
            if(scvf.index() == data.darcyScvfIdx)
            {
                //calculate the free-flow velocity
                const Scalar velocity = -1*(data.velocity * scvf.unitOuterNormal());

                const auto& stokesScvf = data.fvGeometry.scvf(data.stokesScvfIdx);
                const auto& stokesVolVars = data.volVars;

                //Calculate the projected temperature value for the stokes face
                auto temp = [](const auto& elemVolVars, const auto& scv)
                                { return elemVolVars[scv].temperature(); };
                Scalar interfaceTemperature = this->calculateProjection(stokesScvf, *data.stokesContext, element, darcyElemVolVars, temp);

                const bool insideIsUpstream = velocity > 0.0;
                flux += -1*energyFlux_(data.fvGeometry,
                                       stokesVolVars,
                                       stokesScvf,
                                       darcyElemVolVars[scvf.insideScvIdx()],
                                       velocity,
                                       interfaceTemperature,
                                       insideIsUpstream,
                                       diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<freeFlowIdx>& element,
                                   const FVElementGeometry<freeFlowIdx>& fvGeometry,
                                   const ElementVolumeVariables<freeFlowIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<freeFlowIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<freeFlowIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);

        //Calculate the projected temperature value for the stokes face
        auto temp = [](const auto& elemVolVars, const auto& scv)
                    { return elemVolVars[scv].temperature(); };

        Scalar interfaceTemperature = this->calculateProjection(scvf, stokesContext, temp);

        Scalar flux = 0.0;
        for (const auto& data : stokesContext)
        {
            if (scvf.index() == data.stokesScvfIdx)
            {
                const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf() * scvf.directionSign();

                const auto& elemVolVars = *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                const bool insideIsUpstream = velocity > 0.0;
                flux += energyFlux_(fvGeometry,
                                    stokesElemVolVars[scvf.insideScvIdx()],
                                    scvf,
                                    elemVolVars[darcyScvf.insideScvIdx()],
                                    velocity,
                                    interfaceTemperature,
                                    insideIsUpstream,
                                    diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

private:

    /*!
     * \brief Evaluate the mole/mass flux across the interface.
     */
    Scalar massFlux_(const Scalar velocity,
                     const Scalar insideDensity,
                     const Scalar outSideDensity,
                     bool insideIsUpstream) const
    {
        return this->advectiveFlux(insideDensity, outSideDensity, velocity, insideIsUpstream);
    }

    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(const FVElementGeometry<freeFlowIdx>& stokesFvGeometry,
                       const VolumeVariables<freeFlowIdx>& stokesVolVars,
                       const SubControlVolumeFace<freeFlowIdx>& scvfStokes,
                       const VolumeVariables<porousMediumIdx>& darcyVolVars,
                       const Scalar velocity,
                       const Scalar interfaceTemperature,
                       const bool insideIsUpstream,
                       const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const Scalar stokesTerm = stokesVolVars.density(couplingPhaseIdx(freeFlowIdx)) * stokesVolVars.enthalpy(couplingPhaseIdx(freeFlowIdx));
        const Scalar darcyTerm = darcyVolVars.density(couplingPhaseIdx(porousMediumIdx)) * darcyVolVars.enthalpy(couplingPhaseIdx(porousMediumIdx));

        flux += this->advectiveFlux(stokesTerm, darcyTerm, velocity, insideIsUpstream);

        const auto& stokesScv = (*scvs(stokesFvGeometry).begin());

        const Scalar deltaT = interfaceTemperature - stokesVolVars.temperature();
        const Scalar dist = (scvfStokes.center() - stokesScv.center()).two_norm();
        if(diffCoeffAvgType == DiffusionCoefficientAveragingType::ffOnly)
        {
            flux += -1*stokesVolVars.effectiveThermalConductivity() * deltaT / dist ;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Multidomain staggered box coupling only works for DiffusionCoefficientAveragingType = ffOnly");

        return flux;
    }
};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, true, DiscretizationMethod::box>
: public StokesDarcyCouplingDataBoxBase<MDTraits, CouplingManager, enableEnergyBalance>
{
    using ParentType = StokesDarcyCouplingDataBoxBase<MDTraits, CouplingManager, enableEnergyBalance>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using FluidSystem  = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;

    static constexpr auto numComponents = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::ModelTraits>::numFluidComponents();
    static constexpr auto replaceCompEqIdx = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::ModelTraits>::replaceCompEqIdx();
    static constexpr bool useMoles = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::ModelTraits>::useMoles();
    static constexpr auto referenceSystemFormulation = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation();

    static_assert(GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::ModelTraits>::numFluidComponents() == numComponents, "Both submodels must use the same number of components");
    static_assert(getPropValue<SubDomainTypeTag<porousMediumIdx>, Properties::UseMoles>() == useMoles, "Both submodels must either use moles or not");
    static_assert(getPropValue<SubDomainTypeTag<porousMediumIdx>, Properties::ReplaceCompEqIdx>() == replaceCompEqIdx, "Both submodels must use the same replaceCompEqIdx");
    static_assert(GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation() == referenceSystemFormulation,
                  "Both submodels must use the same reference system formulation for diffusion");

    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr bool isFicksLaw = IsFicksLaw<GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::MolecularDiffusionType>>();
    static_assert(isFicksLaw == IsFicksLaw<GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::MolecularDiffusionType>>(),
                  "Both submodels must use the same diffusion law.");

    static_assert(isFicksLaw, "Box-Staggered Coupling only implemented for Fick's law!");

    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

    using MolecularDiffusionType = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::MolecularDiffusionType>;

    using ProjectionMethod = typename CouplingManager::ProjectionMethod;
    static constexpr auto projectionMethod = CouplingManager::projectionMethod;
public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
    using ParentType::couplingCompIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    NumEqVector massCouplingCondition(const Element<porousMediumIdx>& element,
                                      const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                      const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                                      const ElementFluxVariablesCache<porousMediumIdx>& elementFluxVarsCache,
                                      const SubControlVolumeFace<porousMediumIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContextVector(element, scvf);

        NumEqVector flux(0.0);
        for(const auto& data : darcyContext)
        {
            if(scvf.index() == data.darcyScvfIdx)
            {
                // get the free-flow velocity
                const Scalar velocity = -1*(data.velocity * scvf.unitOuterNormal());
                const bool insideIsUpstream = velocity > 0.0;

                const auto& stokesScvf = data.fvGeometry.scvf(data.stokesScvfIdx);
                const auto& stokesVolVars = data.volVars;

                // Calculate the projected massOrMoleFraction value for the stokes face
                const auto& stokesContext = *data.stokesContext;
                auto interfaceMassOrMoleFraction = [this, &stokesScvf, &stokesContext, &element, &darcyElemVolVars](int compIdx)
                {
                    auto value = [&compIdx](const auto& elemVolVars, const auto& scv)
                                    { return massOrMoleFraction(elemVolVars[scv], referenceSystemFormulation, couplingPhaseIdx(porousMediumIdx), compIdx); };

                    return this->calculateProjection(stokesScvf, stokesContext, element, darcyElemVolVars, value);
                };

                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += -1*massFlux_(porousMediumIdx,
                                     freeFlowIdx,
                                     data.fvGeometry,
                                     fvGeometry,
                                     stokesVolVars,
                                     darcyElemVolVars,
                                     stokesScvf,
                                     scvf,
                                     velocity,
                                     interfaceMassOrMoleFraction,
                                     insideIsUpstream,
                                     diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    NumEqVector massCouplingCondition(const Element<freeFlowIdx>& element,
                                      const FVElementGeometry<freeFlowIdx>& fvGeometry,
                                      const ElementVolumeVariables<freeFlowIdx>& stokesElemVolVars,
                                      const ElementFaceVariables<freeFlowIdx>& stokesElemFaceVars,
                                      const SubControlVolumeFace<freeFlowIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);

        //Calculate the projected massOrMoleFraction value for the stokes face
        auto interfaceMassOrMoleFraction = [this, &scvf, &stokesContext](int compIdx)
        {
            auto value = [&compIdx](const auto& elemVolVars, const auto& scv)
                            { return massOrMoleFraction(elemVolVars[scv], referenceSystemFormulation, couplingPhaseIdx(porousMediumIdx), compIdx); };

            return this->calculateProjection(scvf, stokesContext, value);
        };

        NumEqVector flux(0.0);
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf() * scvf.directionSign();
        for (const auto& data : stokesContext)
        {
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto& darcyElemVolVars = *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                const bool insideIsUpstream = velocity > 0.0;
                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += massFlux_(freeFlowIdx,
                                  porousMediumIdx,
                                  fvGeometry,
                                  data.fvGeometry,
                                  stokesElemVolVars[scvf.insideScvIdx()],
                                  darcyElemVolVars,
                                  scvf,
                                  darcyScvf,
                                  velocity,
                                  interfaceMassOrMoleFraction,
                                  insideIsUpstream,
                                  diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<porousMediumIdx>& element,
                                   const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                   const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                                   const ElementFluxVariablesCache<porousMediumIdx>& elementFluxVarsCache,
                                   const SubControlVolumeFace<porousMediumIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContextVector(element, scvf);

        Scalar flux = 0.0;
        for(const auto& data : darcyContext)
        {
            if(scvf.index() == data.darcyScvfIdx)
            {
                // get the free-flow velocity
                const Scalar velocity = -1*(data.velocity * scvf.unitOuterNormal());
                const bool insideIsUpstream = velocity > 0.0;

                const auto& stokesScvf = data.fvGeometry.scvf(data.stokesScvfIdx);
                const auto& stokesVolVars = data.volVars;

                //Calculate the projected massOrMoleFraction value for the stokes face
                const auto& stokesContext = *data.stokesContext;
                auto interfaceMassOrMoleFraction = [this, &stokesScvf, &stokesContext, &element, &darcyElemVolVars](int compIdx)
                {
                    auto value = [&compIdx](const auto& elemVolVars, const auto& scv)
                                    { return massOrMoleFraction(elemVolVars[scv], referenceSystemFormulation, couplingPhaseIdx(porousMediumIdx), compIdx); };

                    return this->calculateProjection(stokesScvf, stokesContext, element, darcyElemVolVars, value);
                };

                //Calculate the projected temperature value for the stokes face
                auto temp = [](const auto& elemVolVars, const auto& scv)
                                { return elemVolVars[scv].temperature(); };
                Scalar interfaceTemperature = this->calculateProjection(stokesScvf, stokesContext, element, darcyElemVolVars, temp);

                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += -1*energyFlux_(porousMediumIdx,
                                       freeFlowIdx,
                                       data.fvGeometry,
                                       fvGeometry,
                                       stokesVolVars,
                                       darcyElemVolVars,
                                       stokesScvf,
                                       scvf,
                                       velocity,
                                       interfaceMassOrMoleFraction,
                                       interfaceTemperature,
                                       insideIsUpstream,
                                       diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<freeFlowIdx>& element,
                                   const FVElementGeometry<freeFlowIdx>& fvGeometry,
                                   const ElementVolumeVariables<freeFlowIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<freeFlowIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<freeFlowIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContextVector(element, scvf);

        //Calculate the projected massOrMoleFraction value for the stokes face
        auto interfaceMassOrMoleFraction = [this, &scvf, &stokesContext](int compIdx)
        {
            auto value = [&compIdx](const auto& elemVolVars, const auto& scv)
                            { return massOrMoleFraction(elemVolVars[scv], referenceSystemFormulation, couplingPhaseIdx(porousMediumIdx), compIdx); };

            return this->calculateProjection(scvf, stokesContext, value);
        };

        //Calculate the projected temperature value for the stokes face
        auto temp = [](const auto& elemVolVars, const auto& scv)
                    { return elemVolVars[scv].temperature(); };

        Scalar interfaceTemperature = this->calculateProjection(scvf, stokesContext, temp);

        Scalar flux = 0.0;
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf() * scvf.directionSign();
        for (const auto& data : stokesContext)
        {
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto& darcyElemVolVars = *(data.elementVolVars);
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);

                const bool insideIsUpstream = velocity > 0.0;
                // Division by scvf.area() is needed, because the final flux results from multiplication with scvf.area()
                flux += energyFlux_(freeFlowIdx,
                                    porousMediumIdx,
                                    fvGeometry,
                                    data.fvGeometry,
                                    stokesElemVolVars[scvf.insideScvIdx()],
                                    darcyElemVolVars,
                                    scvf,
                                    darcyScvf,
                                    velocity,
                                    interfaceMassOrMoleFraction,
                                    interfaceTemperature,
                                    insideIsUpstream,
                                    diffCoeffAvgType)*data.segmentGeometry.volume()/scvf.area();
            }
        }

        return flux;
    }

protected:

    /*!
     * \brief Evaluate the compositional mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j, class Function>
    NumEqVector massFlux_(Dune::index_constant<i> domainI,
                          Dune::index_constant<j> domainJ,
                          const FVElementGeometry<freeFlowIdx>& stokesFvGeometry,
                          const FVElementGeometry<porousMediumIdx>& darcyFvGeometry,
                          const VolumeVariables<freeFlowIdx>& stokesVolVars,
                          const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                          const SubControlVolumeFace<freeFlowIdx>& stokesScvf,
                          const SubControlVolumeFace<porousMediumIdx>& darcyScvf,
                          const Scalar velocity,
                          const Function& interfaceMoleOrMassFraction,
                          const bool insideIsUpstream,
                          const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector flux(0.0);
        NumEqVector diffusiveFlux(0.0);

        const auto& darcyVolVars = darcyElemVolVars[darcyScvf.insideScvIdx()];

        auto moleOrMassFraction = [](const auto& volVars, int phaseIdx, int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto moleOrMassDensity = [](const auto& volVars, int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // treat the advective fluxes
        auto insideTerm = [&](int compIdx)
        { return moleOrMassFraction(stokesVolVars, couplingPhaseIdx(freeFlowIdx), compIdx) * moleOrMassDensity(stokesVolVars, couplingPhaseIdx(freeFlowIdx)); };

        // ToDO interpolate using box basis functions
        auto outsideTerm = [&](int compIdx)
        { return moleOrMassFraction(darcyVolVars, couplingPhaseIdx(porousMediumIdx), compIdx) * moleOrMassDensity(darcyVolVars, couplingPhaseIdx(porousMediumIdx)); };

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(freeFlowIdx, compIdx);
            const int domainJCompIdx = couplingCompIdx(porousMediumIdx, compIdx);
            flux[couplingCompIdx(domainI, compIdx)] += this->advectiveFlux(insideTerm(domainICompIdx), outsideTerm(domainJCompIdx), velocity, insideIsUpstream);
        }

        // treat the diffusive fluxes
        diffusiveFlux += diffusiveMolecularFluxFicksLaw_(domainI,
                                                         domainJ,
                                                         stokesFvGeometry,
                                                         darcyFvGeometry,
                                                         stokesVolVars,
                                                         darcyElemVolVars,
                                                         stokesScvf,
                                                         darcyScvf,
                                                         velocity,
                                                         interfaceMoleOrMassFraction,
                                                         diffCoeffAvgType);

        //convert to correct units if necessary
        if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged && useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] *= 1/FluidSystem<domainI>::molarMass(domainICompIdx);
            }
        }
        if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged && !useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] *= FluidSystem<domainI>::molarMass(domainICompIdx);
            }
        }

        flux += diffusiveFlux;
        // convert to total mass/mole balance, if set be user
        if (replaceCompEqIdx < numComponents)
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);

        return flux;
    }

    Scalar getComponentEnthalpy(const VolumeVariables<freeFlowIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<freeFlowIdx>::componentEnthalpy(volVars.fluidState(), 0, compIdx);
    }

    Scalar getComponentEnthalpy(const VolumeVariables<porousMediumIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<porousMediumIdx>::componentEnthalpy(volVars.fluidState(), phaseIdx, compIdx);
    }

    template<std::size_t i, std::size_t j, class Function>
    NumEqVector diffusiveMolecularFluxFicksLaw_(Dune::index_constant<i> domainI,
                                                Dune::index_constant<j> domainJ,
                                                const FVElementGeometry<freeFlowIdx>& stokesFvGeometry,
                                                const FVElementGeometry<porousMediumIdx>& darcyFvGeometry,
                                                const VolumeVariables<freeFlowIdx>& stokesVolVars,
                                                const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                                                const SubControlVolumeFace<freeFlowIdx>& stokesScvf,
                                                const SubControlVolumeFace<porousMediumIdx>& darcyScvf,
                                                const Scalar velocity,
                                                const Function& interfaceMoleOrMassFraction,
                                                const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector diffusiveFlux(0.0);

        const Scalar rhoStokes = massOrMolarDensity(stokesVolVars, referenceSystemFormulation, couplingPhaseIdx(freeFlowIdx));
        const Scalar rhoDarcy = massOrMolarDensity(darcyElemVolVars[darcyScvf.insideScvIdx()], referenceSystemFormulation, couplingPhaseIdx(porousMediumIdx));
        const Scalar avgDensity = 0.5 * (rhoStokes + rhoDarcy);

        for (int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            const int stokesMainCompIdx = couplingPhaseIdx(freeFlowIdx);
            const int stokesCompIdx = couplingCompIdx(freeFlowIdx, compIdx);
            const int darcyCompIdx = couplingCompIdx(porousMediumIdx, compIdx);

            assert(FluidSystem<freeFlowIdx>::componentName(stokesCompIdx) == FluidSystem<porousMediumIdx>::componentName(darcyCompIdx));

            const Scalar massOrMoleFractionStokes = massOrMoleFraction(stokesVolVars, referenceSystemFormulation, couplingPhaseIdx(freeFlowIdx), stokesCompIdx);

            const Scalar deltaMassOrMoleFrac = interfaceMoleOrMassFraction(darcyCompIdx) - massOrMoleFractionStokes;
            const auto& stokesScv = (*scvs(stokesFvGeometry).begin());
            const Scalar dist = (stokesScvf.center() - stokesScv.center()).two_norm();
            if(diffCoeffAvgType == DiffusionCoefficientAveragingType::ffOnly)
                diffusiveFlux[couplingCompIdx(domainI, compIdx)] += -avgDensity * stokesVolVars.effectiveDiffusionCoefficient(couplingPhaseIdx(freeFlowIdx), stokesMainCompIdx, stokesCompIdx)
                                                                     * deltaMassOrMoleFrac / dist;
            else
                DUNE_THROW(Dune::NotImplemented, "Multidomain staggered box coupling only works for DiffusionCoefficientAveragingType = ffOnly");
        }

        const Scalar cumulativeFlux = std::accumulate(diffusiveFlux.begin(), diffusiveFlux.end(), 0.0);
        diffusiveFlux[couplingCompIdx(domainI, 0)] = -cumulativeFlux;

        return diffusiveFlux;
    }

    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, class Function, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(Dune::index_constant<i> domainI,
                       Dune::index_constant<j> domainJ,
                       const FVElementGeometry<freeFlowIdx>& stokesFvGeometry,
                       const FVElementGeometry<porousMediumIdx>& darcyFvGeometry,
                       const VolumeVariables<freeFlowIdx>& stokesVolVars,
                       const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                       const SubControlVolumeFace<freeFlowIdx>& stokesScvf,
                       const SubControlVolumeFace<porousMediumIdx>& darcyScvf,
                       const Scalar velocity,
                       const Function& interfaceMoleOrMassFraction,
                       const Scalar interfaceTemperature,
                       const bool insideIsUpstream,
                       const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const auto& stokesScv = (*scvs(stokesFvGeometry).begin());
        const auto& darcyVolVars = darcyElemVolVars[darcyScvf.insideScvIdx()];

        const Scalar stokesTerm = stokesVolVars.density(couplingPhaseIdx(freeFlowIdx)) * stokesVolVars.enthalpy(couplingPhaseIdx(freeFlowIdx));
        // ToDO interpolate using box basis functions
        const Scalar darcyTerm = darcyVolVars.density(couplingPhaseIdx(porousMediumIdx)) * darcyVolVars.enthalpy(couplingPhaseIdx(porousMediumIdx));

        flux += this->advectiveFlux(stokesTerm, darcyTerm, velocity, insideIsUpstream);

        const Scalar deltaT = interfaceTemperature - stokesVolVars.temperature();
        const Scalar dist = (stokesScvf.center() - stokesScv.center()).two_norm();
        if(diffCoeffAvgType == DiffusionCoefficientAveragingType::ffOnly)
        {
            flux += -1*stokesVolVars.effectiveThermalConductivity() * deltaT / dist;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Multidomain staggered box coupling only works for DiffusionCoefficientAveragingType = ffOnly");

        auto diffusiveFlux = diffusiveMolecularFluxFicksLaw_(domainI,
                                                             domainJ,
                                                             stokesFvGeometry,
                                                             darcyFvGeometry,
                                                             stokesVolVars,
                                                             darcyElemVolVars,
                                                             stokesScvf,
                                                             darcyScvf,
                                                             velocity,
                                                             interfaceMoleOrMassFraction,
                                                             diffCoeffAvgType);

        for (int compIdx = 0; compIdx < diffusiveFlux.size(); ++compIdx)
        {
            const int stokesCompIdx = couplingCompIdx(freeFlowIdx, compIdx);
            const int darcyCompIdx = couplingCompIdx(porousMediumIdx, compIdx);
            const int domainCompIdx = couplingCompIdx(domainI, compIdx);

            const Scalar componentEnthalpy = diffusiveFlux[domainCompIdx] > 0 ?
                                             getComponentEnthalpy(stokesVolVars, couplingPhaseIdx(freeFlowIdx), stokesCompIdx)
                                           : getComponentEnthalpy(darcyVolVars, couplingPhaseIdx(porousMediumIdx), darcyCompIdx);

            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                flux += diffusiveFlux[domainCompIdx] * componentEnthalpy;
            else
                flux += diffusiveFlux[domainCompIdx] * FluidSystem<domainI>::molarMass(domainCompIdx) * componentEnthalpy;
        }

        return flux;
    }
};

} // end namespace Dumux

#endif // DUMUX_STOKES_DARCY_COUPLINGDATA_HH
