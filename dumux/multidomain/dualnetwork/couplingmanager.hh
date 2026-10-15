// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DualNetworkCoupling
 * \ingroup DarcyDarcyCoupling
 * \brief Coupling manager for equal-dimension boundary coupling
 */

#ifndef DUMUX_MULTIDOMAIN_DUAL_NETWORK_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_DUAL_NETWORK_COUPLINGMANAGER_HH

#include <algorithm>
#include <array>
#include <deque>
#include <iostream>
#include <type_traits>
#include <vector>
#include <tuple>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/common/enumerate.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/dimensionlessnumbers.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/format.hh>
#include <dumux/assembly/coloring.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/assemblytraits.hh>
#include "extendedsourcestencil.hh"

#include "couplingmapper.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup DualNetworkCoupling
 * \ingroup DarcyDarcyCoupling
 * \brief Coupling manager for dual network approach for pore network models
 * \note Concept and algorithms described in Koch et al (2021) https://doi.org/10.1007/s11242-021-01602-5
 */
template<class MDTraits>
class PNMHeatTransferCouplingManager
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using ThisType = PNMHeatTransferCouplingManager<MDTraits>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using Problem = GetPropType<SubDomainTypeTag<i>, Properties::Problem>;
    template<std::size_t i> using PrimaryVariables = GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>;
    template<std::size_t i> using NumEqVector = Dumux::NumEqVector<PrimaryVariables<i>>;
    template<std::size_t i> using GridVolumeVariables = GetPropType<SubDomainTypeTag<i>, Properties::GridVolumeVariables>;
    template<std::size_t i> using ElementVolumeVariables = typename GridVolumeVariables<i>::LocalView;
    template<std::size_t i> using ElementFluxVariablesCache = typename GetPropType<SubDomainTypeTag<i>, Properties::GridFluxVariablesCache>::LocalView;
    template<std::size_t i> using VolumeVariables = typename GridVolumeVariables<i>::VolumeVariables;
    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using FVElementGeometry = typename GridGeometry<i>::LocalView;
    template<std::size_t i> using SubControlVolumeFace = typename GridGeometry<i>::SubControlVolumeFace;
    template<std::size_t i> using SubControlVolume = typename GridGeometry<i>::SubControlVolume;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;
    template<std::size_t i> using Element = typename GridView<i>::template Codim<0>::Entity;
    template<std::size_t i> using ElementSeed = typename GridView<i>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t i> using Indices = typename GetPropType<SubDomainTypeTag<i>, Properties::ModelTraits>::Indices;

    template<std::size_t id> using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    using CouplingMapper = DualNetworkCouplingMapper<Scalar>;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    using CouplingStencil = std::vector<std::size_t>;
    using DimLessNum = DimensionlessNumbers<Scalar>;


public:
    static constexpr auto solidDomainIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto voidDomainIdx = typename MDTraits::template SubDomain<1>::Index();

private:

    //! whether the void subproblem keeps its volume variables / flux variables caches grid-wide
    static constexpr bool voidVolVarsCachingEnabled
        = getPropValue<SubDomainTypeTag<voidDomainIdx>, Properties::EnableGridVolumeVariablesCache>();
    static constexpr bool voidFluxVarsCachingEnabled
        = getPropValue<SubDomainTypeTag<voidDomainIdx>, Properties::EnableGridFluxVariablesCache>();

    struct CouplingContextForOneConnection
    {
        std::size_t connectionGlobalId;
        std::vector<FVElementGeometry<solidDomainIdx>> solidFVGeometry; // only one needed, but FVElementGeometry can't be default constructed
        VolumeVariables<solidDomainIdx> solidVolVars;
        std::vector<FVElementGeometry<voidDomainIdx>> voidFVGeometry;
        std::vector<ElementVolumeVariables<voidDomainIdx>> voidElemVolVars;
        std::vector<ElementFluxVariablesCache<voidDomainIdx>> voidElemFluxVarsCache;
    };

    using CouplingContextForOnePore = std::pair<std::size_t, std::vector<CouplingContextForOneConnection>>;

    /*!
     * \brief The coupling context of one element, holding for each of its (at most two) coupled
     *        pores the data needed to evaluate that pore's convective coupling source
     * \note Building the context is far too expensive to redo for every deflected dof, so it is
     *       kept while the derivatives of one element are computed and only patched up by
     *       updateCouplingContext(). One instance exists per thread, see elementCouplingContext_().
     */
    struct ElementCouplingContext
    {
        auto& data()
        { return data_; }

        template<std::size_t i>
        void bind(Dune::index_constant<i> domainI,
                  const Element<i>& element,
                  const ThisType& couplingManager)
        {
            // do nothing if the context is already bound to this element by this manager
            const auto eIdx = couplingManager.gridGeometry(domainI).elementMapper().index(element);
            if (boundManager_ == &couplingManager && boundDomainId_ == i && boundElementIndex_ == eIdx)
                return;

            boundManager_ = &couplingManager;
            boundDomainId_ = i;
            boundElementIndex_ = eIdx;
            data_.clear();

            // nothing to bind if the element does not couple at all
            // (this check is more expensive, so we do it after the above one)
            constexpr auto domainJ = Dune::index_constant<1-i>{};
            if (couplingManager.couplingStencil(domainI, element, domainJ).empty())
                return;

            // each element has two pores for which we may need to bind a coupling context
            for (int localVertexIdx = 0; localVertexIdx < 2; ++localVertexIdx)
            {
                const auto vIdx = couplingManager.gridGeometry(domainI).vertexMapper().subIndex(element, localVertexIdx, 1);
                if (!couplingManager.isCoupledPore(domainI, vIdx))
                    continue;

                data_.push_back(CouplingContextForOnePore{
                    vIdx, couplingManager.couplingContextForPore(domainI, vIdx)
                });
            }
        }

        const auto& operator[](const std::size_t dofIdx) const
        {
            for (const auto& d : data_)
            {
                if (d.first == dofIdx)
                    return d.second;
            }
            DUNE_THROW(Dune::InvalidStateException, "No connection found");
        }

    private:
        const ThisType* boundManager_ = nullptr;
        std::size_t boundDomainId_ = 0;
        std::size_t boundElementIndex_ = 0;
        Dune::ReservedVector<CouplingContextForOnePore, 2> data_;
    };

    /*!
     * \brief The coupling context of the element that the calling thread is assembling
     * \note Each thread has its own, which is what makes multithreaded assembly possible. The
     *       context also remembers the manager that bound it, so that two managers of the same
     *       type can never see each other's context.
     */
    ElementCouplingContext& elementCouplingContext_() const
    {
        static thread_local ElementCouplingContext context;
        return context;
    }

public:

    //! export traits
    using MultiDomainTraits = MDTraits;
    //! export stencil types
    using CouplingStencils = std::unordered_map<std::size_t, CouplingStencil>;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    template<class HostGridView, class HostGridData, class VoidGridView, class SolidGridView>
    void init(std::shared_ptr<Problem<solidDomainIdx>> solidProblem,
              std::shared_ptr<Problem<voidDomainIdx>> voidProblem,
              const HostGridView& hostGridView,
              const HostGridData& hostGridData,
              const VoidGridView& voidGridView,
              const SolidGridView& solidGridView,
              const SolutionVector& curSol)
    {
        voidGridGeometry_ = &voidProblem->gridGeometry();
        solidGridGeometry_ = &solidProblem->gridGeometry();
        couplingMapper_ = std::make_unique<CouplingMapper>(hostGridView, hostGridData, voidProblem->gridGeometry(), solidProblem->gridGeometry());
        this->setSubProblems(std::make_tuple(solidProblem, voidProblem));
        this->updateSolution(curSol);

        setupExtendedStencil();
    }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param element the coupled element of domain í
     * \param domainJ the domain index of domain j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto eIdx = gridGeometry(domainI).elementMapper().index(element);

        auto getStencils = [this, domainI]() -> const auto&
        {
            return (domainI == solidDomainIdx) ? couplingMapper().solidToVoidStencils() : couplingMapper().voidToSolidStencils();
        };

        const auto& stencils = getStencils();

        return (stencils.count(eIdx)) ? stencils.at(eIdx) : emptyStencil_;
    }

    /*!
     * \brief Returns whether a given solid grain/void pore body is coupled to the other domain
     */
    template<std::size_t i>
    bool isCoupledPore(Dune::index_constant<i> domainI, const std::size_t dofIdx) const
    {
        const auto& isCoupledDof = [&]
        {
            if constexpr(domainI == solidDomainIdx)
                return couplingMapper().isCoupledSolidDof();
            else //voidDomainIdx
                return couplingMapper().isCoupledVoidDof();
        }();
        return isCoupledDof[dofIdx];
    }

    /*!
    * \brief Returns summed conductive flux between void and solid for one void pore body or one solid grain
    */
    template<std::size_t i>
    Scalar conductionSource(Dune::index_constant<i> domainI,
                            const Element<i>& element,
                            const FVElementGeometry<i>& fvGeometry,
                            const ElementVolumeVariables<i>& elemVolVars,
                            const SubControlVolume<i>& scv) const
    {
        const auto& solidSol = this->curSol(solidDomainIdx);
        const auto& voidSol = this->curSol(voidDomainIdx);

        Scalar source = 0.0;
        const auto& connections = [&]
        {
            if constexpr(domainI == solidDomainIdx)
                return couplingMapper().solidToVoidConnections(scv.dofIndex());
            else //voidDomainIdx
                return couplingMapper().voidToSolidConnections(scv.dofIndex());
        }();

        // iterate over all connection throats
        for (const auto& connection : connections)
        {
            const Scalar t = getConnectionTransmissiblity(domainI, connection, elemVolVars, scv);
            const Scalar deltaT = [&]
            {
                if constexpr(domainI == solidDomainIdx)
                    return solidSol[scv.dofIndex()][Indices<solidDomainIdx>::temperatureIdx] - voidSol[connection.voidVertexIdx][Indices<voidDomainIdx>::temperatureIdx];
                else //voidDomainIdx
                    return voidSol[scv.dofIndex()][Indices<voidDomainIdx>::temperatureIdx] - solidSol[connection.solidVertexIdx][Indices<solidDomainIdx>::temperatureIdx];
            }();

            source -= t * deltaT;
        }

        source /= (scv.volume() * fvGeometry.gridGeometry().coordinationNumber()[scv.dofIndex()]);
        return source;
    }

    template<std::size_t i, class Connection, class Context>
    Scalar convectiveHeatFluxForOneConnection(Dune::index_constant<i> domainI,
                                              const Connection& connection,
                                              const Context& context) const
    {
        const auto& solidSol = this->curSol(solidDomainIdx);
        const auto& voidSol = this->curSol(voidDomainIdx);
        const Scalar tSolid = solidSol[connection.solidVertexIdx];
        Scalar resultConvection = 0.0;

        // iterate over all convection void elements
        for (const auto& [localConvectionVoidElementIdx, convectionVoidElementIdx] : enumerate(connection.convectionVoidElementIdx))
        {
            const auto& voidElement = voidGridGeometry_->element(convectionVoidElementIdx);
            const auto& voidFVGeometry = context.voidFVGeometry[localConvectionVoidElementIdx];
            const auto& voidElemVolVars = context.voidElemVolVars[localConvectionVoidElementIdx];
            const auto& voidElemFluxVarsCache = context.voidElemFluxVarsCache[localConvectionVoidElementIdx];

            const Scalar distance = [&, convectionVoidElementIdx = convectionVoidElementIdx, solidVertexIdx = connection.solidVertexIdx]
            {
                const auto& throatCenter = this->problem(voidDomainIdx).spatialParams().throatCenter(convectionVoidElementIdx);
                for (const auto& solidScv : scvs(context.solidFVGeometry[0]))
                {
                    if (solidScv.dofIndex() == solidVertexIdx)
                        return (solidScv.dofPosition() - throatCenter).two_norm();
                }
                DUNE_THROW(Dune::InvalidStateException, "No solid scv found");
            }();

            using VoidFluxVariables = GetPropType<SubDomainTypeTag<voidDomainIdx>, Properties::FluxVariables>;
            VoidFluxVariables fluxVars;
            const auto& scvf = voidFVGeometry.scvf(0/*localScvfIdx*/); //only one scvf per element -> localScvfIdx = 0
            fluxVars.init(this->problem(voidDomainIdx), voidElement, voidFVGeometry, voidElemVolVars, scvf, voidElemFluxVarsCache);

            static constexpr auto phaseIdx = 0;
            const Scalar flux = fluxVars.advectiveFlux(phaseIdx, [phaseIdx = phaseIdx](const auto& volVars){ return volVars.mobility(phaseIdx);});
            const Scalar velocity = flux / voidElemFluxVarsCache[scvf].throatCrossSectionalArea(phaseIdx);

            const auto tFluidMean = [&]
            {
                Scalar result = 0.0;
                const auto numScv = voidFVGeometry.numScv();
                assert(numScv == 2);
                for (const auto& voidScv : scvs(voidFVGeometry))
                    result += 1.0/numScv * voidSol[voidScv.dofIndex()][Indices<voidDomainIdx>::temperatureIdx];
                return result;
            };

            const auto tFluidUpstream = [&]
            {
                const auto upstreamDofIdx = flux > 0.0 ? voidFVGeometry.scv(scvf.insideScvIdx()).dofIndex() : voidFVGeometry.scv(scvf.outsideScvIdx()).dofIndex();
                return voidSol[upstreamDofIdx][Indices<voidDomainIdx>::temperatureIdx];
            };

            enum class FluidTemperatureMode {mean, self, upstream};
            static const auto fluidTemperatureMode = [&]
            {
                static const auto mode = getParam<std::string>("DualNetwork.FluidTemperatureMode", "mean");
                std::cout << "Using FluidTemperatureMode " << mode << std::endl;
                if (mode == "mean")
                    return FluidTemperatureMode::mean;
                else if (mode == "self")
                    return FluidTemperatureMode::self;
                else if (mode == "upstream")
                    return FluidTemperatureMode::upstream;
                else
                    DUNE_THROW(Dune::IOError, mode << " is not a valid FluidTemperatureMode");
            }();

            const Scalar tFluid = [&, voidVertexIdx = connection.voidVertexIdx]
            {
                if (fluidTemperatureMode == FluidTemperatureMode::mean)
                    return tFluidMean();
                else if (fluidTemperatureMode == FluidTemperatureMode::self)
                    return voidSol[voidVertexIdx][Indices<voidDomainIdx>::temperatureIdx];
                else
                    return tFluidUpstream();
            }();

            const Scalar meanKinematicViscosity = 0.5*(voidElemVolVars[0].viscosity(phaseIdx)/voidElemVolVars[0].density(phaseIdx)
                                                + voidElemVolVars[1].viscosity(phaseIdx)/voidElemVolVars[1].density(phaseIdx));
            const Scalar characteristicLength = 2.0*voidElemFluxVarsCache[scvf].throatInscribedRadius();
            using std::abs;
            const Scalar Re = DimLessNum::reynoldsNumber(abs(velocity), characteristicLength, meanKinematicViscosity);

            static const Scalar fixedLambda = getParam<Scalar>("DualNetwork.FixedConvectionLambda", -1.0);
            static const Scalar lambaFactor = getParam<Scalar>("DualNetwork.ConvectionLambaFactor", 0.9);
            static const Scalar lambaExponent = getParam<Scalar>("DualNetwork.ConvectionLambaExponent", 0.4);
            using std::pow;
            const Scalar lambda = fixedLambda > 0.0 ? fixedLambda : 1.0 + lambaFactor*pow(Re, lambaExponent); //approximation: see eq.30 in Koch et al (2021) https://doi.org/10.1007/s11242-021-01602-5
            const Scalar selfA = connection.connectionArea / connection.convectionVoidElementIdx.size();

            const auto neighborA = [&]
            {
                // get the area associated to the other void dof
                for (const auto& voidScv : scvs(voidFVGeometry))
                {
                    if (voidScv.dofIndex() != connection.voidVertexIdx)
                    {
                        const auto& otherConnections = couplingMapper().voidToSolidConnections(voidScv.dofIndex());
                        for (const auto& otherConn : otherConnections)
                        {
                            if (otherConn.solidVertexIdx == connection.solidVertexIdx)
                            {
                                if (otherConn.voidVertexIdx != voidScv.dofIndex())
                                    DUNE_THROW(Dune::InvalidStateException, "Void dofs don't match");

                                return otherConn.connectionArea/otherConn.convectionVoidElementIdx.size();
                            }
                        }
                    }
                }
                DUNE_THROW(Dune::InvalidStateException, "No neighbor area found");
            };

            static const bool useAvgA = getParam<bool>("DualNetwork.UseAverageConvectionArea", false);
            const Scalar A = useAvgA ? 0.5*(selfA + neighborA()) : selfA;

            const Scalar deltaT = (domainI == voidDomainIdx) ? (tSolid - tFluid) : (tFluid - tSolid);

            static const int verbose = getParam<int>("DualNetwork.SourceVerboseForDof", -1);
            if (verbose >= 0 && (connection.voidVertexIdx == verbose || connection.solidVertexIdx == verbose))
            {
                std::cout << "\n" << std::endl;
                const auto domain = (domainI == solidDomainIdx) ? "solid" : "void";
                std::cout << "At " << domain << ", dof " << verbose << ", flow elem " << convectionVoidElementIdx << ", globalId " << connection.id << std::endl;
                std::cout << std::scientific << std::setprecision(15) << "velocity " << velocity << ", Re "  << Re << ", delTaT " << deltaT << ", result " << lambda*A/distance*deltaT << std::endl;
                std::cout <<  std::scientific << std::setprecision(15) << "A " << A << ", distance " << distance << std::endl;
                std::cout << std::endl;
            }

            resultConvection += lambda*A/distance * deltaT;
        }

        return resultConvection;
    }

    /*!
    * \brief Returns summed conductive heat fluxes for one void pore body coupled to solid grains (or the other way around)
    *        that occur due to convection in void throats
    */
    template<std::size_t i>
    Scalar convectionSource(Dune::index_constant<i> domainI,
                            const Element<i>& element,
                            const FVElementGeometry<i>& fvGeometry,
                            const ElementVolumeVariables<i>& elemVolVars,
                            const SubControlVolume<i>& scv) const
    {
        static const int verbose = getParam<int>("DualNetwork.SourceVerboseForDof", -1);
        if (scv.dofIndex() == verbose)
            std::cout << "Start Source at elemn " << fvGeometry.gridGeometry().elementMapper().index(element) << " *******************************" <<  std::endl;

        // iterate over all connection throats
        const auto& connections = [&]
        {
            if constexpr (domainI == solidDomainIdx)
                return couplingMapper().solidToVoidConnections(scv.dofIndex());
            else
                return couplingMapper().voidToSolidConnections(scv.dofIndex());
        }();

        Scalar source = 0.0;
        bindCouplingContext(domainI, element);
        const auto& contextForPore = elementCouplingContext_()[scv.dofIndex()];

        for (const auto& [localConnectionIdx, connection] : enumerate(connections))
            source += convectiveHeatFluxForOneConnection(domainI, connection, contextForPore[localConnectionIdx]);

        if (scv.dofIndex() == verbose)
            std::cout << std::scientific << std::setprecision(15) << "total conv source  " << source << "\n\n ********************" << std::endl;

        source /= (scv.volume() * fvGeometry.gridGeometry().coordinationNumber()[scv.dofIndex()]);

        return source;
    }

    //TODO: this should not be in the general coupling manager
    template<std::size_t i, std::size_t j>
    Scalar sourceWithFixedTransmissibility(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           const FVElementGeometry<i>& fvGeometry,
                                           const ElementVolumeVariables<i>& elemVolVars,
                                           const SubControlVolume<i>& scv,
                                           Dune::index_constant<j> domainJ) const
    {
        const auto& voidSol = this->curSol(voidDomainIdx);
        const auto& solidSol = this->curSol(solidDomainIdx);

        Scalar source = 0.0;

        // iterate over all connection throats
        for (const auto& connection : couplingMapper().voidToSolidConnections(scv.dofIndex()))
        {
            const Scalar t = this->problem(voidDomainIdx).getInternalReferenceHeatTransmissibilityCoupling();
            const Scalar deltaT = [&]
            {
                if constexpr(domainI == solidDomainIdx)
                    return solidSol[scv.dofIndex()][Indices<solidDomainIdx>::temperatureIdx] - voidSol[connection.voidVertexIdx][Indices<voidDomainIdx>::temperatureIdx];
                else //voidDomainIdx
                    return voidSol[scv.dofIndex()][Indices<voidDomainIdx>::temperatureIdx] - solidSol[connection.solidVertexIdx][Indices<solidDomainIdx>::temperatureIdx];
            }();

            source -= t * deltaT;
        }

        source /= (scv.volume() * fvGeometry.gridGeometry().coordinationNumber()[scv.dofIndex()]);
        return source;
    }

    template<std::size_t i, class Connection>
    Scalar getConnectionTransmissiblity(Dune::index_constant<i> domainI,
                                        const Connection& connection,
                                        const ElementVolumeVariables<i>& elemVolVars,
                                        const SubControlVolume<i>& scv) const
    {
        static constexpr bool isVoid = (domainI == voidDomainIdx);

        const auto poreRadiusVoid = [&]
        {
            static const bool useExactPoreRadiusVoid = getParam<bool>("DualNetwork.UseExactPoreRadiusVoid", false);
            if (useExactPoreRadiusVoid)
            {
                using std::sqrt;
                static const Scalar R = getParam<Scalar>("DualNetwork.SphereRadius", 50e-6);
                static const Scalar overlapFactor = getParam<Scalar>("DualNetwork.OverlapFactor");
                static const Scalar dx = overlapFactor*R;
                static const Scalar r = sqrt(3.0) * dx - R;
                return r;
            }
            else
                return gridGeometry(voidDomainIdx).poreInscribedRadius(connection.voidVertexIdx);
        }();

        const auto poreRadiusSolid = [&]
        {
            static const bool useExactPoreRadiusSolid = getParam<bool>("DualNetwork.UseExactPoreRadiusSolid", false);
            if (useExactPoreRadiusSolid)
            {
                static const Scalar R = getParam<Scalar>("DualNetwork.SphereRadius", 50e-6);
                return R;
            }
            else
                return this->problem(solidDomainIdx).spatialParams().poreExtendedRadius(connection.solidVertexIdx);
        }();

        const Scalar fluidThermalConductivity = [&]
        {
            if constexpr (isVoid)
                return elemVolVars[scv].effectiveThermalConductivity();
            else
            {
                const auto& voidElement = voidGridGeometry_->element(connection.someVoidElementIdx);
                auto voidFVGeometry = localView(*voidGridGeometry_);
                voidFVGeometry.bindElement(voidElement);

                for (const auto& voidScv : scvs(voidFVGeometry))
                {
                    if (voidScv.dofIndex() == connection.voidVertexIdx)
                        return volVars(voidDomainIdx, voidElement, voidScv).effectiveThermalConductivity();
                }

                DUNE_THROW(Dune::InvalidStateException, "No scv found");
            }
        }();

        const Scalar solidThermalConductivity = [&]
        {
            if constexpr (!isVoid) //solid
                return elemVolVars[scv].effectiveThermalConductivity();
            else
            {
                const auto& solidElement = solidGridGeometry_->element(connection.someSolidElementIdx);
                auto solidFVGeometry = localView(*solidGridGeometry_);
                solidFVGeometry.bindElement(solidElement);

                for (const auto& solidScv : scvs(solidFVGeometry))
                {
                    if (solidScv.dofIndex() == connection.solidVertexIdx)
                        return volVars(solidDomainIdx, solidElement, solidScv).effectiveThermalConductivity();
                }

                DUNE_THROW(Dune::InvalidStateException, "No scv found");
            }
        }();

        const Scalar kappa = fluidThermalConductivity / solidThermalConductivity;
        static const Scalar Nu = getParam<Scalar>("DualNetwork.Nu", 1.0);
        static const Scalar Bi = getParam<Scalar>("DualNetwork.Bi", 1.0);

        static const bool useExactConnectionLength = getParam<bool>("DualNetwork.UseExactConnectionLength", false);
        const Scalar length = useExactConnectionLength ? poreRadiusSolid + poreRadiusVoid : connection.connectionLength;

        static const bool useExactConnectionAreaSphere = getParam<bool>("DualNetwork.UseExactConnectionAreaSphere", false);
        static const Scalar connectionAreaShapeFactor = getParam<Scalar>("DualNetwork.ConnectionAreaShapeFactor", 0.9);
        const Scalar area = [&]()
        {
            if (useExactConnectionAreaSphere)
            {
                static const Scalar R = getParam<Scalar>("DualNetwork.SphereRadius", 50e-6);
                static const Scalar overlapFactor = getParam<Scalar>("DualNetwork.OverlapFactor");
                static const auto dx = overlapFactor*R;
                static const auto h = R - dx;
                static const auto interfacialArea = 4*M_PI*R*R - 6*(2*M_PI*R*h);
                assert(std::abs(length - std::sqrt(3.0) * dx) < 1e-14); // TODO remove
                return interfacialArea/8.0*connectionAreaShapeFactor;
            }
            else
                return connection.connectionArea*connectionAreaShapeFactor;
        }();
        //distance weighted harmonic mean of thermal conductivities
        const Scalar meanThermalConductivity = (fluidThermalConductivity*Bi*(poreRadiusSolid + poreRadiusVoid))/(poreRadiusSolid*kappa + poreRadiusVoid*Bi/Nu);
        return area / length * meanThermalConductivity;
    }

    // \}

    /*!
    * \brief Return the volume variables of domain i for a given element and scv
    */
    template<std::size_t i>
    VolumeVariables<i> volVars(Dune::index_constant<i> domainI,
                                const Element<i>& element,
                                const SubControlVolume<i>& scv) const
    {
        VolumeVariables<i> volVars;
        const auto elemSol = elementSolution(element, this->curSol(domainI), gridGeometry(domainI));
        volVars.update(elemSol, this->problem(domainI), element, scv);
        return volVars;
    }

    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ)
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        typename LocalAssemblerI::LocalResidual::ElementResidualVector residual;

        const auto& element = localAssemblerI.element();
        const auto& fvGeometry = localAssemblerI.fvGeometry();
        const auto& curElemVolVars = localAssemblerI.curElemVolVars();

        residual.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
        {
            auto couplingSource = this->problem(domainI).source(element, fvGeometry, curElemVolVars, scv);
            couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
            residual[scv.indexInElement()] = couplingSource;
        }

        return residual;
    }

    /*!
     * \brief Build the coupling context for one pore of domain i
     *
     * The context holds everything needed to evaluate the convective coupling source of that
     * pore: for each of its connections the adjacent solid grain and the void throats carrying
     * the flow. It is returned by value and built from the current solution only, so that it
     * can be assembled concurrently on several threads.
     */
    template<std::size_t i>
    std::vector<CouplingContextForOneConnection>
    couplingContextForPore(Dune::index_constant<i> domainI, const std::size_t dofIdx) const
    {
        const auto& connections = [&]
        {
            if constexpr (domainI == solidDomainIdx)
                return couplingMapper().solidToVoidConnections(dofIdx);
            else
                return couplingMapper().voidToSolidConnections(dofIdx);
        }();

        std::vector<CouplingContextForOneConnection> contextForPore;
        contextForPore.reserve(std::distance(connections.begin(), connections.end()));

        auto solidFVGeometry = localView(*solidGridGeometry_);
        auto voidFVGeometry = localView(*voidGridGeometry_);
        auto voidElemVolVars = localView(gridVariables(voidDomainIdx).curGridVolVars());
        auto voidElemFluxVarsCache = localView(gridVariables(voidDomainIdx).gridFluxVarsCache());

        for (const auto& connection : connections)
        {
            CouplingContextForOneConnection context;
            context.connectionGlobalId = connection.id;

            const auto& solidElement = solidGridGeometry_->element(connection.someSolidElementIdx);
            solidFVGeometry.bindElement(solidElement);
            context.solidFVGeometry.push_back(solidFVGeometry);

            for (const auto& solidScv : scvs(solidFVGeometry))
                if (solidScv.dofIndex() == connection.solidVertexIdx)
                    context.solidVolVars = volVars(solidDomainIdx, solidElement, solidScv);

            const auto numConvectionVoidElems = connection.convectionVoidElementIdx.size();
            context.voidFVGeometry.reserve(numConvectionVoidElems);
            context.voidElemVolVars.reserve(numConvectionVoidElems);
            context.voidElemFluxVarsCache.reserve(numConvectionVoidElems);
            for (const auto voidElemIdx : connection.convectionVoidElementIdx)
            {
                const auto voidElem = voidGridGeometry_->element(voidElemIdx);
                voidFVGeometry.bindElement(voidElem);
                voidElemVolVars.bindElement(voidElem, voidFVGeometry, this->curSol(voidDomainIdx));
                voidElemFluxVarsCache.bindElement(voidElem, voidFVGeometry, voidElemVolVars);

                context.voidFVGeometry.push_back(voidFVGeometry);
                context.voidElemVolVars.push_back(voidElemVolVars);
                context.voidElemFluxVarsCache.push_back(voidElemFluxVarsCache);
            }

            contextForPore.push_back(std::move(context));
        }

        return contextForPore;
    }

    /*!
     * \brief Bind the coupling context of the calling thread to the given element
     * \note TODO remove Assembler
     */
    template<std::size_t i, class Assembler = int>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler = 0) const
    { elementCouplingContext_().bind(domainI, element, *this); }

    /*!
     * \brief Update the coupling context
     */
    template<std::size_t i, class LocalAssemblerI, std::size_t j, class PriVars>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PriVars& priVars,
                               int pvIdxJ)
    {
        this->curSol(domainJ)[dofIdxGlobalJ][pvIdxJ] = priVars[pvIdxJ];

        // patch up this thread's coupling context so that it reflects the deflected dof
        for (auto& context : elementCouplingContext_().data())
        {
            for (auto& connInfo : context.second)
            {
                const auto& staticConnectionInfo = couplingMapper().connectionInfo()[connInfo.connectionGlobalId];

                // when the solid dof is deflected, just update the solid volVars
                if constexpr (domainJ == solidDomainIdx)
                {
                    const auto& solidElement = gridGeometry(domainJ).element(staticConnectionInfo.someSolidElementIdx);
                    for (const auto& scv : scvs(connInfo.solidFVGeometry.front()))
                    {
                        if (scv.dofIndex() == dofIdxGlobalJ)
                        {
                            const auto elemSol = elementSolution(solidElement, this->curSol(domainJ), gridGeometry(domainJ));
                            connInfo.solidVolVars.update(elemSol, this->problem(domainJ), solidElement, scv);
                        }
                    }
                }
                else // deflection of void dofs
                {
                    assert(staticConnectionInfo.convectionVoidElementIdx.size() == connInfo.voidFVGeometry.size());
                    for (int voidElemLocalIdx = 0; voidElemLocalIdx < staticConnectionInfo.convectionVoidElementIdx.size(); ++voidElemLocalIdx)
                    {
                        const auto eIdx = staticConnectionInfo.convectionVoidElementIdx[voidElemLocalIdx];
                        const auto& element = gridGeometry(voidDomainIdx).element(eIdx);
                        const auto& fvGeometry = connInfo.voidFVGeometry[voidElemLocalIdx];
                        auto& elemVolVars = connInfo.voidElemVolVars[voidElemLocalIdx];

                        for (const auto& scv : scvs(fvGeometry))
                        {
                            // writing into the grid-wide volume variables cache (if the void
                            // domain has one) is safe under multithreaded assembly: the assembly
                            // coloring keeps elements writing the same void element's cached data
                            // in different colors, see computeColorsForAssembly()
                            if (scv.dofIndex() == dofIdxGlobalJ)
                            {
                                const auto elemSol = elementSolution(element, this->curSol(voidDomainIdx), gridGeometry(voidDomainIdx));
                                getVolVarAccess_(domainJ, gridVars_(voidDomainIdx).curGridVolVars(), elemVolVars, scv).update(elemSol, this->problem(voidDomainIdx), element, scv);
                            }
                        }

                        for (const auto& scvf : scvfs(fvGeometry))
                        {
                            if constexpr (voidFluxVarsCachingEnabled)
                                gridVars_(voidDomainIdx).gridFluxVarsCache().cache(eIdx, scvf.index()).update(problem(voidDomainIdx), element, fvGeometry, elemVolVars, scvf);
                            else
                                connInfo.voidElemFluxVarsCache[voidElemLocalIdx][scvf].update(this->problem(voidDomainIdx), element, fvGeometry, elemVolVars, scvf);
                        }
                    }
                }
            }
        }
    }

    template<std::size_t i>
    const auto& gridGeometry(Dune::index_constant<i> domainI) const
    {
        if constexpr (i == voidDomainIdx)
            return *voidGridGeometry_;
        else
            return *solidGridGeometry_;
    }

    const CouplingMapper& couplingMapper() const
    { return *couplingMapper_; }

    /*!
     * \brief set the pointers to the grid variables
     * \param gridVariables A tuple of shared pointers to the grid variables
     */
    void setGridVariables(GridVariablesTuple&& gridVariables)
    { gridVariables_ = gridVariables; }

    /*!
     * \brief set a pointer to one of the grid variables
     * \param gridVariables a pointer to the grid variables
     * \param domainIdx the domain index of the grid variables
     */
    template<class  GridVariables, std::size_t i>
    void setGridVariables(std::shared_ptr<GridVariables> gridVariables, Dune::index_constant<i> domainIdx)
    { std::get<i>(gridVariables_) = gridVariables; }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainIdx) const
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    /*!
     * \brief extend the jacobian pattern of the diagonal block of domain i
     *        by those entries that are not already in the uncoupled pattern
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {
        extendedSourceStencil_.extendJacobianPattern(*this, domainI, pattern);
    }

    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (per default this is not the case)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \note This is the same for box and cc
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector&,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {
        // std::cout << "calling evalAdditionalDomainDerivatives " << std::endl;
        extendedSourceStencil_.evalAdditionalDomainDerivatives(*this, domainI, localAssemblerI, this->curSol(domainI), A, gridVariables); //TODO: call without curSol, but currently protected .../dumux/dumux/multidomain/couplingmanager.hh
    }

    //! Return a reference to an empty stencil
    std::vector<std::size_t>& emptyStencil()
    { return emptyStencil_; }


    //! Return a reference to an empty stencil
    template<std::size_t i>
    const std::vector<std::size_t>& emptyStencil(Dune::index_constant<i> domainI) const
    { return emptyStencil_; }

    template<std::size_t i>
    const auto& gridView(Dune::index_constant<i> domainI) const
    {
        return gridGeometry(domainI).gridView();
    }


    void setupExtendedStencil()
    {
        extendedSourceStencil_.clear();

        for (const auto& element : elements(gridView(voidDomainIdx)))
        {
            const auto eIdx = gridGeometry(voidDomainIdx).elementMapper().index(element);
            const auto vIdx0 = gridGeometry(voidDomainIdx).vertexMapper().subIndex(element, 0, 1);
            const auto vIdx1 = gridGeometry(voidDomainIdx).vertexMapper().subIndex(element, 1, 1);

            for (const auto& is : intersections(gridView(voidDomainIdx), element))
            {
                if (is.neighbor())
                {
                    const auto& outsideElement = is.outside();
                    for (int i = 0; i < 2; ++i)
                    {
                        const auto outsideVertexIdx = gridGeometry(voidDomainIdx).vertexMapper().subIndex(outsideElement, i, 1);
                        if (outsideVertexIdx != vIdx0 && outsideVertexIdx != vIdx1)
                            extendedSourceStencil_.stencil()[eIdx].push_back(outsideVertexIdx);
                    }
                }
            }

            removeDuplicates_(extendedSourceStencil_.stencil()[eIdx]);
        }
    }


    /*!
     * \brief Compute colors for multithreaded assembly
     *
     * Two elements of the same subdomain may only be assembled concurrently if they never
     * write to the same memory. Assembling one element writes to
     *  -# its own matrix rows and the grid variables caches of its own stencil (the regular
     *     box stencil, which Dumux::Detail::addNeighborColors takes care of),
     *  -# the solution vector entry of every dof it deflects for the numeric differentiation,
     *     that is the dofs of the other subdomain in its coupling stencil and, in the void
     *     domain, additionally the dofs of its extended source stencil, and
     *  -# the cached void volume variables and void flux variables caches that its coupling
     *     context reads back, see updateCouplingContext().
     * Each of those makes a pair of elements conflicting, so none of them may share a color.
     */
    void computeColorsForAssembly()
    {
        Dune::Timer timer;
        using namespace Dune::Hybrid;
        static constexpr auto numDomains = Dune::index_constant<2>{};

        // (1) for each dof of a subdomain, the elements of that subdomain writing it
        auto connectedElements = std::make_tuple(
            Dumux::Detail::computeConnectedElements(gridGeometry(solidDomainIdx)),
            Dumux::Detail::computeConnectedElements(gridGeometry(voidDomainIdx))
        );

        // the extended source stencil makes a void element write dofs outside its regular stencil
        for (const auto& element : elements(gridView(voidDomainIdx)))
        {
            const auto eIdx = gridGeometry(voidDomainIdx).elementMapper().index(element);
            for (const auto dofIdx : extendedSourceStencil_.stencil().at(eIdx))
                addUnique_(std::get<voidDomainIdx>(connectedElements)[dofIdx], eIdx);
        }

        // (2) for each dof of a subdomain, the elements of the *other* subdomain deflecting it
        std::array<std::vector<std::vector<std::size_t>>, 2> coupledDofWriters;
        forEach(integralRange(numDomains), [&](const auto i)
        {
            constexpr auto j = Dune::index_constant<1-i>{};
            coupledDofWriters[j()].resize(gridGeometry(j).numDofs());
            for (const auto& element : elements(gridView(i)))
            {
                const auto eIdx = gridGeometry(i).elementMapper().index(element);
                for (const auto dofIdx : couplingStencil(i, element, j))
                    addUnique_(coupledDofWriters[j()][dofIdx], eIdx);
            }
        });

        // (3) for each void element, the elements of each subdomain writing its cached data.
        //     Without grid-wide caches in the void domain, updateCouplingContext() only ever
        //     touches the calling thread's own coupling context, so there is nothing to avoid.
        std::array<std::vector<std::vector<std::size_t>>, 2> voidCacheWriters;
        if constexpr (voidVolVarsCachingEnabled || voidFluxVarsCachingEnabled)
        {
            std::vector<std::size_t> writeSet;
            forEach(integralRange(numDomains), [&](const auto i)
            {
                voidCacheWriters[i()].resize(gridView(voidDomainIdx).size(0));
                for (const auto& element : elements(gridView(i)))
                {
                    const auto eIdx = gridGeometry(i).elementMapper().index(element);
                    voidCacheWriteSet_(i, element, writeSet);
                    for (const auto voidEIdx : writeSet)
                        voidCacheWriters[i()][voidEIdx].push_back(eIdx);
                }
            });
        }

        // greedily color each subdomain such that no conflicting elements share a color
        forEach(integralRange(numDomains), [&](const auto i)
        {
            constexpr auto j = Dune::index_constant<1-i>{};
            const auto& gg = gridGeometry(i);
            const auto& connected = std::get<i>(connectedElements);

            std::vector<int> colors(gg.gridView().size(0), -1);
            std::vector<int> neighborColors; neighborColors.reserve(200);
            std::vector<bool> colorUsed; colorUsed.reserve(200);
            std::vector<std::size_t> ownWriteSet;

            auto& elementSets = std::get<i>(elementSets_);
            elementSets.clear();

            for (const auto& element : elements(gg.gridView()))
            {
                const auto eIdx = gg.elementMapper().index(element);
                neighborColors.clear();

                // (1) everyone sharing a dof of our own stencil
                Dumux::Detail::addNeighborColors(gg, element, colors, connected, neighborColors);

                // (1b) everyone writing a dof of our extended source stencil (void domain only)
                if constexpr (i == voidDomainIdx)
                    for (const auto dofIdx : extendedSourceStencil_.stencil().at(eIdx))
                        for (const auto nIdx : connected[dofIdx])
                            neighborColors.push_back(colors[nIdx]);

                // (2) everyone deflecting a coupled dof that we deflect
                for (const auto dofIdx : couplingStencil(i, element, j))
                    for (const auto nIdx : coupledDofWriters[j()][dofIdx])
                        neighborColors.push_back(colors[nIdx]);

                // (3) everyone writing cached void data that we write
                if constexpr (voidVolVarsCachingEnabled || voidFluxVarsCachingEnabled)
                {
                    voidCacheWriteSet_(i, element, ownWriteSet);
                    for (const auto voidEIdx : ownWriteSet)
                        for (const auto nIdx : voidCacheWriters[i()][voidEIdx])
                            neighborColors.push_back(colors[nIdx]);
                }

                // assign the smallest color not taken by any conflicting element
                const auto color = Dumux::Detail::smallestAvailableColor(neighborColors, colorUsed);
                colors[eIdx] = color;

                if (static_cast<std::size_t>(color) < elementSets.size())
                    elementSets[color].push_back(element.seed());
                else
                    elementSets.push_back(std::vector<ElementSeed<i>>{ element.seed() });
            }
        });

        std::cout << Fmt::format("Colored dual network in {} seconds:\n", timer.elapsed());
        forEach(integralRange(numDomains), [&](const auto i)
        {
            std::cout << Fmt::format(
                "-- {} elements in subdomain {} with {} colors\n",
                gridView(i).size(0), i(), std::get<i>(elementSets_).size()
            );
        });
    }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainI the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainI, AssembleElementFunc&& assembleElement) const
    {
        if (std::get<i>(elementSets_).empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel: elements of the same color are guaranteed
        // not to write the same data (see computeColorsForAssembly), so each color can be
        // assembled by multiple threads while the colors themselves are done one after another
        const auto& grid = gridView(domainI).grid();
        for (const auto& elements : std::get<i>(elementSets_))
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t n)
            {
                const auto element = grid.entity(elements[n]);
                assembleElement(element);
            });
        }
    }

protected:

    template<std::size_t i>
    VolumeVariables<i>& getVolVarAccess_(Dune::index_constant<i> domainIdx, GridVolumeVariables<i>& gridVolVars, ElementVolumeVariables<i>& elemVolVars, const SubControlVolume<i>& scv)
    {
        if constexpr (getPropValue<SubDomainTypeTag<i>, Properties::EnableGridVolumeVariablesCache>())
            return gridVolVars.volVars(scv);
        else
            return elemVolVars[scv];
    }

    /*!
     * \brief Return a reference to the grid variables of a sub problem
     * \param domainIdx The domain index
     */
    template<std::size_t i>
    GridVariables<i>& gridVars_(Dune::index_constant<i> domainIdx)
    {
        if (std::get<i>(gridVariables_))
            return *std::get<i>(gridVariables_);
        else
            DUNE_THROW(Dune::InvalidStateException, "The gridVariables pointer was not set. Use setGridVariables() before calling this function");
    }

    //! Removes duplicate entries from the coupling stencils
    void removeDuplicates_(std::vector<std::size_t>& stencil)
    {
        std::sort(stencil.begin(), stencil.end());
        stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
    }

    //! Appends an element index to a list of element indices unless it is already in there
    static void addUnique_(std::vector<std::size_t>& elems, const std::size_t eIdx)
    {
        if (std::find(elems.begin(), elems.end(), eIdx) == elems.end())
            elems.push_back(eIdx);
    }

    /*!
     * \brief The void elements whose cached volume variables and flux variables caches are
     *        written while the residual of the given element of domain i is assembled
     * \note These are exactly the void elements that the coupling context of this element
     *       reads and that updateCouplingContext() therefore has to deflect.
     */
    template<std::size_t i>
    void voidCacheWriteSet_(Dune::index_constant<i> domainI,
                            const Element<i>& element,
                            std::vector<std::size_t>& writeSet) const
    {
        writeSet.clear();

        // each element has two pores, each of which may be coupled to the other domain
        for (int localVertexIdx = 0; localVertexIdx < 2; ++localVertexIdx)
        {
            const auto vIdx = gridGeometry(domainI).vertexMapper().subIndex(element, localVertexIdx, 1);
            if (!isCoupledPore(domainI, vIdx))
                continue;

            const auto& connections = [&]
            {
                if constexpr (domainI == solidDomainIdx)
                    return couplingMapper().solidToVoidConnections(vIdx);
                else
                    return couplingMapper().voidToSolidConnections(vIdx);
            }();

            for (const auto& connection : connections)
                for (const auto voidEIdx : connection.convectionVoidElementIdx)
                    writeSet.push_back(voidEIdx);
        }

        std::sort(writeSet.begin(), writeSet.end());
        writeSet.erase(std::unique(writeSet.begin(), writeSet.end()), writeSet.end());
    }


    bool convectiveHeatTransfer_;

    std::vector<std::size_t> emptyStencil_;
    std::unique_ptr<const CouplingMapper> couplingMapper_;
    const GridGeometry<voidDomainIdx>* voidGridGeometry_;
    const GridGeometry<solidDomainIdx>* solidGridGeometry_;


    /*!
     * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
     */
    GridVariablesTuple gridVariables_;

    //! the extended source stencil object
    PoreNetwork::PNMHeatExtendedSourceStencil<ThisType> extendedSourceStencil_;

    //! for each subdomain, the element seeds grouped by color for multithreaded assembly
    std::tuple<
        std::deque<std::vector<ElementSeed<solidDomainIdx>>>,
        std::deque<std::vector<ElementSeed<voidDomainIdx>>>
    > elementSets_;
};

} // end namespace Dumux::PoreNetwork

namespace Dumux {

/*!
 * \brief Multithreaded assembly is supported for the dual network coupling manager: the
 *        per-element coupling data is derived from the solution independently on each thread
 *        and the remaining shared data written during assembly (the coupled solution vector
 *        and the void domain's grid variables caches) is kept free of races by the assembly
 *        coloring, see PNMHeatTransferCouplingManager::computeColorsForAssembly().
 */
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<PoreNetwork::PNMHeatTransferCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
