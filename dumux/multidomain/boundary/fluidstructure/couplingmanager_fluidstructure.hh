// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Coupling manager for fluid-structure interaction
 */
#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FLUIDSTRUCTURE_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FLUIDSTRUCTURE_COUPLINGMANAGER_HH

#include <memory>
#include <tuple>
#include <vector>
#include <array>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/assembly/coloring.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/multidomain/couplingmanager.hh>

#include <dumux/multidomain/boundary/fluidstructure/couplingmanager_structuremesh.hh>
#include <dumux/multidomain/boundary/fluidstructure/couplingmanager_cvfe_simple.hh>

namespace Dumux {

template<class MDTraits>
class FluidStructureCouplingManager
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using Indices
        = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>::VolumeVariables::Indices;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;

    using StructureMeshMotionMDTraits = MultiDomainTraits<SubDomainTypeTag<0>, SubDomainTypeTag<1>>;
    using FluidMDTraits = MultiDomainTraits<SubDomainTypeTag<2>, SubDomainTypeTag<3>>;

    template<std::size_t id>
    using SubSolutionVector = std::decay_t<decltype(std::declval<typename MDTraits::SolutionVector>()[Dune::index_constant<id>()])>;
    using SolutionVectors = typename MDTraits::template TupleOfSharedPtr<SubSolutionVector>;

    using CouplingStencil = std::vector<std::size_t>;
public:
    static constexpr auto fluidMomentumIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto fluidMassIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto structureIdx = typename MDTraits::template SubDomain<2>::Index();
    static constexpr auto meshMotionIdx = typename MDTraits::template SubDomain<3>::Index();

    //! export traits
    using MultiDomainTraits = MDTraits;

    FluidStructureCouplingManager()
    : structureMeshMotionCouplingManager_()
    , fluidCouplingManager_()
    {
        using namespace Dune::Hybrid;
        forEach(solutionVectors_, [&](auto&& solutionVector)
        {
            solutionVector = std::make_shared<typename std::decay_t<decltype(solutionVector)>::element_type>();
        });
    };

    FluidStructureCouplingManager(std::shared_ptr<GridGeometry<structureIdx>> structureGG,
                                  std::shared_ptr<GridGeometry<meshMotionIdx>> meshMotionGG,
                                  std::shared_ptr<GridGeometry<fluidMomentum>> fluidMomentumGG,
                                  std::shared_ptr<GridGeometry<fluidMass>> fluidMassGG)
    : structureMeshMotionCouplingManager_(structureGG, meshMotionGG)
    , fluidCouplingManager_(fluidMomentumGG, fluidMassGG)
    {
        // compute coupling stencils for fluid mesh motion and fluid structure
    }

    void init(std::shared_ptr<Problem<structureIdx>> structureProblem,
              std::shared_ptr<Problem<meshMotionIdx>> meshMotionProblem,
              std::shared_ptr<Problem<fluidMomentum>> fluidMomentumProblem,
              std::shared_ptr<Problem<fluidMass>> fluidMassProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(structureProblem, meshMotionProblem, fluidMomentumProblem, fluidMassProblem));

        using FFSol = typename CVFEFreeFlowCouplingManagerSimple<FluidMDTraits>::SolutionVectorStorage;
        fluidCouplingManager_.init(
            fluidMomentumProblem, fluidMassProblem,
            FFSol{ std::get<fluidMomentumIndex>(this->curSol()), std::get<fluidwMassIndex>(this->curSol()) }
        );

        using SMSol = typename StructureMeshCouplingManager<StructureMeshMotionMDTraits>::SolutionVectorStorage;
        structureMeshMotionCouplingManager_.init(
            structureProblem, meshMotionProblem,
            SMSol{ std::get<structureIdx>(this->curSol()), std::get<meshMotionIdx>(this->curSol()) }
        );
    }

    //! Update the solution vector before assembly
    void updateSolution(const typename MDTraits::SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(solutionVectors_)), [&](const auto id)
        {
            *std::get<id>(solutionVectors_) = curSol[id];
        });
    }

    /*!
     * \brief Return the coupling element stencil in domain J for a given element in domain I
     */
    template<std::size_t i, class Entity, std::size_t j>
    const auto& couplingStencil(Dune::index_constant<i> domainI,
                                const Entity& entity,
                                Dune::index_constant<j> domainJ) const
    {
        // if the domains are coupled according to the map, forward to sub-coupling manager
        if constexpr (domainI == fluidMassIdx && domainJ == fluidMomentumIdx)
            return fluidCouplingManager_.couplingStencil(domainI, entity, domainJ);
        else if constexpr (domainI == fluidMomentumIdx && domainJ == fluidMassIdx)
            return fluidCouplingManager_.couplingStencil(domainI, entity, domainJ);
        else if constexpr (domainI == structureIdx && domainJ == meshMotionIdx)
            return structureMeshMotionCouplingManager_.couplingStencil(domainI, entity, domainJ);
        else if constexpr (domainI == meshMotionIdx && domainJ == structureIdx)
            return structureMeshMotionCouplingManager_.couplingStencil(domainI, entity, domainJ);
        else
            return emptyStencil_;
    }

    // coupling quantity interfaces

    // coupling quantity forward interfaces
    auto displacement(typename GridGeometry<meshMotionIdx>::LocalView::Element const& element,
                      typename GridGeometry<meshMotionIdx>::SubControlVolume const& scv) const
    { return structureMeshMotionCouplingManager_.displacement(element, scv); }

    Scalar pressure(const Element<fluidMomentumIdx>& element,
                    const FVElementGeometry<fluidMomentumIdx>& fvGeometry,
                    const SubControlVolumeFace<fluidMomentumIdx>& scvf,
                    const bool considerPreviousTimeStep = false) const
    { return fluidCouplingManager_.pressure(element, fvGeometry, scvf, considerPreviousTimeStep); }

    Scalar pressure(const Element<fluidMomentumIdx>& element,
                    const FVElementGeometry<fluidMomentumIdx>& fvGeometry,
                    const SubControlVolume<fluidMomentumIdx>& scv,
                    const bool considerPreviousTimeStep = false) const
    { return fluidCouplingManager_.pressure(element, fvGeometry, scv, considerPreviousTimeStep); }

    auto faceVelocity(const Element<fluidMassIdx>& element,
                      const SubControlVolumeFace<fluidMassIdx>& scvf) const
    { fluidCouplingManager_.faceVelocity(element, scvf); }

    auto elementVelocity(const FVElementGeometry<fluidMassIdx>& fvGeometry) const
    { fluidCouplingManager_.elementVelocity(fvGeometry); }

protected:
    SolutionVectors& curSol()
    { return solutionVectors_; }

    const SolutionVectors& curSol() const
    { return solutionVectors_; }

private:
    StructureMeshCouplingManager<StructureMeshMotionMDTraits> structureMeshMotionCouplingManager_;
    CVFEFreeFlowCouplingManagerSimple<FluidMDTraits> fluidCouplingManager_;

    SolutionVectors solutionVectors_;

    CouplingStencil emptyStencil_;
};

// we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<FluidStructureCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
