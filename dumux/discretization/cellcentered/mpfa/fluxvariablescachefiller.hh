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
 * \brief The flux variables cache filler class
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux
{

//! forward declaration of properties
namespace Properties
{
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
NEW_PROP_TAG(ThermalConductivityModel);
};

/*!
 * \ingroup ImplicitModel
 * \brief Helper class to fill the flux var caches
 */
template<class TypeTag>
class CCMpfaFluxVariablesCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool doAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static constexpr bool doDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static constexpr bool doHeatConduction = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);

    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

    enum ProcessIndices : unsigned int
    {
        advectionIdx,
        diffusionIdx,
        heatConductionIdx
    };

public:
    //! The constructor. Sets the problem pointer
    CCMpfaFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param doSubCaches Array of bools indicating which sub caches have to be updated
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              const std::array<bool, 3>& doSubCaches = std::array<bool, 3>({true, true, true}))
    {
        // Set pointers
        elementPtr_ = &element;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        scvfPtr_ = &scvf;

        // prepare interaction volume and fill caches of all the scvfs connected to it
        const auto& globalFvGeometry = problem().model().globalFvGeometry();
        if (globalFvGeometry.isInBoundaryInteractionVolume(scvf))
        {
            bIv_ = std::make_unique<BoundaryInteractionVolume>(globalFvGeometry.boundaryInteractionVolumeSeed(scvf),
                                                               problem(),
                                                               fvGeometry,
                                                               elemVolVars);

            // fill the caches for all the scvfs in the interaction volume
            fillCachesInInteractionVolume_(fluxVarsCacheContainer, scvfFluxVarsCache, boundaryInteractionVolume(), doSubCaches);
        }
        else
        {
            iv_ = std::make_unique<InteractionVolume>(globalFvGeometry.interactionVolumeSeed(scvf),
                                                      problem(),
                                                      fvGeometry,
                                                      elemVolVars);

            // fill the caches for all the scvfs in the interaction volume
            fillCachesInInteractionVolume_(fluxVarsCacheContainer, scvfFluxVarsCache, interactionVolume(), doSubCaches);
        }
    }

    /*!
     * \brief function to update the flux variables caches during derivative calculation
     *
     * \copydoc fill
     */
    template<class FluxVariablesCacheContainer>
    void update(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                FluxVariablesCache& scvfFluxVarsCache,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        // array of bool with which we indicate the sub-caches which have to be
        // filled. During update, we only update solution-dependent quantities.
        static const std::array<bool, 3> doSubCaches = []()
        {
            std::array<bool, 3> doCaches;
            doCaches[ProcessIndices::advectionIdx] = doAdvection && soldependentAdvection;
            doCaches[ProcessIndices::diffusionIdx] = doDiffusion && soldependentDiffusion;
            doCaches[ProcessIndices::heatConductionIdx] = doHeatConduction && soldependentHeatConduction;
            return doCaches;
        } ();

        // forward to fill routine
        fill(fluxVarsCacheContainer, scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf, doSubCaches);
    }

    static bool isSolutionIndependent()
    {
        static const bool isSolDependent = (doAdvection && soldependentAdvection) ||
                                           (doDiffusion && soldependentDiffusion) ||
                                           (doHeatConduction && soldependentHeatConduction);
        return !isSolDependent;
    }

    const InteractionVolume& interactionVolume() const
    { return *iv_.get(); }

    const BoundaryInteractionVolume& boundaryInteractionVolume() const
    { return *bIv_.get(); }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvfPtr_; }

    InteractionVolume& interactionVolume()
    { return *iv_.get(); }

    BoundaryInteractionVolume& boundaryInteractionVolume()
    { return *bIv_.get(); }

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVariablesCacheContainer, class InteractionVolumeType>
    void fillCachesInInteractionVolume_(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                                        FluxVariablesCache& scvfFluxVarsCache,
                                        InteractionVolumeType& iv,
                                        const std::array<bool, 3>& doSubCaches)
    {
        // First we upate the interior boundary data and set the update status.
        // We store pointers to the other flux var caches and the elements they are embedded in simultaneously.
        // This way we have to obtain this data only once and can use it again in the sub-cache fillers.
        const auto numOtherScvfs = iv.globalLocalScvfPairedData().size()-1;
        std::vector<FluxVariablesCache*> otherFluxVarCaches(numOtherScvfs);
        std::vector<Element> otherElements(numOtherScvfs);

        scvfFluxVarsCache.updateInteriorBoundaryData(iv, scvFace());
        scvfFluxVarsCache.setUpdateStatus(true);

        const auto curScvfIdx = scvFace().index();
        unsigned int otherScvfIdx = 0;
        for (const auto& dataPair : iv.globalLocalScvfPairedData())
        {
            const auto& scvfJ = *dataPair.first;
            if (curScvfIdx == scvfJ.index())
                continue;

            // get the element scvfJ is embedded in
            const auto scvfJInsideScvIndex = scvfJ.insideScvIdx();
            otherElements[otherScvfIdx] =  scvfJInsideScvIndex == scvFace().insideScvIdx() ?
                                           element() :
                                           problem().model().globalFvGeometry().element(scvfJInsideScvIndex);

            // get the corresponding flux var cache
            otherFluxVarCaches[otherScvfIdx] = &fluxVarsCacheContainer[scvfJ];
            otherFluxVarCaches[otherScvfIdx]->updateInteriorBoundaryData(iv, scvfJ);
            otherFluxVarCaches[otherScvfIdx]->setUpdateStatus(true);
            otherScvfIdx++;
        }

        //! Maybe update the advective quantities
        if (doSubCaches[ProcessIndices::advectionIdx])
            fillAdvection(fluxVarsCacheContainer, scvfFluxVarsCache, iv, otherFluxVarCaches, otherElements);

        //! Maybe update the diffusive quantities
        if (doSubCaches[ProcessIndices::diffusionIdx])
            fillDiffusion(fluxVarsCacheContainer, scvfFluxVarsCache, iv, otherFluxVarCaches, otherElements);

        //! Maybe update quantities related to heat conduction
        if (doSubCaches[ProcessIndices::heatConductionIdx])
            fillHeatConduction(fluxVarsCacheContainer, scvfFluxVarsCache, iv, otherFluxVarCaches, otherElements);
    }

    //! method to fill the advective quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  FluxVariablesCache& scvfFluxVarsCache,
                  InteractionVolumeType& iv,
                  const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                  const std::vector<Element> otherElements)
    {
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        using AdvectionFiller = typename AdvectionType::CacheFiller;

        static constexpr auto AdvectionMethod = AdvectionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, AdvectionMethod>;

        // maybe solve the local system subject to K (if AdvectionType uses mpfa)
        if (AdvectionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getAdvectionLambda());

        // fill the caches of all scvfs within this interaction volume
        AdvectionFiller::fill(scvfFluxVarsCache, problem(), element(), fvGeometry(), elemVolVars(), scvFace(), *this);

        unsigned int otherScvfIdx = 0;
        const auto curScvfIdx = scvFace().index();
        for (const auto& dataPair : iv.globalLocalScvfPairedData())
        {
            const auto& scvfJ = *dataPair.first;
            if (curScvfIdx == scvfJ.index())
                continue;

            // fill corresponding cache
            AdvectionFiller::fill(*otherFluxVarCaches[otherScvfIdx],
                                  problem(),
                                  otherElements[otherScvfIdx],
                                  fvGeometry(),
                                  elemVolVars(),
                                  scvfJ,
                                  *this);
            otherScvfIdx++;
        }
    }

    //! do nothing if advection is not enabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool advectionEnabled = doAdvection>
    typename std::enable_if<!advectionEnabled>::type
    fillAdvection(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  FluxVariablesCache& scvfFluxVarsCache,
                  InteractionVolumeType& iv,
                  const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                  const std::vector<Element> otherElements)
    {}

    //! method to fill the diffusive quantities
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool diffusionEnabled = doDiffusion>
    typename std::enable_if<diffusionEnabled>::type
    fillDiffusion(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  FluxVariablesCache& scvfFluxVarsCache,
                  InteractionVolumeType& iv,
                  const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                  const std::vector<Element> otherElements)
    {
        using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
        using DiffusionFiller = typename DiffusionType::CacheFiller;

        static constexpr auto DiffusionMethod = DiffusionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, DiffusionMethod>;

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                if (phaseIdx == compIdx)
                    continue;

                // solve the local system subject to the diffusion tensor (if uses mpfa)
                if (DiffusionMethod == DiscretizationMethods::CCMpfa)
                    iv.solveLocalSystem(LambdaFactory::getDiffusionLambda(phaseIdx, compIdx));

                // fill the caches of all scvfs within this interaction volume
                DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element(), fvGeometry(), elemVolVars(), scvFace(), *this);

                unsigned int otherScvfIdx = 0;
                const auto curScvfIdx = scvFace().index();
                for (const auto& dataPair : iv.globalLocalScvfPairedData())
                {
                    const auto& scvfJ = *dataPair.first;
                    if (curScvfIdx == scvfJ.index())
                        continue;

                    // fill corresponding cache
                    DiffusionFiller::fill(*otherFluxVarCaches[otherScvfIdx],
                                          phaseIdx,
                                          compIdx,
                                          problem(),
                                          otherElements[otherScvfIdx],
                                          fvGeometry(),
                                          elemVolVars(),
                                          scvfJ,
                                          *this);
                    otherScvfIdx++;
                }
            }
        }
    }

    //! do nothing if diffusion is not enabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool diffusionEnabled = doDiffusion>
    typename std::enable_if<!diffusionEnabled>::type
    fillDiffusion(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  FluxVariablesCache& scvfFluxVarsCache,
                  InteractionVolumeType& iv,
                  const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                  const std::vector<Element> otherElements)
    {}

    //! method to fill the quantities related to heat conduction
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       FluxVariablesCache& scvfFluxVarsCache,
                       InteractionVolumeType& iv,
                       const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                       const std::vector<Element> otherElements)
    {
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
        using HeatConductionFiller = typename HeatConductionType::CacheFiller;

        static constexpr auto HeatConductionMethod = HeatConductionType::myDiscretizationMethod;
        using LambdaFactory = TensorLambdaFactory<TypeTag, HeatConductionMethod>;

        // maybe solve the local system subject to fourier coefficient
        if (HeatConductionMethod == DiscretizationMethods::CCMpfa)
            iv.solveLocalSystem(LambdaFactory::getHeatConductionLambda());

        // fill the caches of all scvfs within this interaction volume
        HeatConductionFiller::fill(scvfFluxVarsCache, problem(), element(), fvGeometry(), elemVolVars(), scvFace(), *this);

        unsigned int otherScvfIdx = 0;
        const auto curScvfIdx = scvFace().index();
        for (const auto& dataPair : iv.globalLocalScvfPairedData())
        {
            const auto& scvfJ = *dataPair.first;
            if (curScvfIdx == scvfJ.index())
                continue;

            // fill corresponding cache
            HeatConductionFiller::fill(*otherFluxVarCaches[otherScvfIdx],
                                       problem(),
                                       otherElements[otherScvfIdx],
                                       fvGeometry(),
                                       elemVolVars(),
                                       scvfJ,
                                       *this);
            otherScvfIdx++;
        }
    }

    //! do nothing if heat conduction is disabled
    template<class FluxVariablesCacheContainer, class InteractionVolumeType, bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<!heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                       FluxVariablesCache& scvfFluxVarsCache,
                       InteractionVolumeType& iv,
                       const std::vector<FluxVariablesCache*>& otherFluxVarCaches,
                       const std::vector<Element> otherElements)
    {}

    const Problem* problemPtr_;
    const Element* elementPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;
    const SubControlVolumeFace* scvfPtr_;

    // We store pointers to an inner and boundary interaction volume
    // these are updated during the filling of the caches and the
    // physics-related caches have access to them
    std::unique_ptr<InteractionVolume> iv_;
    std::unique_ptr<BoundaryInteractionVolume> bIv_;
};

} // end namespace

#endif
