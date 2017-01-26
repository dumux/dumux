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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH
#define DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class PorousMediumFluxVariablesCacheImplementation;

namespace Properties
{
// forward declaration
NEW_PROP_TAG(NumPhases);
}

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store flux stencils and data required for flux calculation
 */
template<class TypeTag>
using PorousMediumFluxVariablesCache = PorousMediumFluxVariablesCacheImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

// specialization for the Box Method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;
    using TransmissibilityVector = std::vector<IndexType>;

    using CoordScalar = typename GridView::ctype;
    static const int dim = GridView::dimension;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

public:

    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeValues_); // do we need the shapeValues for the flux?
    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

    const std::vector<ShapeValue>& shapeValues() const
    { return shapeValues_; }

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

    // The stencil info is obsolete for the box method.
    // It is here for compatibility with cc methods
    const Stencil& stencil() const
    {
        return stencil_;
    }

private:
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;

    Stencil stencil_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        stencil_ = FluxVariables::computeStencil(problem, element, fvGeometry, scvf);
        tij_ = AdvectionType::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
    }

    const Stencil& stencil() const
    { return stencil_; }

    const Scalar& tij() const
    { return tij_; }

private:
    Stencil stencil_;
    Scalar tij_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Classes building up the porous medium flux variables cache for mpfa methods
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// forward declaration of the base class of the mpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class MpfaPorousMediumFluxVariablesCache;

//! Base class for the advective cache in mpfa methods
template<class TypeTag>
class MpfaAdvectionCache
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    static constexpr bool facetCoupling = GET_PROP_VALUE(TypeTag, MpfaFacetCoupling);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    // We always use the dynamic types here to be compatible on the boundary
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;
    using PositionVector = typename BoundaryInteractionVolume::PositionVector;

public:
    //! update cached objects
    template<class InteractionVolume, class LocalFaceData>
    void updateAdvection(const InteractionVolume& interactionVolume,
                         const SubControlVolumeFace &scvf,
                         const LocalFaceData& scvfLocalFaceData)
    {
        // update the quantities that are equal for all phases
        advectionVolVarsStencil_ = interactionVolume.volVarsStencil();
        advectionVolVarsPositions_ = interactionVolume.volVarsPositions();
        advectionTij_ = interactionVolume.getTransmissibilities(scvfLocalFaceData);

        // we will need the neumann flux transformation only on interior boundaries with facet coupling
        if (enableInteriorBoundaries && facetCoupling)
            advectionCij_ = interactionVolume.getNeumannFluxTransformationCoefficients(scvfLocalFaceData);

        // The neumann fluxes always have to be set per phase
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            phaseNeumannFluxes_[phaseIdx] = interactionVolume.getNeumannFlux(scvfLocalFaceData, phaseIdx);
    }

    //! Returns the volume variables indices necessary for flux computation
    //! This includes all participating boundary volume variables. Since we
    //! do not allow mixed BC for the mpfa this is the same for all phases.
    const Stencil& advectionVolVarsStencil() const
    { return advectionVolVarsStencil_; }

    //! Returns the position on which the volume variables live. This is
    //! necessary as we need to evaluate gravity also for the boundary volvars
    const PositionVector& advectionVolVarsPositions() const
    { return advectionVolVarsPositions_; }

    //! Returns the transmissibilities associated with the volume variables
    //! All phases flow through the same rock, thus, tij are equal for all phases
    const CoefficientVector& advectionTij() const
    { return advectionTij_; }

    //! Returns the vector of coefficients with which the vector of neumann boundary conditions
    //! has to be multiplied in order to transform them on the scvf this cache belongs to
    const CoefficientVector& advectionCij() const
    { return advectionCij_; }

    //! If the useTpfaBoundary property is set to false, the boundary conditions
    //! are put into the local systems leading to possible contributions on all faces
    Scalar advectionNeumannFlux(unsigned int phaseIdx) const
    { return phaseNeumannFluxes_[phaseIdx]; }

private:
    // Quantities associated with advection
    Stencil advectionVolVarsStencil_;
    PositionVector advectionVolVarsPositions_;
    CoefficientVector advectionTij_;
    CoefficientVector advectionCij_;
    std::array<Scalar, numPhases> phaseNeumannFluxes_;
};

//! Base class for the diffusive cache in mpfa methods
template<class TypeTag>
class MpfaDiffusionCache
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    // We always use the dynamic types here to be compatible on the boundary
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;
    using PositionVector = typename BoundaryInteractionVolume::PositionVector;

public:
    //! The constructor. Initializes the Neumann flux to zero
    MpfaDiffusionCache() { componentNeumannFluxes_.fill(0.0); }

    // update cached objects for the diffusive fluxes
    template<typename InteractionVolume, class LocalFaceData>
    void updateDiffusion(const InteractionVolume& interactionVolume,
                         const SubControlVolumeFace &scvf,
                         const LocalFaceData& scvfLocalFaceData,
                         unsigned int phaseIdx,
                         unsigned int compIdx)
    {
        diffusionVolVarsStencils_[phaseIdx][compIdx] = interactionVolume.volVarsStencil();
        diffusionTij_[phaseIdx][compIdx] = interactionVolume.getTransmissibilities(scvfLocalFaceData);

        if (enableInteriorBoundaries)
            diffusionCij_[phaseIdx][compIdx] = interactionVolume.getNeumannFluxTransformationCoefficients(scvfLocalFaceData);

        //! For compositional models, we associate neumann fluxes with the phases (main components)
        //! This is done in the AdvectionCache. However, in single-phasic models we solve the phase AND
        //! the component mass balance equations. Thus, in this case we have diffusive neumann contributions.
        //! we assume compIdx = eqIdx
        if (numPhases == 1 && phaseIdx != compIdx)
            componentNeumannFluxes_[compIdx] = interactionVolume.getNeumannFlux(scvfLocalFaceData, compIdx);

    }

    //! Returns the volume variables indices necessary for diffusive flux
    //! computation. This includes all participating boundary volume variables
    //! and it can be different for the phases & components.
    const Stencil& diffusionVolVarsStencil(unsigned int phaseIdx,
                                           unsigned int compIdx) const
    { return diffusionVolVarsStencils_[phaseIdx][compIdx]; }

    //! Returns the transmissibilities associated with the volume variables
    //! This can be different for the phases & components.
    const CoefficientVector& diffusionTij(unsigned int phaseIdx,
                                          unsigned int compIdx) const
    { return diffusionTij_[phaseIdx][compIdx]; }

    //! Returns the vector of coefficients with which the vector of neumann boundary conditions
    //! has to be multiplied in order to transform them on the scvf this cache belongs to
    const CoefficientVector& diffusionCij(unsigned int phaseIdx,
                                          unsigned int compIdx) const
    { return diffusionCij_[phaseIdx][compIdx]; }

    //! If the useTpfaBoundary property is set to false, the boundary conditions
    //! are put into the local systems leading to possible contributions on all faces
    Scalar componentNeumannFlux(unsigned int compIdx) const
    {
        assert(numPhases == 1);
        return componentNeumannFluxes_[compIdx];
    }

private:
    // Quantities associated with molecular diffusion
    std::array< std::array<Stencil, numComponents>, numPhases> diffusionVolVarsStencils_;
    std::array< std::array<CoefficientVector, numComponents>, numPhases> diffusionTij_;
    std::array< std::array<CoefficientVector, numComponents>, numPhases> diffusionCij_;

    // diffusive neumann flux for single-phasic models
    std::array<Scalar, numComponents> componentNeumannFluxes_;
};

//! Base class for the heat conduction cache in mpfa methods
template<class TypeTag>
class MpfaHeatConductionCache
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    // We always use the dynamic types here to be compatible on the boundary
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;

    static constexpr int energyEqIdx = GET_PROP_TYPE(TypeTag, Indices)::energyEqIdx;
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);
public:
    // update cached objects for heat conduction
    template<typename InteractionVolume, class LocalFaceData>
    void updateHeatConduction(const InteractionVolume& interactionVolume,
                              const SubControlVolumeFace &scvf,
                              const LocalFaceData& scvfLocalFaceData)
    {
        heatConductionVolVarsStencil_ = interactionVolume.volVarsStencil();
        heatConductionTij_ = interactionVolume.getTransmissibilities(scvfLocalFaceData);
        heatNeumannFlux_ = interactionVolume.getNeumannFlux(scvfLocalFaceData, energyEqIdx);

        if (enableInteriorBoundaries)
            heatConductionCij_ = interactionVolume.getNeumannFluxTransformationCoefficients(scvfLocalFaceData);
    }

    //! Returns the volume variables indices necessary for heat conduction flux
    //! computation. This includes all participating boundary volume variables
    //! and it can be different for the phases & components.
    const Stencil& heatConductionVolVarsStencil() const
    { return heatConductionVolVarsStencil_; }

    //! Returns the transmissibilities associated with the volume variables
    //! This can be different for the phases & components.
    const CoefficientVector& heatConductionTij() const
    { return heatConductionTij_; }

    //! Returns the vector of coefficients with which the vector of neumann boundary conditions
    //! has to be multiplied in order to transform them on the scvf this cache belongs to
    const CoefficientVector& heatConductionCij() const
    { return heatConductionCij_; }

    //! If the useTpfaBoundary property is set to false, the boundary conditions
    //! are put into the local systems leading to possible contributions on all faces
    Scalar heatNeumannFlux() const
    { return heatNeumannFlux_; }

private:
    // Quantities associated with heat conduction
    Stencil heatConductionVolVarsStencil_;
    CoefficientVector heatConductionTij_;
    CoefficientVector heatConductionCij_;
    Scalar heatNeumannFlux_;
};

// specialization of the flux variables cache for cell centered mpfa methods
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
       : public MpfaPorousMediumFluxVariablesCache<TypeTag,
                                                   GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                   GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                   GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>
{
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ParentType = MpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                   GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                   GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

public:
    //! the constructor
    PorousMediumFluxVariablesCacheImplementation()
    : ParentType(),
      isUpdated_(false),
      interiorBoundaryInfoSelf_(false, -1)
    {}

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. This is used to only update
    //! the cache once after solving the local system. When visiting an scvf
    //! of the same interaction region again, the update is skipped.
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

    //! maybe update data on interior Dirichlet boundaries
    template<class InteractionVolume>
    void updateInteriorBoundaryData(const InteractionVolume& interactionVolume,
                                    const SubControlVolumeFace &scvf)
    {
        if (enableInteriorBoundaries)
        {
            interiorBoundaryData_ = interactionVolume.interiorBoundaryData();

            // check if the actual scvf is on an interior Dirichlet boundary
            const auto scvfIdx = scvf.index();
            unsigned int indexInData = 0;
            for (auto&& data : interiorBoundaryData_)
            {
                if (data.scvfIndex() == scvfIdx)
                {
                    interiorBoundaryInfoSelf_ = std::make_pair(true, indexInData);
                    break;
                }

                indexInData++;
            }
        }
    }

    bool isInteriorBoundary() const
    { return interiorBoundaryInfoSelf_.first; }

    const std::vector<InteriorBoundaryData>& interiorBoundaryData() const
    { return interiorBoundaryData_; }

    const InteriorBoundaryData& interiorBoundaryDataSelf() const
    {
        assert(interiorBoundaryInfoSelf_.first && "Trying to obtain interior boundary data on a face that is not marked as such");
        assert(interiorBoundaryInfoSelf_.second != -1 && "The index to the interior boundary data of this face has not been set");
        return interiorBoundaryData_[interiorBoundaryInfoSelf_.second];
    }

private:
    // indicates whether or not this cache has been fully updated
    bool isUpdated_;

    // if this face is an interior Dirichlet boundary itself, store additional data
    std::pair<bool, unsigned int> interiorBoundaryInfoSelf_;

    // contains all the interior Dirichlet boundary data within the stencil of this face
    std::vector<InteriorBoundaryData> interiorBoundaryData_;
};

// specialization for the case of pure advection
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, false, false> : public MpfaAdvectionCache<TypeTag> {};

// specialization for the case of advection & diffusion
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, true, false> : public MpfaAdvectionCache<TypeTag>,
                                                                       public MpfaDiffusionCache<TypeTag> {};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, false, true> : public MpfaAdvectionCache<TypeTag>,
                                                                       public MpfaHeatConductionCache<TypeTag> {};

// specialization for the case of advection, diffusion & heat conduction
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, true, true> : public MpfaAdvectionCache<TypeTag>,
                                                                      public MpfaDiffusionCache<TypeTag>,
                                                                      public MpfaHeatConductionCache<TypeTag> {};

} // end namespace

#endif
