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
                const SubControlVolumeFace &scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeValues_); // do we need the shapeValues for the flux?

        // The stencil info is obsolete for the box method.
        // It is here for compatibility with cc methods
        stencil_ = Stencil(0);
    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

    const std::vector<ShapeValue>& shapeValues() const
    { return shapeValues_; }

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

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
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, element, fvGeometry, scvf);
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

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    // We always use the dynamic types here to be compatible on the boundary
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using TransmissibilityVector = typename BoundaryInteractionVolume::Vector;
    using PositionVector = typename BoundaryInteractionVolume::PositionVector;

public:
    MpfaAdvectionCache() { phaseNeumannFluxes_.fill(0.0); }

    //! update cached objects
    template<typename InteractionVolume>
    void updateAdvection(const SubControlVolumeFace &scvf,
                         const InteractionVolume& interactionVolume)
    {
        const auto& localFaceData = interactionVolume.getLocalFaceData(scvf);
        volVarsStencil_ = interactionVolume.volVarsStencil();
        volVarsPositions_ = interactionVolume.volVarsPositions();
        tij_ = interactionVolume.getTransmissibilities(localFaceData);
    }

    //! update cached neumann boundary flux
    template<typename InteractionVolume>
    void updatePhaseNeumannFlux(const SubControlVolumeFace &scvf,
                                const InteractionVolume& interactionVolume,
                                const unsigned int phaseIdx)
    {
        const auto& localFaceData = interactionVolume.getLocalFaceData(scvf);
        phaseNeumannFluxes_[phaseIdx] = interactionVolume.getNeumannFlux(localFaceData);
    }

    //! Returns the volume variables indices necessary for flux computation
    //! This includes all participating boundary volume variables. Since we
    //! do not allow mixed BC for the mpfa this is the same for all phases.
    const Stencil& advectionVolVarsStencil(const unsigned int phaseIdx) const
    { return volVarsStencil_; }

    //! Returns the position on which the volume variables live. This is
    //! necessary as we need to evaluate gravity also for the boundary volvars
    const PositionVector& advectionVolVarsPositions(const unsigned int phaseIdx) const
    { return volVarsPositions_; }

    //! Returns the transmissibilities associated with the volume variables
    //! All phases flow through the same rock, thus, tij are equal for all phases
    const TransmissibilityVector& advectionTij(const unsigned int phaseIdx) const
    { return tij_; }

    //! If the useTpfaBoundary property is set to false, the boundary conditions
    //! are put into the local systems leading to possible contributions on all faces
    Scalar advectionNeumannFlux(const unsigned int phaseIdx) const
    { return phaseNeumannFluxes_[phaseIdx]; }

private:
    // Quantities associated with advection
    Stencil volVarsStencil_;
    PositionVector volVarsPositions_;
    TransmissibilityVector tij_;
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

    // We always use the dynamic types here to be compatible on the boundary
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using TransmissibilityVector = typename BoundaryInteractionVolume::Vector;
    using PositionVector = typename BoundaryInteractionVolume::PositionVector;

public:
    // update cached objects for the diffusive fluxes
    template<typename InteractionVolume>
    void updateDiffusion(const SubControlVolumeFace &scvf,
                         const InteractionVolume& interactionVolume,
                         const unsigned int phaseIdx,
                         const unsigned int compIdx)
    {
        const auto& localFaceData = interactionVolume.getLocalFaceData(scvf);
        diffusionVolVarsStencils_[phaseIdx][compIdx] = interactionVolume.volVarsStencil();
        diffusionTij_[phaseIdx][compIdx] = interactionVolume.getTransmissibilities(localFaceData);
    }

    //! Returns the volume variables indices necessary for diffusive flux
    //! computation. This includes all participating boundary volume variables
    //! and it can be different for the phases & components.
    const Stencil& diffusionVolVarsStencil(const unsigned int phaseIdx,
                                           const unsigned int compIdx) const
    { return diffusionVolVarsStencils_[phaseIdx][compIdx]; }

    //! Returns the transmissibilities associated with the volume variables
    //! This can be different for the phases & components.
    const TransmissibilityVector& diffusionTij(const unsigned int phaseIdx,
                                               const unsigned int compIdx) const
    { return diffusionTij_[phaseIdx][compIdx]; }

private:
    // Quantities associated with molecular diffusion
    std::array< std::array<Stencil, numComponents>, numPhases> diffusionVolVarsStencils_;
    std::array< std::array<TransmissibilityVector, numComponents>, numPhases> diffusionTij_;
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
    using TransmissibilityVector = typename BoundaryInteractionVolume::Vector;
public:
    MpfaHeatConductionCache() : heatNeumannFlux_(0.0) {}

    // update cached objects for heat conduction
    template<typename InteractionVolume>
    void updateHeatConduction(const SubControlVolumeFace &scvf,
                              const InteractionVolume& interactionVolume)
    {
        const auto& localFaceData = interactionVolume.getLocalFaceData(scvf);
        heatConductionVolVarsStencil_ = interactionVolume.volVarsStencil();
        heatConductionTij_ = interactionVolume.getTransmissibilities(localFaceData);
    }

    // update cached objects for heat conduction
    template<typename InteractionVolume>
    void updateHeatNeumannFlux(const SubControlVolumeFace &scvf,
                               const InteractionVolume& interactionVolume)
    {
        const auto& localFaceData = interactionVolume.getLocalFaceData(scvf);
        heatNeumannFlux_ = interactionVolume.getNeumannFlux(localFaceData);
    }

    //! Returns the volume variables indices necessary for heat conduction flux
    //! computation. This includes all participating boundary volume variables
    //! and it can be different for the phases & components.
    const Stencil& heatConductionVolVarsStencil() const
    { return heatConductionVolVarsStencil_; }

    //! Returns the transmissibilities associated with the volume variables
    //! This can be different for the phases & components.
    const TransmissibilityVector& heatConductionTij() const
    { return heatConductionTij_; }

    //! If the useTpfaBoundary property is set to false, the boundary conditions
    //! are put into the local systems leading to possible contributions on all faces
    Scalar heatNeumannFlux() const
    { return heatNeumannFlux_; }

private:
    // Quantities associated with heat conduction
    Stencil heatConductionVolVarsStencil_;
    TransmissibilityVector heatConductionTij_;
    Scalar heatNeumannFlux_;
};

// specialization of the flux variables cache for cell centered mpfa methods
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
       : public MpfaPorousMediumFluxVariablesCache<TypeTag,
                                                   GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                   GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                   GET_PROP_VALUE(TypeTag, EnableEnergyBalance)> {};

// specialization for the case of pure advection
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, false, false>
             : public MpfaAdvectionCache<TypeTag>
{
    using AdvectionCache = MpfaAdvectionCache<TypeTag>;

public:
    //! the constructor
    MpfaPorousMediumFluxVariablesCache()
    : AdvectionCache(),
      isUpdated_(false)
    {}

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. Allows an update of the cache specific
    //! to processes that have solution dependent parameters, e.g. only updating
    //! the diffusion transmissibilities leaving the advective ones untouched
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

private:
    bool isUpdated_;
};

// specialization for the case of advection & diffusion
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, true, false>
             : public MpfaAdvectionCache<TypeTag>,
               public MpfaDiffusionCache<TypeTag>
{
    using AdvectionCache = MpfaAdvectionCache<TypeTag>;
    using DiffusionCache = MpfaDiffusionCache<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    // the constructor
    MpfaPorousMediumFluxVariablesCache()
    : AdvectionCache(),
      DiffusionCache(),
      isUpdated_(false)
    {}

    //! For compositional problems, neumann fluxes are not associated with a phase anymore
    //! TODO: How to implement neumann fluxes for !useTpfa
    template<typename InteractionVolume>
    void updatePhaseNeumannFlux(const SubControlVolumeFace &scvf,
                                const InteractionVolume& interactionVolume,
                                const unsigned int phaseIdx) {}

    //! TODO: How to implement neumann fluxes for !useTpfa
    Scalar advectionNeumannFlux(const unsigned int phaseIdx) const
    { return 0.0; }

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. Allows an update of the cache specific
    //! to processes that have solution dependent parameters, e.g. only updating
    //! the diffusion transmissibilities leaving the advective ones untouched
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

private:
    bool isUpdated_;
};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, false, true>
             : public MpfaAdvectionCache<TypeTag>,
               public MpfaHeatConductionCache<TypeTag>
{
    using AdvectionCache = MpfaAdvectionCache<TypeTag>;
    using HeatConductionCache = MpfaHeatConductionCache<TypeTag>;

public:
    // the constructor
    MpfaPorousMediumFluxVariablesCache()
    : AdvectionCache(),
      HeatConductionCache(),
      isUpdated_(false)
    {}

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. Allows an update of the cache specific
    //! to processes that have solution dependent parameters, e.g. only updating
    //! the diffusion transmissibilities leaving the advective ones untouched
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

private:
    bool isUpdated_;
};

// specialization for the case of advection, diffusion & heat conduction
template<class TypeTag>
class MpfaPorousMediumFluxVariablesCache<TypeTag, true, true, true>
             : public MpfaAdvectionCache<TypeTag>,
               public MpfaDiffusionCache<TypeTag>,
               public MpfaHeatConductionCache<TypeTag>
{
    using AdvectionCache = MpfaAdvectionCache<TypeTag>;
    using DiffusionCache = MpfaDiffusionCache<TypeTag>;
    using HeatConductionCache = MpfaHeatConductionCache<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    // the constructor
    MpfaPorousMediumFluxVariablesCache()
    : AdvectionCache(),
      DiffusionCache(),
      HeatConductionCache(),
      isUpdated_(false)
    {}

    //! TODO: How to implement neumann fluxes for !useTpfa when diffusion/heat conduction is active?
    template<typename InteractionVolume>
    void updatePhaseNeumannFlux(const SubControlVolumeFace &scvf,
                                const InteractionVolume& interactionVolume,
                                const unsigned int phaseIdx) {}


    //! This method is needed for compatibility reasons
    //! TODO: How to implement neumann fluxes for !useTpfa when diffusion is active?
    Scalar advectionNeumannFlux(const unsigned int phaseIdx) const
    { return 0.0; }

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. Allows an update of the cache specific
    //! to processes that have solution dependent parameters, e.g. only updating
    //! the diffusion transmissibilities leaving the advective ones untouched
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

private:
    bool isUpdated_;
};

} // end namespace

#endif
