// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief Base class for the flux variables.
 */
#ifndef DUMUX_POROUSMEDIUM_FLUXVARIABLESCACHE_HH
#define DUMUX_POROUSMEDIUM_FLUXVARIABLESCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class PorousMediumFluxVariablesCacheImplementation;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup PorousmediumflowModels
 * \brief The flux variables cache classes for porous media.
 *
 * Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 * there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 * cache class are provided for different combinations of processes.
 */
template<class TypeTag>
using PorousMediumFluxVariablesCache = PorousMediumFluxVariablesCacheImplementation<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

//! We only store discretization-related quantities for control-volume finite element methods. Thus, we need no
//! physics-dependent specialization and simply inherit from the physics-independent implementation.
template<class TypeTag, class DM>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CVFE<DM>>
: public CVFEFluxVariablesCache<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::GridGeometry>>
{
public:
    //! export type used for scalar values
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
};

// the following classes choose the cache type: empty if the law disabled and the law's cache if it's enabled
// if advections is disabled the advection type is still instantiated if we use std::conditional_t and has to be a full type
// in order to prevent that instead of std::conditional_t we use this helper type is only dependent on the advection type
// if advection is enabled otherwise its an empty cache type
template<class TypeTag, bool EnableAdvection> class AdvectionCacheChooser : public FluxVariablesCaching::EmptyAdvectionCache {};
template<class TypeTag> class AdvectionCacheChooser<TypeTag, true> : public GetPropType<TypeTag, Properties::AdvectionType>::Cache {};
template<class TypeTag, bool EnableMolecularDiffusion> class DiffusionCacheChooser : public FluxVariablesCaching::EmptyDiffusionCache {};
template<class TypeTag> class DiffusionCacheChooser<TypeTag, true> : public GetPropType<TypeTag, Properties::MolecularDiffusionType>::Cache {};
template<class TypeTag, bool EnableEnergyBalance> class EnergyCacheChooser : public FluxVariablesCaching::EmptyHeatConductionCache {};
template<class TypeTag> class EnergyCacheChooser<TypeTag, true> : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache {};


// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
: public AdvectionCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableAdvection()>
, public DiffusionCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableMolecularDiffusion()>
, public EnergyCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance()>
{
public:
    //! export type used for scalar values
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
};

//! Specialization of the flux variables cache for the cell centered finite volume mpfa scheme.
//! Stores data which is commonly used by all the different types of processes.
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
: public AdvectionCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableAdvection()>
, public DiffusionCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableMolecularDiffusion()>
, public EnergyCacheChooser<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance()>
{
    using GridIndexType = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::IndexSet::IndexType;

    using MpfaHelper = typename GetPropType<TypeTag, Properties::GridGeometry>::MpfaHelper;
    static constexpr bool considerSecondary = MpfaHelper::considerSecondaryIVs();
public:
    //! export type used for scalar values
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    //! Returns whether or not this cache has been updated
    bool isUpdated() const { return isUpdated_; }

    //! Returns whether or not this cache is associated with a secondary interaction volume
    //! Specialization for deactivated secondary interaction volumes
    template< bool doSecondary = considerSecondary, std::enable_if_t<!doSecondary, int> = 0>
    constexpr bool usesSecondaryIv() const { return false; }

    //! Returns whether or not this cache is associated with a secondary interaction volume
    //! Specialization for activated secondary interaction volumes
    template< bool doSecondary = considerSecondary, std::enable_if_t<doSecondary, int> = 0>
    bool usesSecondaryIv() const { return usesSecondaryIv_; }

    //! Returns the index of the iv (this scvf is embedded in) in its container
    GridIndexType ivIndexInContainer() const { return ivIndexInContainer_; }
    //! Returns interaction volume-local face index
    unsigned int ivLocalFaceIndex() const { return ivLocalFaceIdx_; }
    //! Returns index of the face among "outside" faces of iv-local "positive" face
    unsigned int indexInOutsideFaces() const { return idxInOutsideFaces_; }

    //! Sets the update status. When set to true, consecutive updates will be skipped
    void setUpdateStatus(bool status) { isUpdated_ = status; }
    //! Sets if this cache is associated with a secondary iv
    void setSecondaryIvUsage(bool status) { usesSecondaryIv_ = status; }
    //! Sets the index of the iv (this scvf is embedded in) in its container
    void setIvIndexInContainer(GridIndexType ivIndex) { ivIndexInContainer_ = ivIndex; }
    //! Sets the iv-local face index
    void setIvLocalFaceIndex(unsigned int idx) { ivLocalFaceIdx_ = idx; }
    //! Sets the index of the face among the "positive" face's outside scvfs
    void setIndexInOutsideFaces(unsigned int idx) { idxInOutsideFaces_ = idx; }

private:

    bool isUpdated_ = false;       // returns true if cache has been fully updated
    bool usesSecondaryIv_ = false; // returns true if scvf is embedded in secondary interaction volume

    GridIndexType ivIndexInContainer_;  // index of the iv (this scvf is embedded in) in its container
    unsigned int ivLocalFaceIdx_;       // the interaction volume-local face index of this scvf
    unsigned int idxInOutsideFaces_;    // index of scvf among outside scvfs of iv-local "positive" face (only surface grids)
};

} // end namespace Dumux

#endif
