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
 * \ingroup CCMpfaDiscretization
 * \brief Base class for interaction volumes of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*
 * \ingroup CCMpfaDiscretization
 * \brief Type Traits to retrieve types associated with an implementation of a Dumux::CCMpfaInteractionVolumeBase.
 *        You have to specialize this class for every implementation of Dumux::CCMpfaInteractionVolumeBase.
 *
 * \code
 * //! export the interaction volume type to be used on boundaries etc.
 * using SecondaryInteractionVolume = ...;
 * //! export the type used for local indices
 * using LocalIndexType = ...;
 * //! export the type used for indices on the grid
 * using GridIndexType = ...;
 * //! export the type for the interaction volume index set
 * using IndexSet = ...;
 * //! export the type used for global coordinates
 * using GlobalPosition = ...;
 * //! export the type of interaction-volume local scvs
 * using LocalScvType = ...;
 * //! export the type of interaction-volume local scvfs
 * using LocalScvfType = ...;
 * //! export the type of used for the iv-local face data
 * using LocalFaceData = ...;
 * //! export the type of face data container (use dynamic container here)
 * using LocalFaceDataContainer = ...;
 * //! export the type used for iv-local matrices
 * using Matrix = ...;
 * //! export the type used for iv-local vectors
 * using Vector = ...;
 * //! export the type used for the iv-stencils
 * using Stencil = ...;
 * //! export the data handle type for this iv
 * using DataHandle = InteractionVolumeDataHandle< TypeTag, InteractionVolumeType >;
 * \endcode
 */
template< class InteractionVolume >
struct CCMpfaInteractionVolumeTraits {};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Base class for the interaction volumes of mpfa methods. It defines
 *        the interface and actual implementations should derive from this class.
 */
template<class TypeTag, class Implementation>
class CCMpfaInteractionVolumeBase
{
    // Curiously recurring template pattern
    Implementation & asImp() { return static_cast<Implementation&>(*this); }
    const Implementation & asImp() const { return static_cast<const Implementation&>(*this); }

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // Types required in the assembly of the local eq system
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    //! state the traits type publicly
    using Traits = CCMpfaInteractionVolumeTraits< Implementation >;

    //! Prepares everything for the assembly
    void setUpLocalScope(const typename Traits::IndexSet& indexSet,
                         const Problem& problem,
                         const FVElementGeometry& fvGeometry)
    { asImp().setUpLocalScope(); }

    //! returns the number of "primary" scvfs of this interaction volume
    std::size_t numFaces() const { return asImp().numFaces(); }

    //! returns the number of intermediate unknowns within this interaction volume
    std::size_t numUnknowns() const { return asImp().numUnknowns(); }

    //! returns the number of (in this context) known solution values within this interaction volume
    std::size_t numKnowns() const { return asImp().numKnowns(); }

    //! returns the number of scvs embedded in this interaction volume
    std::size_t numScvs() const { return asImp().numScvs(); }

    //! returns the number of scvfs embedded in this interaction volume
    std::size_t numScvfs() const { return asImp().numScvfs(); }

    //! returns a reference to the container with the local face data
    const typename Traits::LocalFaceDataContainer& localFaceData() const { asImp().localFaceData(); }

    //! returns the cell-stencil of this interaction volume
    const typename Traits::Stencil& stencil() const { return asImp().stencil(); }

    //! returns the local scvf entity corresponding to a given iv-local scvf idx
    const typename Traits::LocalScvfType& localScvf(typename Traits::LocalIndexType ivLocalScvfIdx) const
    { return asImp().localScvf(ivLocalScvfIdx); }

    //! returns the local scv entity corresponding to a given iv-local scv idx
    const typename Traits::LocalScvType& localScv(typename Traits::LocalIndexType ivLocalScvIdx) const
    { return asImp().localScv(ivLocalScvIdx); }

    //! returns the element in which the scv with the given local idx is embedded in
    const Element& element(typename Traits::LocalIndexType ivLocalScvIdx) const
    { return asImp().element(); }

    //! returns the number of interaction volumes living around a vertex
    template< class NodalIndexSet >
    static std::size_t numInteractionVolumesAtVertex(const NodalIndexSet& nodalIndexSet)
    { return Implementation::numInteractionVolumesAtVertex(nodalIndexSet); }

    //! adds the iv index sets living around a vertex to a given container
    //! and stores the the corresponding index in a map for each scvf
    template<class IvIndexSetContainer, class ScvfIndexMap, class NodalIndexSet>
    static void addInteractionVolumeIndexSets(IvIndexSetContainer& ivIndexSetContainer,
                                              ScvfIndexMap& scvfIndexMap,
                                              const NodalIndexSet& nodalIndexSet)
    { Implementation::addInteractionVolumeIndexSets(ivIndexSetContainer, scvfIndexMap, nodalIndexSet); }
};

} // end namespace Dumux

#endif
