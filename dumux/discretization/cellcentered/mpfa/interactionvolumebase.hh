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
 * \brief Base class for interaction volumes of mpfa methods. Defines the interface.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH

#include <dumux/discretization/cellcentered/mpfa/interactionvolumeseed.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of mpfa methods.
 *        It defines the interface. Actual implementations should derive from this class.
 */
template<class TypeTag, typename Traits>
class CCMpfaInteractionVolumeBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    // some types to be exported
    using BoundaryInteractionVolume = typename Traits::BoundaryInteractionVolume;
    using LocalIndexType = typename Traits::LocalIndexType;
    using LocalIndexSet = typename Traits::LocalIndexSet;
    using LocalIndexPair = typename Traits::LocalIndexPair;
    using GlobalIndexType = typename Traits::GlobalIndexType;
    using GlobalIndexSet = typename Traits::GlobalIndexSet;
    using Vector = typename Traits::Vector;
    using PositionVector = typename Traits::PositionVector;
    using Seed = typename Traits::Seed;

    //! solves the local equation system for the computation of the transmissibilities
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a solveLocalSystem() method."); }

    //! returns the dof indices in the stencil of the interaction volume
    const GlobalIndexSet& stencil() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a stencil() method."); }

    //! returns the indices of the volvars in the stencil of the interaction volume
    const GlobalIndexSet& volVarsStencil() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a volVarsStencil() method."); }

    //! returns the positions corresponding to the volvars in the stencil of the interaction volume (cell centers or scvf.ipGlobal() on boundary)
    const PositionVector& volVarsPositions() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a volVarsPositions() method."); }

    //! returns a list of global scvf indices that are connected to this interaction volume
    GlobalIndexSet globalScvfs() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a globalScvfs() method."); }

    //! returns the local index of an scvf in the IV and a boolean whether or not it is on the negative side of the local scvf (flux has to be inverted)
    LocalIndexPair getLocalIndexPair(const SubControlVolumeFace& scvf) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getLocalIndexPair() method."); }

    //! returns the transmissibilities corresponding to a local scvf
    Vector getTransmissibilities(const LocalIndexPair& localIndexPair) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getTransmissibilities() method."); }

    //! returns the neumann flux corresponding to a local scvf
    Scalar getNeumannFlux(LocalIndexPair& localIndexPair) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getNeumannFlux() method."); }

    //! returns the local index in a vector for a given global index
    template<typename IdxType1, typename IdxType2>
    LocalIndexType findLocalIndex(const std::vector<IdxType1>& vector, const IdxType2 globalIdx) const
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);
        assert(it != vector.end() && "could not find local index in the vector for the given global index!");
        return std::distance(vector.begin(), it);
    }
};

} // end namespace

#endif
