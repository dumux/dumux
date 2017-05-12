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

#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
//! Base class for the interaction volume traits
template<class TypeTag>
class CCMpfaInteractionVolumeTraitsBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! Data types for the indices and index sets
    using LocalIndexType = std::uint8_t;
    using LocalIndexSet = std::vector<LocalIndexType>;
    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using GlobalIndexSet = std::vector<GlobalIndexType>;

    //! for network grids this means that we assume the tensors
    //! to be given in world coordinates! If a transformation of
    //! given data has to be performed, it has to be done in the
    //! spatial parameters method where the permeability is returned
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
};

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of mpfa methods.
 *        It defines the interface. Actual implementations should derive from this class.
 */
template<class TypeTag, typename T>
class CCMpfaInteractionVolumeBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);

    // export types from the traits class
    using LocalIndexType = typename T::LocalIndexType;
    using LocalIndexSet = typename T::LocalIndexSet;
    using GlobalIndexType = typename T::GlobalIndexType;
    using GlobalIndexSet = typename T::GlobalIndexSet;
    using PositionVector = typename T::PositionVector;
    using Tensor = typename T::Tensor;
    using Matrix = typename T::Matrix;
    using Vector = typename T::Vector;
    using DataHandle = typename T::DataHandle;
    using Seed = typename T::Seed;

public:
    // state the traits class for other classes to export types
    using Traits = T;

    struct LocalFaceData
    {
        LocalIndexType localScvfIndex;
        LocalIndexType localScvIndex;
        bool isOutside;

        //! default constructor
        LocalFaceData() = default;

        //! Constructor fully initializing the members
        LocalFaceData(const LocalIndexType faceIndex,
                      const LocalIndexType scvIndex,
                      bool isOut)
        : localScvfIndex(faceIndex),
          localScvIndex(scvIndex),
          isOutside(isOut) {}
    };

    using GlobalLocalFaceDataPair = std::pair<const SubControlVolumeFace*, LocalFaceData>;

    //! solves the local equation system for the computation of the transmissibilities
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a solveLocalSystem() method."); }

    //! returns the indices of the volvars in the stencil of the interaction volume
    const GlobalIndexSet& volVarsStencil() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a volVarsStencil() method."); }

    //! returns the positions corresponding to the volvars in the stencil of the interaction volume (cell centers or scvf.ipGlobal() on boundary)
    const PositionVector& volVarsPositions() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a volVarsPositions() method."); }

    //! returns the local index of an scvf in the IV and a boolean whether or not it is on the negative side of the local scvf (flux has to be inverted)
    LocalFaceData getLocalFaceData(const SubControlVolumeFace& scvf) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getLocalFaceData() method."); }

    //! returns the transmissibilities corresponding to a local scvf
    Vector getTransmissibilities(const LocalFaceData& localFaceData) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getTransmissibilities() method."); }

    //! returns the neumann flux corresponding to a local scvf
    Scalar getNeumannFlux(const LocalFaceData& localFaceData, unsigned int eqIdx) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getNeumannFlux() method."); }

    //! returns the local index in a vector for a given global index
    template<typename IdxType1, typename IdxType2>
    LocalIndexType findIndexInVector(const std::vector<IdxType1>& vector, const IdxType2 globalIdx) const
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);
        assert(it != vector.end() && "could not find local index in the vector for the given global index!");
        return std::distance(vector.begin(), it);
    }

    //! returns GlobalLocalFaceDataPair objects for the scvfs involved in this interaction volume
    const std::vector<GlobalLocalFaceDataPair>& globalLocalScvfPairedData() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a globalLocalScvfPairedData() method."); }
};

} // end namespace

#endif
