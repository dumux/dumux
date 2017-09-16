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
 * \brief Base class for interaction volumes of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEBASE_HH

#include <dune/common/dynmatrix.hh>

#include <dumux/discretization/cellcentered/mpfa/interactionvolumedatahandle.hh>

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volume traits. The types stated here
 *        have to be defined in interaction volume traits. It is recommended
 *        that different implementations inherit from this class and overwrite the
 *        desired types or publicly state the typedef of this base class.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeTraitsBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    using LocalIndexType = std::uint8_t;

    //! The dynamic types are used e.g. by the mpfa-o method.
    //! To be compatible with schemes using both dynamic and static
    //! array types (e.g. L-method using mpfa-o interaction volumes
    //! on the boudaries), other classes interacting with the interaction
    //! volumes (e.g. flux vars cache) export the dynamic types. If your
    //! scheme is fully static on the entire grid, overwrite these traits.
    using DynamicLocalIndexContainer = std::vector<LocalIndexType>;
    using DynamicGlobalIndexContainer = std::vector<typename GridView::IndexSet::IndexType>;
    using DynamicMatrix = Dune::DynamicMatrix<Scalar>;
    using DynamicVector = typename DynamicMatrix::row_type;

    //! The data handle type. Uses the dynamic types as well per default...
    using DataHandle = InteractionVolumeDataHandle<TypeTag>;

    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using ScvBasis = std::array<GlobalPosition, dim>;

    //! for network grids this means that we assume the tensors
    //! to be given in world coordinates! If a transformation of
    //! given data has to be performed, it has to be done in the
    //! spatial parameters method where the tensor is returned or
    //! in the volume variables where it is stored
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
};

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of mpfa methods. It defines
 *        the interface and actual implementations should derive from this class.
 */
template<class TypeTag, typename T>
class CCMpfaInteractionVolumeBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using LocalIndexSet = typename T::DynamicLocalIndexContainer;
    using LocalIndexType = typename LocalIndexSet::value_type;
    using GlobalIndexSet = typename T::DynamicGlobalIndexContainer;
    using GlobalIndexType = typename GlobalIndexSet::value_type;
    using Vector = typename T::Vector;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // state the traits type publicly
    using Traits = T;

    struct LocalFaceData
    {
        LocalIndexType localScvfIndex;
        LocalIndexType localScvIndex;
        bool isOutside;

        //! Constructor
        LocalFaceData(LocalIndexType faceIndex, LocalIndexType scvIndex, bool isOut)
        : localScvfIndex(faceIndex),
          localScvIndex(scvIndex),
          isOutside(isOut) {}
    };

    struct DirichletData
    {
        GlobalIndexType volVarIndex;
        GlobalPosition ipGlobal;

        DirichletData(const GlobalIndexType index, const GlobalPosition& ip)
        : volVarIndex(index)
        , ipGlobal(ip)
        {}
    };

    using GlobalLocalFaceDataPair = std::pair<const SubControlVolumeFace*, LocalFaceData>;
    using DirichletDataContainer = std::vector<DirichletData>;

    //! solves the local equation system for the computation of the transmissibilities
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a solveLocalSystem() method."); }

    //! returns the indices of the volvars in the stencil of the interaction volume
    const GlobalIndexSet& volVarsStencil() const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a volVarsStencil() method."); }

    //! returns the local index of an scvf in the IV and a boolean whether or not it is on the negative side of the local scvf (flux has to be inverted)
    LocalFaceData getLocalFaceData(const SubControlVolumeFace& scvf) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getLocalFaceData() method."); }

    //! returns the transmissibilities corresponding to a local scvf
    Vector getTransmissibilities(const LocalFaceData& localFaceData) const
    { DUNE_THROW(Dune::NotImplemented, "Actual interaction volume implementation does not provide a getTransmissibilities() method."); }

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
