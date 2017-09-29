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
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using DataHandle = typename T::DataHandle;
    using IndexSet = typename T::IndexSet;
    using LocalIndexContainer = typename T::DynamicLocalIndexContainer;
    using LocalIndexType = typename LocalIndexContainer::value_type;
    using GlobalIndexContainer = typename T::DynamicGlobalIndexContainer;
    using GlobalIndexType = typename GlobalIndexContainer::value_type;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // state the traits type publicly
    using Traits = T;

    class LocalFaceData
    {
        LocalIndexType ivLocalScvfIndex_;          //! the iv-local scvf index this scvf maps to
        LocalIndexType ivLocalInsideScvIndex_;     //! the iv-local index of the scvfs' inside scv
        LocalIndexType ivLocalOutsideScvfIndex_;   //! the index of this scvf in the iv-local outside faces
        LocalIndexType scvfLocalOutsideScvfIndex_; //! the index of this scvf in the scvf-local outside faces
        GlobalIndexType globalScvfIndex_;          //! the index of the corresponding global scvf
        bool isOutside_;                           //! indicates if this face maps to the iv-local index from "outside"

      public:
        //! Constructor for "inside" faces
        LocalFaceData(LocalIndexType faceIndex,
                      LocalIndexType scvIndex,
                      GlobalIndexType globalScvfIndex)
        : ivLocalScvfIndex_(faceIndex),
          ivLocalInsideScvIndex_(scvIndex),
          globalScvfIndex_(globalScvfIndex),
          isOutside_(false) {}

        //! Constructor for "outside" faces
        LocalFaceData(LocalIndexType faceIndex,
                      LocalIndexType scvIndex,
                      LocalIndexType indexInIvOutsideFaces,
                      LocalIndexType indexInScvfOutsideFaces,
                      GlobalIndexType globalScvfIndex)
        : ivLocalScvfIndex_(faceIndex),
          ivLocalInsideScvIndex_(scvIndex),
          ivLocalOutsideScvfIndex_(indexInIvOutsideFaces),
          scvfLocalOutsideScvfIndex_(indexInScvfOutsideFaces),
          globalScvfIndex_(globalScvfIndex),
          isOutside_(true) {}

        //! The index of the scvf within the inside faces
        LocalIndexType ivLocalScvfIndex() const { return ivLocalScvfIndex_; }
        LocalIndexType ivLocalInsideScvIndex() const { return ivLocalInsideScvIndex_; }
        LocalIndexType ivLocalOutsideScvfIndex() const { assert(isOutside_); return ivLocalOutsideScvfIndex_; }
        LocalIndexType scvfLocalOutsideScvfIndex() const { assert(isOutside_); return scvfLocalOutsideScvfIndex_; }
        GlobalIndexType globalScvfIndex() const { return globalScvfIndex_; }
        bool isOutside() const { return isOutside_; }
    };

    struct DirichletData
    {
        GlobalIndexType volVarIndex_;
        GlobalPosition ipGlobal_;

      public:
        DirichletData(const GlobalIndexType index, const GlobalPosition& ip)
        : volVarIndex_(index)
        , ipGlobal_(ip)
        {}

        const GlobalPosition& ipGlobal() const { return ipGlobal_; }
        GlobalIndexType volVarIndex() const { return volVarIndex_; }
    };

    using DirichletDataContainer = std::vector<DirichletData>;
    using LocalFaceDataContainer = std::vector<LocalFaceData>;

    //! Sets up the local scope (geometries etc) for a given iv index set!
    void setUpLocalScope(const IndexSet& indexSet,
                         const Problem& problem,
                         const FVElementGeometry& fvGeometry)
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a setUpLocalScope() method."); }

    //! sets the sizes of the corresponding matrices in the data handle
    void prepareDataHandle(DataHandle& dataHandle)
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a prepareDataHandle() method."); }

    //! solves for the transmissibilities subject to a given tensor
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor,
                          const Problem& problem,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          DataHandle& dataHandle)
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a solveLocalSystem() method."); }

    //! obtain the local data object for a given global scvf
    const LocalFaceData& getLocalFaceData(const SubControlVolumeFace& scvf) const
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a getLocalFaceData() method."); }

    //!returns a reference to the container with the local face data
    const LocalFaceDataContainer& localFaceData() const
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a localFaceData() method."); }

    //! returns a reference to the container with the data on Dirichlet boundaries
    const DirichletDataContainer& dirichletData() const
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a dirichletData() method."); }

    //! returns the indices of the volvars in the stencil of the interaction volume
    const GlobalIndexContainer& volVarsStencil() const
    { DUNE_THROW(Dune::NotImplemented, "Interaction volume implementation does not provide a volVarsStencil() method."); }

    //! returns the local index in a vector for a given global index
    template<typename IdxType1, typename IdxType2>
    LocalIndexType findIndexInVector(const std::vector<IdxType1>& vector, const IdxType2 globalIdx) const
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);
        assert(it != vector.end() && "could not find local index in the vector for the given global index!");
        return std::distance(vector.begin(), it);
    }
};

} // end namespace

#endif
