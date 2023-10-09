// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Abstract interface for interaction volumes in mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUME_HH

#include <variant>
#include <optional>
#include <functional>
#include <memory>

#include <dune/common/fmatrix.hh>

#include "dualgridindexset.hh"

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Abstract interface for flux computations in mpfa methods.
 *        Computes and stores the data required for the assembly of
 *        fluxes across sub-control volume faces within an interaction volume.
 */
template<class GridGeometry, class S>
class CCMpfaFluxes
{
    using GridView = typename GridGeometry::GridView;
    using GridIndex = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    class TensorId
    {
        friend CCMpfaFluxes;
        int id = 0;
    public:
        TensorId() = default;
        TensorId(TensorId&&) = default;
        TensorId(const TensorId&) = default;
        TensorId(int i) : id(i) {}
        TensorId& operator=(const TensorId&) = default;
    };

    class FluxId
    {
        friend CCMpfaFluxes;
        TensorId tensorId{};
        int id = 0;
    public:
        FluxId() = default;
        FluxId(FluxId&&) = default;
        FluxId(const FluxId&) = default;
        FluxId(TensorId t, int i) : tensorId(std::move(t)), id(i) {}
        FluxId& operator=(const FluxId&) = default;
    };

    virtual ~CCMpfaFluxes() = default;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using Scalar = S;
    using Vector = Dune::FieldVector<S, dimWorld>;
    using Tensor = Dune::FieldMatrix<S, dimWorld, dimWorld>;
    using TensorVariant = std::variant<Scalar, Tensor>;
    using TensorAccessor = std::function<TensorVariant(const SubControlVolume&)>;
    using ForceAccessor = std::function<Vector(const SubControlVolume&)>;
    using ValueAccesor = std::function<Scalar(const SubControlVolume&)>;
    using BoundaryValueAccesor = std::function<Scalar(const SubControlVolumeFace&)>;

    /*!
     * \brief Get a copy of this fluxes instance.
     */
    std::unique_ptr<CCMpfaFluxes> clone() const
    { return clone_(); }

    /*!
     * \brief Register a tensor to be used in flux computations.
     * \param tensor Functor that returns a tensor (T in the above equation) per scv
     * \return An identifier of type `TensorId` for the registered tensor.
     */
    TensorId registerTensor(TensorAccessor&& tensor)
    { return {registerTensor_(std::move(tensor))}; }

    /*!
     * \brief Update a previously registered tensor.
     * \param id The id of the registered tensor.
     * \param tensor Functor that returns a tensor (T in the above equation) per scv.
     */
    void updateTensor(const TensorId& id, TensorAccessor&& tensor)
    { updateTensor_(id.id, std::move(tensor)); }

    /*!
     * \brief Register the values to be used for flux computations with a registered tensor.
     * \param tensorId The id of the tensor to be used in this flux expression
     * \param values Functor to obtain the values per scv
     * \param boundaryValues (optional) Functor to obtain the values at dirichlet boundary faces
     * \param forces (optional) Functor to obtain a force per scv occuring in the flux expression.
     * \return An identifier of type `FluxId` for the registered fluxes.
     */
    FluxId registerValuesFor(const TensorId& tensorId,
                             const ValueAccesor& values,
                             const std::optional<BoundaryValueAccesor>& boundaryValues = {},
                             const std::optional<ForceAccessor>& forces = {})
    { return {tensorId, registerValuesFor_(tensorId.id, values, boundaryValues, forces)}; }

    /*!
     * \brief Update registered values to be used for flux computations.
     * \param fluxId The id of the flux computation.
     * \param values Functor to obtain the values per scv
     * \param boundaryValues (optional) Functor to obtain the values at dirichlet boundary faces
     * \param forces (optional) Functor to obtain a force per scv occuring in the flux expression.
     * \return An identifier of type `FluxId` for the registered fluxes.
     */
    void updateValuesFor(const FluxId& fluxId,
                         const ValueAccesor& values,
                         const std::optional<BoundaryValueAccesor>& boundaryValues = {},
                         const std::optional<ForceAccessor>& forces = {})
    { updateValuesFor_(fluxId.id, fluxId.tensorId.id, values, boundaryValues, forces); }


    //! Compute the flux for the given id & face
    Scalar computeFluxFor(const FluxId& id, const SubControlVolumeFace& scvf) const
    { return computeFluxFor_(id.id, id.tensorId.id, scvf); }

    // TODO: Transmissibility visitor or export?

private:
    virtual std::unique_ptr<CCMpfaFluxes> clone_() const = 0;
    virtual int registerTensor_(TensorAccessor&&) = 0;
    virtual void updateTensor_(const int, TensorAccessor&&) = 0;
    virtual int registerValuesFor_(const int,
                                   const ValueAccesor&,
                                   const std::optional<BoundaryValueAccesor>&,
                                   const std::optional<ForceAccessor>&) = 0;
    virtual void updateValuesFor_(const int,
                                  const int,
                                  const ValueAccesor&,
                                  const std::optional<BoundaryValueAccesor>&,
                                  const std::optional<ForceAccessor>&) = 0;
    virtual Scalar computeFluxFor_(const int, const int, const SubControlVolumeFace&) const = 0;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Abstract interface for interaction volumes in mpfa methods.
 *        Allows extraction of a transmissibilities instance.
 */
template<class GridGeometry, class Scalar>
class CCMpfaInteractionVolume
{
    using GridView = typename GridGeometry::GridView;

public:
    struct Size
    {
        unsigned int numScvs;
        unsigned int numFluxScvfs;
        unsigned int numAuxiliaryScvfs;
        unsigned int numBoundaryScvfs;

        unsigned int numTotalScvfs() const
        { return numFluxScvfs + numAuxiliaryScvfs; }
    };

    virtual ~CCMpfaInteractionVolume() = default;

    using GridIndex = typename GridView::IndexSet::IndexType;
    using GridIndexVisitor = std::function<void(const GridIndex&)>;

    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using DirichletBoundaryPredicate = std::function<bool(const Element&, const SubControlVolumeFace&)>;
    using Fluxes = CCMpfaFluxes<GridGeometry, Scalar>;

    //! Return a transmissibility computation instance for this interaction volume.
    std::unique_ptr<Fluxes> fluxes(const FVElementGeometry& fvGeometry,
                                   const DirichletBoundaryPredicate& dirichletPredicate) const
    { return fluxes_(fvGeometry, dirichletPredicate); }

    //! Visit the indices of the scvs involved in flux computations in this interaction volume.
    void visitGridScvIndices(const GridIndexVisitor& v) const
    { visitGridScvIndices_(v); }

    //! Visit the indices of the all scvfs involved in flux computations.
    void visitGridScvfIndices(const GridIndexVisitor& v) const
    { visitGridScvfIndices_(v); }

    //! Visit the indices of the scvfs whose fluxes are computed in this interaction volume.
    void visitFluxGridScvfIndices(const GridIndexVisitor& v) const
    { visitFluxGridScvfIndices_(v); }

    //! Return true if this fluxes across this scvf are computed in this interaction volume.
    bool isFluxScvf(const SubControlVolumeFace& scvf) const
    { return isFluxScvf_(scvf); }

    //! Return the size info on this interaction volume
    const Size& size() const
    {
        if (!size_.has_value())
            DUNE_THROW(Dune::InvalidStateException, "Interaction volume implementation did not set the iv sizes");
        return *size_;
    }

protected:
    std::optional<Size> size_ = std::nullopt;

private:
    virtual bool isFluxScvf_(const SubControlVolumeFace&) const = 0;
    virtual void visitGridScvIndices_(const GridIndexVisitor&) const = 0;
    virtual void visitGridScvfIndices_(const GridIndexVisitor&) const = 0;
    virtual void visitFluxGridScvfIndices_(const GridIndexVisitor&) const = 0;
    virtual std::unique_ptr<Fluxes> fluxes_(const FVElementGeometry&, const DirichletBoundaryPredicate&) const = 0;
};

template<class GridGeometry, class Scalar>
class CCMpfaInteractionVolumeFactory
{
    using GridView = typename GridGeometry::GridView;

public:
    virtual ~CCMpfaInteractionVolumeFactory() = default;

    using DualGridNodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using InteractionVolume = CCMpfaInteractionVolume<GridGeometry, Scalar>;
    using InteractionVolumesVisitor = std::function<void(std::unique_ptr<InteractionVolume>&&)>;

    void visitInteractionVolumesAt(const DualGridNodalIndexSet& ni, const InteractionVolumesVisitor& v) const
    { visitInteractionVolumesAt_(ni, v); }

private:
    virtual void visitInteractionVolumesAt_(const DualGridNodalIndexSet&, const InteractionVolumesVisitor&) const = 0;
};

} // end namespace Dumux

#endif
