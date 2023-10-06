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
    class Id {
        friend CCMpfaFluxes;
        Id(int i) : id{i} {}
        int id;
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

    /*!
     * \brief Register and prepare a new flux computation of the form \f$ f = -T(gradu + f) \f$
     * \param tensor Functor that returns a tensor (T in the above equation) per scv
     * \param values Functor that returns the values per scv (u in the above equation)
     * \param forces (optional) Functor that returns a force vector per scv (e.g. gravity; f in the above equation)
     * \return An identifier of type `Id` for thie registered flux computation.
     */
    Id add(const TensorAccessor& tensor,
           const ValueAccesor& values,
           const std::optional<ForceAccessor>& forces = {})
    { add_(tensor, values, forces); }

    //! Compute the flux for the given id & face
    Scalar computeFluxFor(const Id& id, const SubControlVolumeFace& scvf)
    { return computeFluxFor_(id.id, scvf); }

    // TODO: Transmissibility visitor or export?

private:
    virtual int add_(const TensorAccessor&, const ValueAccesor&, const std::optional<ForceAccessor>&) = 0;
    virtual Scalar computeFluxFor_(const int, const SubControlVolumeFace&) const = 0;
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
