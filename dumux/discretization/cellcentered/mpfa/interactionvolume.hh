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
#include <functional>

#include <dune/common/fmatrix.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Abstract interface for transmissibilities in mpfa methods.
 *        Computes and stores transmissibilities within interaction volumes.
 */
template<class GridGeometry, class S>
class CCMpfaTransmissibilities
{
    using GridView = typename GridGeometry::GridView;
    using GridIndex = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    class Id {
        friend CCMpfaTransmissibilities;
        Id(int i) : id{i} {}
        int id;
    };

    virtual ~CCMpfaTransmissibilities() = default;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using Scalar = S;
    using Tensor = Dune::FieldMatrix<S, dimWorld, dimWorld>;
    using TensorVariant = std::variant<Scalar, Tensor>;
    using TensorAccessor = std::function<TensorVariant(const SubControlVolume&)>;
    using DofAccessor = std::function<Scalar(const SubControlVolume&)>;

    //! Compute transmissibilities and return a unique identifier for the transmissibilities
    Id computeTransmissibilities(const TensorAccessor& f)
    { return Id{computeTransmissibilities_(f)}; }

    //! Compute the flux for the given face & transmissibilities id
    Scalar computeFlux(const SubControlVolumeFace& scvf, const DofAccessor& a, const Id& id)
    { return computeFlux_(scvf, a, id.id); }
    // TODO: Transmissibility visitor or export?

private:
    virtual int computeTransmissibilities_(const TensorAccessor&) = 0;
    virtual Scalar computeFlux_(const SubControlVolumeFace&, const DofAccessor&, int) const = 0;
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
    using GridIndex = typename GridView::IndexSet::IndexType;
    using GridIndexVisitor = std::function<void(const GridIndex&)>;

    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using DirichletBoundaryPredicate = std::function<bool(const Element&, const SubControlVolumeFace&)>;
    using Transmissibilities = CCMpfaTransmissibilities<GridGeometry, Scalar>;

    virtual ~CCMpfaInteractionVolume() = default;

    //! Return a transmissibility computation instance for this interaction volume.
    std::unique_ptr<Transmissibilities> transmissibilities(const FVElementGeometry& fvGeometry,
                                                           const DirichletBoundaryPredicate& dirichletPredicate) const
    { return transmissibilities_(fvGeometry, dirichletPredicate); }

    //! Visit the indices of the scvs involved in flux computations in this interaction volume.
    void visitGridScvIndices(const GridIndexVisitor& v) const
    { visitGridScvIndices_(v); }

    //! Visit the indices of the scvfs whose fluxes are computed in this interaction volume.
    void visitGridScvfIndices(const GridIndexVisitor& v) const
    { visitGridScvfIndices_(v); }

private:
    virtual void visitGridScvIndices_(const GridIndexVisitor&) const = 0;
    virtual void visitGridScvfIndices_(const GridIndexVisitor&) const = 0;
    virtual std::unique_ptr<Transmissibilities> transmissibilities_(const FVElementGeometry&,
                                                                    const DirichletBoundaryPredicate&) const = 0;
};

} // end namespace Dumux

#endif
