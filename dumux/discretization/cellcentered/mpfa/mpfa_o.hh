// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Implementation of abstract interfaces for the mpfa-o methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_HH

#include <optional>
#include <functional>

#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>

#include "interactionvolume.hh"
#include "dualgridindexset.hh"


// temporarily reuse old implementations
#include "omethod/interactionvolume.hh"
#include "omethod/localassembler.hh"
#include "omethod/ivassembler.hh"


namespace Dumux::CCMpfaO {

#ifndef DOXYGEN
namespace Detail {

template<int dim>
struct SpatialParams {
    template<class... Args>
    Dune::FieldVector<double, dim> gravity(Args&&...) const
    {
        Dune::FieldVector<double, dim> g;
        g = 0.0;
        g[dim-1] = -9.81;
    }
};

// Todo: delete once obsolete when refactoring has advanced
template<int dim, typename DirichletPredicate>
struct ProblemFacade
{
    ProblemFacade(const DirichletPredicate& p)
    : p_{p}
    {}

    template<typename Element, typename Scvf>
    Dumux::BoundaryTypes<1> boundaryTypes(const Element& e, const Scvf& scvf) const
    {
        Dumux::BoundaryTypes<1> bctypes;
        bctypes.setAllNeumann();
        if (p_(e, scvf))
            bctypes.setAllDirichlet();
        return bctypes;
    }

    SpatialParams<dim> spatialParams() const
    { return SpatialParams<dim>{}; }

private:
    const DirichletPredicate& p_;
};

template<int dim, int dimWorld, class Scalar>
class MatrixHandle
{
    static const bool isSurfaceGrid = dim < dimWorld;
    using DimVector = Dune::FieldVector<Scalar, dim>;

    using AMatrix = Dune::DynamicMatrix<Scalar>;
    using BMatrix = Dune::DynamicMatrix<Scalar>;
    using CMatrix = Dune::DynamicMatrix<Scalar>;
    using TMatrix = Dune::DynamicMatrix<Scalar>;
    using OutsideTij = std::vector<std::vector<Dune::DynamicVector<Scalar>>>; // TODO: size?
    using FaceOmegas = Dumux::ReservedVector<DimVector, 2>;
    using OmegaStorage = std::vector<FaceOmegas>;
    using FaceScalars = std::vector<Scalar>;

public:
    MatrixHandle(const bool withForces = false)
    {
        if constexpr (isSurfaceGrid)
            tijOutside_.emplace();
        if (withForces)
        {
            forces_.emplace();
            deltaForces_.emplace();
        }
    }

    const CMatrix& CAInverse() const { return CA_; }
    CMatrix& CAInverse() { return CA_; }

    const AMatrix& AInverse() const { return A_; }
    AMatrix& AInverse() { return A_; }

    const BMatrix& AInverseB() const { return AB_; }
    BMatrix& AInverseB() { return AB_; }

    const TMatrix& T() const { return T_; }
    TMatrix& T() { return T_; }

    const OmegaStorage& omegas() const { return wijk_; }
    OmegaStorage& omegas() { return wijk_; }

    const OutsideTij& tijOutside() const
    {
        if constexpr (!isSurfaceGrid)
            DUNE_THROW(Dune::InvalidStateException, "Outside tij only available on surface grids");
        return *tijOutside_;
    }

    OutsideTij& tijOutside()
    {
        if constexpr (!isSurfaceGrid)
            DUNE_THROW(Dune::InvalidStateException, "Outside tij only available on surface grids");
        return *tijOutside_;
    }

    const FaceScalars& forces() const
    {
        if (!forces_.has_value())
            DUNE_THROW(Dune::InvalidStateException, "No forces have been registered");
        return *forces_;
    }

    FaceScalars& forces()
    {
        if (!forces_.has_value())
            DUNE_THROW(Dune::InvalidStateException, "No forces have been registered");
        return *forces_;
    }

    const FaceScalars& deltaForces() const
    {
        if (!deltaForces_.has_value())
            DUNE_THROW(Dune::InvalidStateException, "No delta forces have been registered");
        return *deltaForces_;
    }

    FaceScalars& deltaForces()
    {
        if (!deltaForces_.has_value())
            DUNE_THROW(Dune::InvalidStateException, "No delta forces have been registered");
        return *deltaForces_;
    }

protected:
    OmegaStorage wijk_;
    TMatrix T_;
    AMatrix A_;
    BMatrix AB_;
    CMatrix CA_;
    std::optional<OutsideTij> tijOutside_;
    std::optional<FaceScalars> forces_;
    std::optional<FaceScalars> deltaForces_;
};

} // namespace Detail
#endif // DOXYGEN

template<class GridGeometry, class Scalar>
class Transmissibilities : public CCMpfaTransmissibilities<GridGeometry, Scalar>
{
    using ParentType = CCMpfaTransmissibilities<GridGeometry, Scalar>;
    using GridView = typename GridGeometry::GridView;
    using Handle = Detail::MatrixHandle<GridView::dimension, GridView::dimensionworld, Scalar>;

public:
    using typename ParentType::TensorAccessor;
    using typename ParentType::ForceAccessor;
    using typename ParentType::DofAccessor;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using DualGridNodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;

    template<class DirichletBoundaryPredicate>
    Transmissibilities(const DualGridNodalIndexSet& indexSet,
                       const FVElementGeometry& fvGeometry,
                       const DirichletBoundaryPredicate& dirichletPredicate)
    : indexSet_{indexSet}
    , fvGeometry_{fvGeometry}
    {
        iv_.bind(
            indexSet_,
            Detail::ProblemFacade<GridView::dimensionworld, DirichletBoundaryPredicate>{dirichletPredicate},
            fvGeometry_
        );

        // test assembly
        computeTransmissibilities_([] (const auto&) { return double{1.0}; }, {});
        handles_.clear();
    }

private:
    int computeTransmissibilities_(const TensorAccessor& t, const std::optional<ForceAccessor>& f) override
    {
        if (f.has_value())
            DUNE_THROW(Dune::NotImplemented, "Force assembly");
        handles_.emplace_back();
        assembleMatrices(handles_.back(), iv_, t, fvGeometry_);
        return handles_.size() - 1;
    }

    Scalar computeFlux_(const SubControlVolumeFace& scvf, const DofAccessor& u, int id) const override
    { DUNE_THROW(Dune::NotImplemented, ""); }

    const DualGridNodalIndexSet& indexSet_;
    const FVElementGeometry& fvGeometry_;
    CCMpfaOInteractionVolume<CCMpfaODefaultInteractionVolumeTraits<GridView, Scalar>> iv_;
    std::vector<Handle> handles_;
};

template<class GridGeometry, class Scalar>
class InteractionVolume : public CCMpfaInteractionVolume<GridGeometry, Scalar>
{
    using ParentType = CCMpfaInteractionVolume<GridGeometry, Scalar>;
    using GridView = typename GridGeometry::GridView;

public:
    using DualGridNodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using typename ParentType::GridIndexVisitor;
    using typename ParentType::Transmissibilities;
    using typename ParentType::DirichletBoundaryPredicate;
    using typename ParentType::FVElementGeometry;
    using typename ParentType::SubControlVolumeFace;

    InteractionVolume(const DualGridNodalIndexSet& indexSet)
    : indexSet_{indexSet}
    {
        this->size_.emplace(typename ParentType::Size{
            static_cast<unsigned int>(indexSet_.gridScvIndices().size()),  // numScvs
            static_cast<unsigned int>(indexSet_.gridScvfIndices().size()), // numTotalScvfs
            static_cast<unsigned int>(0),                                  // numAuxiliaryScvfs
            static_cast<unsigned int>(indexSet_.numBoundaryScvfs())        // numBoundaryScvfs
        });
    }

private:
    bool isFluxScvf_(const SubControlVolumeFace&) const override
    { return true; }

    void visitGridScvIndices_(const GridIndexVisitor& v) const override
    {
        for (const auto scvIdx : indexSet_.gridScvIndices())
            v(scvIdx);
    }

    void visitGridScvfIndices_(const GridIndexVisitor& v) const override
    {
        for (const auto scvfIdx : indexSet_.gridScvfIndices())
            v(scvfIdx);
    }

    void visitFluxGridScvfIndices_(const GridIndexVisitor& v) const override
    { visitGridScvfIndices_(v); }

    std::unique_ptr<CCMpfaTransmissibilities<GridGeometry, Scalar>> transmissibilities_(
        const FVElementGeometry& fvGeometry,
        const DirichletBoundaryPredicate& dirichletPredicate
    ) const
    {
        return std::make_unique<CCMpfaO::Transmissibilities<GridGeometry, Scalar>>(
            indexSet_, fvGeometry, dirichletPredicate
        );
    }

    const DualGridNodalIndexSet& indexSet_;
};

template<class GridGeometry, class Scalar>
class InteractionVolumeFactory : public CCMpfaInteractionVolumeFactory<GridGeometry, Scalar>
{
    using ParentType = CCMpfaInteractionVolumeFactory<GridGeometry, Scalar>;

public:
    using typename ParentType::DualGridNodalIndexSet;
    using typename ParentType::InteractionVolumesVisitor;

private:
    void visitInteractionVolumesAt_(const DualGridNodalIndexSet& ni, const InteractionVolumesVisitor& v) const
    { v(std::make_unique<InteractionVolume<GridGeometry, Scalar>>(ni)); }
};

} // end namespace Dumux::CCMpfaO

#endif
