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

#include <functional>
#include <dune/common/fmatrix.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class CCMpfaTransmissibilities
{
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using GridIndex = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dim = GridView::dimension;

public:
    class Id {
        friend CCMpfaTransmissibilities;
        Id(int i) : id{i} {}
        int id;
    };

    virtual ~CCMpfaTransmissibilities() = default;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using ScalarAccessor = std::function<Scalar(const GridIndex&)>;
    using TensorAccessor = std::function<Dune::FieldMatrix<Scalar, dim, dim>(const GridIndex&)>;
    using DofAccessor = std::function<Scalar(const GridIndex&)>;

    Id computeTransmissibilities(const ScalarAccessor& f) { return Id{computeTransmissibilities_(f)}; }
    Id computeTransmissibilities(const TensorAccessor& f) { return Id{computeTransmissibilities_(f)}; }
    Scalar computeFlux(const DofAccessor& a, const Id& id) { return computeFlux_(a, id.id); }
    // TODO: Transmissibility visitor or export?

protected:
    using DirichletBoundaryPredicate = std::function<bool(const Element&, const SubControlVolumeFace&)>;

private:
    virtual void bind_(const FVElementGeometry&, const DirichletBoundaryPredicate&) = 0;
    virtual int computeTransmissibilities_(const ScalarAccessor& f) = 0;
    virtual int computeTransmissibilities_(const TensorAccessor& f) = 0;
    virtual Scalar computeFlux_(const DofAccessor& a, int) const = 0;
};

// TODO: Scalar must come from outside
template<class GridGeometry, class Scalar>
class CCMpfaInteractionVolume
{
    using GridView = typename GridGeometry::GridView;

public:
    using Transmissibilities = CCMpfaTransmissibilities<GridGeometry, Scalar>;

    virtual ~CCMpfaInteractionVolume() = default;

    using GridIndex = typename GridView::IndexSet::IndexType;
    using GridIndexVisitor = std::function<void(const GridIndex&)>;

    std::unique_ptr<Transmissibilities> transmissibilities() const
    { return transmissibilities_(); }

    void visitGridScvIndices(const GridIndexVisitor& v) const
    { visitGridScvIndices_(v); }

    void visitGridScvfIndices(const GridIndexVisitor& v) const
    { visitGridScvfIndices_(v); }

private:
    virtual void visitGridScvIndices_(const GridIndexVisitor&) const = 0;
    virtual void visitGridScvfIndices_(const GridIndexVisitor&) const = 0;
    virtual std::unique_ptr<Transmissibilities> transmissibilities_() const = 0;
};

} // end namespace Dumux

#endif
