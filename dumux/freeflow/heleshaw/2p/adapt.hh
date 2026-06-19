// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Adaptive grid refinement utilities for the Hele-Shaw Darcy-Cahn-Hilliard
// two-phase model. Single-domain Box version — no velocity transfer needed.
//
//  - Indicator : refines where |Δφ| across faces exceeds a threshold
//  - DataTransfer : mass-lumped L² projection for φ and μ (conserves ∫φ dV),
//                   plain interpolation for p (not an extensive quantity)
//
#ifndef DUMUX_FREEFLOW_HELESHAW_2P_ADAPT_HH
#define DUMUX_FREEFLOW_HELESHAW_2P_ADAPT_HH

#include <limits>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/adaptive/griddatatransfer.hh>

namespace Dumux {

/*!
 * \brief Refinement indicator driven by jumps of φ across element faces.
 * Returns +1 to refine, -1 to coarsen, 0 to keep.
 */
template<class TypeTag>
class HeleShawTwoPGridAdaptIndicator
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView     = typename GridGeometry::GridView;
    using Element      = typename GridView::template Codim<0>::Entity;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using Indices      = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int indicatorIdx = Indices::phaseFieldIdx;

public:
    HeleShawTwoPGridAdaptIndicator(std::shared_ptr<const GridGeometry> gg,
                                   const std::string& paramGroup = "")
    : gg_(gg)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , delta_(gg->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    , minElementSize_(getParamFromGroup<Scalar>(paramGroup, "Adaptive.MinElementSize",
                                                std::numeric_limits<Scalar>::max()))
    {}

    void calculate(const SolutionVector& sol,
                   Scalar refineTol = 0.05, Scalar coarsenTol = 0.001)
    {
        refineBound_  = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        delta_.assign(gg_->gridView().size(0), 0.0);

        if (minLevel_ >= maxLevel_)
            return;

        for (const auto& element : elements(gg_->gridView()))
        {
            const auto eIdxI  = gg_->elementMapper().index(element);
            const auto geomI  = element.geometry();
            const auto elemSolI = elementSolution(element, sol, *gg_);
            const Scalar cI = evalSolution(element, geomI, *gg_, elemSolI, geomI.center())[indicatorIdx];

            using std::max, std::abs;
            for (const auto& is : intersections(gg_->gridView(), element))
            {
                if (!is.neighbor()) continue;
                const auto& outside = is.outside();
                const auto eIdxJ = gg_->elementMapper().index(outside);

                if (element.level() > outside.level()
                    || (element.level() == outside.level() && eIdxI < eIdxJ))
                {
                    const auto geomJ   = outside.geometry();
                    const auto elemSolJ = elementSolution(outside, sol, *gg_);
                    const Scalar cJ = evalSolution(outside, geomJ, *gg_, elemSolJ, geomJ.center())[indicatorIdx];
                    const Scalar loc = abs(cI - cJ);
                    delta_[eIdxI] = max(delta_[eIdxI], loc);
                    delta_[eIdxJ] = max(delta_[eIdxJ], loc);
                }
            }
        }

        refineBound_  = refineTol;   // φ ∈ [-1,1] so global scale = 1
        coarsenBound_ = coarsenTol;

        for (const auto& e : elements(gg_->gridView(), Dune::Partitions::interior))
            if (operator()(e) > 0)
                enforce2to1_(e);
    }

    int operator()(const Element& e) const
    {
        const auto idx = gg_->elementMapper().index(e);
        if (e.hasFather() && delta_[idx] < coarsenBound_)                   return -1;
        if (e.level() < (int)maxLevel_ && diameter_(e) > minElementSize_
            && delta_[idx] > refineBound_)                                   return  1;
        return 0;
    }

private:
    void enforce2to1_(const Element& e, std::size_t depth = 1)
    {
        for (const auto& is : intersections(gg_->gridView(), e))
        {
            if (!is.neighbor()) continue;
            const auto& nb = is.outside();
            if (nb.partitionType() == Dune::GhostEntity) continue;
            if (diameter_(nb) > minElementSize_
                && diameter_(nb) > 2.0 * diameter_(e)
                && nb.level() < (int)maxLevel_)
            {
                delta_[gg_->elementMapper().index(nb)] = std::numeric_limits<Scalar>::max();
                if (depth < maxLevel_) enforce2to1_(nb, depth + 1);
            }
        }
    }

    Scalar diameter_(const Element& e) const { return diameter(e.geometry()); }

    std::shared_ptr<const GridGeometry> gg_;
    Scalar refineBound_, coarsenBound_;
    std::vector<Scalar> delta_;
    std::size_t minLevel_, maxLevel_;
    Scalar minElementSize_;
};


/*!
 * \brief Data transfer for Box discretization of the Hele-Shaw 2p model.
 *
 * Primary variables: [0] p, [1] φ, [2] μ
 *   - p  (idx 0): plain interpolation (not an extensive quantity)
 *   - φ, μ (idx ≥ 1): mass-lumped L² projection preserving ∫φ dV and ∫μ dV
 */
template<class TypeTag>
class HeleShawTwoPGridDataTransfer
: public GridDataTransfer<GetPropType<TypeTag, Properties::Grid>>
{
    using Grid         = GetPropType<TypeTag, Properties::Grid>;
    using ParentType   = GridDataTransfer<Grid>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using Problem      = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Extrusion    = Extrusion_t<GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector   = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename Grid::template Codim<0>::Entity;
    using ElementSolution = std::decay_t<decltype(elementSolution(
        std::declval<Element>(), std::declval<SolutionVector>(),
        std::declval<GridGeometry>()))>;

    static_assert(GridGeometry::discMethod == DiscretizationMethods::box,
        "HeleShawTwoPGridDataTransfer requires the Box discretization");

    static constexpr int massProjBegin = 1; // φ, μ get mass-conservative projection
    static constexpr int massProjEnd   = PrimaryVariables::size();

    struct AdaptedValues
    {
        AdaptedValues() : associatedMass(0.0) {}
        ElementSolution u;
        PrimaryVariables associatedMass;
    };
    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

public:
    HeleShawTwoPGridDataTransfer(
        std::shared_ptr<const Problem>       problem,
        std::shared_ptr<GridGeometry>        gridGeometry,
        std::shared_ptr<const GridVariables> gridVariables,
        SolutionVector&                      sol)
    : ParentType()
    , problem_(problem), gg_(gridGeometry)
    , gv_(gridVariables), sol_(sol)
    , map_(gridGeometry->gridView().grid(), 0)
    {}

    void store(const Grid& grid) override
    {
        map_.resize();
        auto fvGeometry = localView(*gg_);

        for (auto level = grid.maxLevel(); level >= 0; --level)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                auto& av = map_[element];
                av.u = ElementSolution(element, sol_, *gg_);

                if (element.isLeaf())
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables vv;
                        vv.update(av.u, *problem_, element, scv);
                        const auto vol = Extrusion::volume(fvGeometry, scv);
                        for (int j = massProjBegin; j < massProjEnd; ++j)
                            av.associatedMass[j] += vol * vv.priVar(j);
                    }
                }

                if (element.level() > 0)
                    map_[element.father()].associatedMass += av.associatedMass;
            }
        }
    }

    void reconstruct(const Grid& grid) override
    {
        gg_->update(grid.leafGridView());
        map_.resize();
        sol_.resize(gg_->numDofs());

        std::vector<PrimaryVariables> massCoeff(gg_->numDofs(), PrimaryVariables(0.0));
        std::vector<PrimaryVariables> assocMass(gg_->numDofs(), PrimaryVariables(0.0));

        auto fvGeometry = localView(*gg_);
        for (const auto& element : elements(grid.leafGridView(), Dune::Partitions::interior))
        {
            fvGeometry.bindElement(element);
            const auto src     = sourceElement_(element);
            const auto& stored = map_[src];
            const auto srcGeom = src.geometry();
            const auto srcVol  = volume(srcGeom, Extrusion{});

            ElementSolution elemSol(element, sol_, *gg_);
            if (element.isNew())
                for (const auto& scv : scvs(fvGeometry))
                    elemSol[scv.localDofIndex()] = evalSolution(
                        src, srcGeom, *gg_, stored.u, scv.dofPosition());
            else
                elemSol = stored.u;

            for (const auto& scv : scvs(fvGeometry))
            {
                sol_[scv.dofIndex()] = elemSol[scv.localDofIndex()];
                const auto vol = Extrusion::volume(fvGeometry, scv);
                for (int j = massProjBegin; j < massProjEnd; ++j)
                {
                    massCoeff[scv.dofIndex()][j] += vol;
                    assocMass[scv.dofIndex()][j] += vol / srcVol * stored.associatedMass[j];
                }
            }
        }

        for (std::size_t i = 0; i < gg_->numDofs(); ++i)
            for (int j = massProjBegin; j < massProjEnd; ++j)
                if (massCoeff[i][j] > 0.0)
                    sol_[i][j] = assocMass[i][j] / massCoeff[i][j];

        map_.resize(typename PersistentContainer::Value());
        map_.shrinkToFit();
        map_.fill(typename PersistentContainer::Value());
    }

private:
    Element sourceElement_(const Element& e) const
    {
        if (!e.isNew()) return e;
        if (!e.hasFather())
            DUNE_THROW(Dune::InvalidStateException, "New element has no father");
        auto ancestor = e.father();
        while (ancestor.isNew() && ancestor.level() > 0)
            ancestor = ancestor.father();
        return ancestor;
    }

    std::shared_ptr<const Problem>       problem_;
    std::shared_ptr<GridGeometry>        gg_;
    std::shared_ptr<const GridVariables> gv_;
    SolutionVector&                      sol_;
    PersistentContainer                  map_;
};

} // end namespace Dumux

#endif
