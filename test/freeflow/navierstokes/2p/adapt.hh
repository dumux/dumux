// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
// Loosely based on Dumux::TwoPGridAdaptIndicator published under the above license.
// SPDX-FileCopyrightText: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Adaptive grid refinement for the two-phase Cahn-Hilliard / Navier-Stokes
// system. Works for both the plain mass model (NavierStokesMassTwoP, primary
// variables [p, φ, μ]) and its vapor extension (NavierStokesMassTwoPVapor,
// primary variables [p, φ, μ, c_v]).
//
//   - φ (phaseFieldIdx)         : drives the refinement indicator and is the
//                                 conserved quantity of the Cahn-Hilliard equation.
//   - μ (chemicalPotentialIdx)  : Lagrange-multiplier-like; reprojected for smoothness.
//   - c_v (vaporIdx, if present): transported mass concentration; mass-conservative
//                                 reprojection preserves ∫ c_v dV under AMR.
//   - p (pressureIdx)           : not extensive; receives plain interpolation.
//
#ifndef DUMUX_TEST_FREEFLOW_2P_ADAPT_HH
#define DUMUX_TEST_FREEFLOW_2P_ADAPT_HH

#include <limits>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

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
 * \brief Refinement indicator based on jumps of the phase field across faces.
 *
 * Returns +1 to refine, -1 to coarsen, 0 to keep. φ ∈ [-1, 1], so absolute and
 * relative tolerances coincide and we use globalDelta = 1 as the reference scale.
 */
template<class TypeTag>
class TwoPhaseCahnHilliardGridAdaptIndicator
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int indicatorIdx = Indices::phaseFieldIdx;

public:
    TwoPhaseCahnHilliardGridAdaptIndicator(std::shared_ptr<const GridGeometry> gridGeometry,
                                           const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , maxLocalConcentrationDelta_(gridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    , minElementSize_(getParamFromGroup<Scalar>(paramGroup, "Adaptive.MinElementSize", std::numeric_limits<Scalar>::max()))
    {}

    void calculate(const SolutionVector& sol, Scalar refineTol = 0.05, Scalar coarsenTol = 0.001)
    {
        refineBound_ = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        maxLocalConcentrationDelta_.assign(gridGeometry_->gridView().size(0), 0.0);

        if (minLevel_ >= maxLevel_)
            return;

        if (coarsenTol > refineTol)
            DUNE_THROW(Dune::InvalidStateException, "Refine tolerance must be higher than coarsen tolerance");

        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            const auto eIdxI = gridGeometry_->elementMapper().index(element);
            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, *gridGeometry_);
            const Scalar cI = evalSolution(element, geometry, *gridGeometry_, elemSol, geometry.center())[indicatorIdx];
            using std::max, std::abs;

            for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
            {
                if (!intersection.neighbor())
                    continue;

                const auto& outside = intersection.outside();
                const auto eIdxJ = gridGeometry_->elementMapper().index(outside);

                // visit each face only once (finer side wins at hanging nodes)
                if (element.level() > outside.level()
                    || (element.level() == outside.level() && eIdxI < eIdxJ))
                {
                    const auto outsideGeometry = outside.geometry();
                    const auto elemSolJ = elementSolution(outside, sol, *gridGeometry_);
                    const Scalar cJ = evalSolution(outside, outsideGeometry, *gridGeometry_, elemSolJ, outsideGeometry.center())[indicatorIdx];
                    const Scalar localdelta = abs(cI - cJ);
                    maxLocalConcentrationDelta_[eIdxI] = max(maxLocalConcentrationDelta_[eIdxI], localdelta);
                    maxLocalConcentrationDelta_[eIdxJ] = max(maxLocalConcentrationDelta_[eIdxJ], localdelta);
                }
            }
        }

        const Scalar globalDelta = 1.0; // φ ∈ [-1, 1]
        refineBound_ = refineTol*globalDelta;
        coarsenBound_ = coarsenTol*globalDelta;

        for (const auto& element : elements(gridGeometry_->gridView(), Dune::Partitions::interior))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);
    }

    int operator() (const Element& element) const
    {
        const auto eIdx = gridGeometry_->elementMapper().index(element);
        // Enforce a uniform minimum refinement level everywhere, not just the interface band:
        // the φ-gradient indicator otherwise coarsens the bulk down to the macro mesh, leaving
        // the (PQ1Bubble) velocity field in the bubble wake / return flow badly under-resolved.
        if (element.level() < int(minLevel_))
            return 1;
        if (element.hasFather() && element.level() > int(minLevel_)
            && maxLocalConcentrationDelta_[eIdx] < coarsenBound_)
            return -1;
        if (element.level() < maxLevel_
            && size_(element) > minElementSize_
            && maxLocalConcentrationDelta_[eIdx] > refineBound_)
            return 1;
        return 0;
    }

private:
    // Diameter-based 2:1 closure (works for non-uniformly sized grids).
    bool checkNeighborsRefine_(const Element& element, std::size_t level = 1)
    {
        for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            if (!intersection.neighbor())
                continue;

            const auto& outside = intersection.outside();
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            const auto outsideSize = size_(outside);
            if (outsideSize > minElementSize_
                && outsideSize > 2.0*size_(element)
                && outside.level() < maxLevel_)
            {
                maxLocalConcentrationDelta_[gridGeometry_->elementMapper().index(outside)]
                    = std::numeric_limits<Scalar>::max();
                if (level < maxLevel_)
                    checkNeighborsRefine_(outside, ++level);
            }
        }
        return true;
    }

    Scalar size_(const Element& element) const
    { return diameter(element.geometry()); }

    std::shared_ptr<const GridGeometry> gridGeometry_;

    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector<Scalar> maxLocalConcentrationDelta_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
    Scalar minElementSize_;
};


/*!
 * \brief Mass-side data transfer for the box discretization of the
 *        two-phase Cahn-Hilliard / Navier-Stokes mass model with vapor.
 *
 * Strategy per primary variable:
 *   - j = 0  (pressure)  : plain interpolation via evalSolution.
 *   - j ≥ 1  (φ, μ, c_v) : mass-lumped L² projection so that the per-element
 *                          integral ∫_E u_j dV is preserved across adaptation.
 *
 * For the box scheme, vertex DOFs persist across levels, so the element
 * solution can be reconstructed at any non-leaf element directly from the
 * leaf-grid sol_ vector.
 */
template<class TypeTag>
class TwoPhaseCahnHilliardMassGridDataTransfer
: public GridDataTransfer<GetPropType<TypeTag, Properties::Grid>>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using ParentType = GridDataTransfer<Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename Grid::template Codim<0>::Entity;
    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element>(),
                                                                  std::declval<SolutionVector>(),
                                                                  std::declval<GridGeometry>()))>;

    struct AdaptedValues
    {
        AdaptedValues() : associatedMass(0.0) {}
        ElementSolution u;
        PrimaryVariables associatedMass; // ∫_E u_j dV for j ≥ 1
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    static_assert(isBox, "TwoPhaseCahnHilliardMassGridDataTransfer requires the box discretization.");

    // Mass-projected component range: everything except pressure (j=0).
    static constexpr int massProjBegin = 1;
    static constexpr int massProjEnd   = PrimaryVariables::size();

public:
    TwoPhaseCahnHilliardMassGridDataTransfer(
        std::shared_ptr<const Problem> problem,
        std::shared_ptr<GridGeometry> gridGeometry,
        std::shared_ptr<const GridVariables> gridVariables,
        SolutionVector& sol
    ) : ParentType()
    , problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , adaptionMap_(gridGeometry->gridView().grid(), 0)
    {}

    void store(const Grid& grid) override
    {
        adaptionMap_.resize();
        auto fvGeometry = localView(*gridGeometry_);

        // Iterate from the finest level down so child masses are summed into
        // their fathers before we visit them.
        for (auto level = grid.maxLevel(); level >= 0; --level)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                auto& adaptedValues = adaptionMap_[element];

                // Box: vertex DOFs live on every level, so this works for both
                // leaf and non-leaf elements (leaf vertex values via global sol_).
                adaptedValues.u = ElementSolution(element, sol_, *gridGeometry_);

                if (element.isLeaf())
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u, *problem_, element, scv);
                        const auto scvVolume = Extrusion::volume(fvGeometry, scv);
                        for (int j = massProjBegin; j < massProjEnd; ++j)
                            adaptedValues.associatedMass[j] += scvVolume * volVars.priVar(j);
                    }
                }

                // Carry the mass up so a coarsened-up father contains its
                // (former) children's integrated mass.
                if (element.level() > 0)
                {
                    auto& fatherValues = adaptionMap_[element.father()];
                    if (&adaptedValues != &fatherValues)
                        fatherValues.associatedMass += adaptedValues.associatedMass;
                }
            }
        }
    }

    void reconstruct(const Grid& grid) override
    {
        gridGeometry_->update(grid.leafGridView());
        reconstruct_();
    }

private:
    void reconstruct_()
    {
        adaptionMap_.resize();
        sol_.resize(gridGeometry_->numDofs());

        // Mass-lumped L² projection accumulators, per DOF and component:
        //   sol_[i][j] = (Σ scvVolume * <source u_j>) / (Σ scvVolume)
        // NOTE: this global lumped redistribution both conserves ∫φ AND smooths; the
        // smoothing is essential to tame the per-step re-meshing perturbation (plain
        // interpolation without it gives ~5x the seam pressure noise). A consistent-mass
        // (not lumped) L² projection would be identity on FE fields and conservative, but
        // is heavier; the lumped version is the best simple trade-off.
        std::vector<PrimaryVariables> massCoeff(gridGeometry_->numDofs(), PrimaryVariables(0.0));
        std::vector<PrimaryVariables> associatedMass(gridGeometry_->numDofs(), PrimaryVariables(0.0));

        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView().grid().leafGridView(), Dune::Partitions::interior))
        {
            fvGeometry.bindElement(element);

            // Source on the *old* grid that holds our stored data:
            //   - surviving leaf      → the element itself (elementVolume ≡ sourceVolume)
            //   - new (refined) leaf  → youngest ancestor that existed before adaptation
            const auto source = sourceElement_(element);
            const auto& storedValues = adaptionMap_[source];
            const auto sourceGeometry = source.geometry();
            const auto sourceVolume = volume(sourceGeometry, Extrusion{});

            // elemSol carrying the plain-interpolation values for this leaf.
            ElementSolution elemSol(element, sol_, *gridGeometry_);
            if (element.isNew())
            {
                for (const auto& scv : scvs(fvGeometry))
                    elemSol[scv.localDofIndex()] = evalSolution(
                        source, sourceGeometry, *gridGeometry_, storedValues.u, scv.dofPosition()
                    );
            }
            else
            {
                elemSol = storedValues.u;
            }

            for (const auto& scv : scvs(fvGeometry))
            {
                // Pressure (and all other components) get plain interpolation here;
                // mass-projected components are overwritten below.
                sol_[scv.dofIndex()] = elemSol[scv.localDofIndex()];

                const auto scvVolume = Extrusion::volume(fvGeometry, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                for (int j = massProjBegin; j < massProjEnd; ++j)
                {
                    massCoeff[dofIdxGlobal][j]      += scvVolume;
                    associatedMass[dofIdxGlobal][j] += scvVolume / sourceVolume * storedValues.associatedMass[j];
                }
            }
        }

        // The lumped L2 projection conserves the per-component integral but applies an
        // ~h^2/4 Laplacian smoothing each remesh, AND it maps phi/mu DIFFERENTLY from the pressure
        // (which gets plain interpolation above). That inconsistency leaves the modified pressure
        // p* = p - mu*phi out of sync with the smoothed mu*phi right after a remesh and feeds
        // single-node interfacial pressure spikes. Default OFF: map all mass variables the same way
        // (plain interpolation = the EXACT, non-smoothing, integral-preserving consistent projection
        // for nested refinement). Verified: mass drift stays ~0.0008% over t=3 (lumping not needed
        // for conservation), the remesh flicker halves, and the QoI are unchanged.
        static const bool useLumped = getParamFromGroup<bool>(
            problem_->paramGroup(), "Adaptive.LumpedMassProjection", false);
        if (useLumped)
            for (std::size_t i = 0; i < gridGeometry_->numDofs(); ++i)
                for (int j = massProjBegin; j < massProjEnd; ++j)
                    if (massCoeff[i][j] > 0.0)
                        sol_[i][j] = associatedMass[i][j] / massCoeff[i][j];

        adaptionMap_.resize(typename PersistentContainer::Value());
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill(typename PersistentContainer::Value());
    }

    Element sourceElement_(const Element& element) const
    {
        if (!element.isNew())
            return element;

        if (!element.hasFather())
            DUNE_THROW(Dune::InvalidStateException,
                "New element does not have a father; cannot reconstruct solution.");

        auto ancestor = element.father();
        while (ancestor.isNew() && ancestor.level() > 0)
            ancestor = ancestor.father();
        return ancestor;
    }

    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<GridGeometry> gridGeometry_;
    std::shared_ptr<const GridVariables> gridVariables_;
    SolutionVector& sol_;
    PersistentContainer adaptionMap_;
};


/*!
 * \brief Momentum-side (velocity) data transfer for control-volume FE
 *        discretizations with element-local DOFs (e.g. pq1bubble).
 *
 * Velocity components are not extensive conserved quantities at DOF level —
 * the momentum balance is enforced by the PDE — so plain basis-function
 * interpolation is the correct transfer here.
 *
 * Subtlety: pq1bubble has both vertex DOFs (codim=dim, persistent across
 * levels) and a centroid bubble DOF (codim=0, per-element). For non-leaf
 * elements the centroid DOF has no entry in the leaf-grid dof mapper, so
 * the standard `ElementSolution(element, sol_, gg)` constructor crashes.
 * For non-leaf elements we therefore construct the element solution
 * manually: vertex DOFs are read from sol_ via the (persistent) vertex
 * mapper, and centroid DOFs are seeded with zero — the next Newton solve
 * recomputes the bubble.
 */
template<class TypeTag>
class TwoPhaseCahnHilliardVelocityGridDataTransfer
: public GridDataTransfer<GetPropType<TypeTag, Properties::Grid>>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using ParentType = GridDataTransfer<Grid>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename Grid::template Codim<0>::Entity;
    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element>(),
                                                                  std::declval<SolutionVector>(),
                                                                  std::declval<GridGeometry>()))>;

    // associatedMomentum = ∫_E u dV (per velocity component), carried conservatively across
    // adaptation so the transferred velocity (which becomes xOld and enters the momentum
    // time-derivative term) preserves momentum - otherwise each remesh injects a spurious
    // momentum change that the next solve registers as a jump in the rise velocity.
    struct StoredValue
    {
        StoredValue() : associatedMomentum(0.0) {}
        ElementSolution u;
        PrimaryVariables associatedMomentum;
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, StoredValue>;

    static constexpr int dim = Grid::dimension;

public:
    TwoPhaseCahnHilliardVelocityGridDataTransfer(
        std::shared_ptr<const Problem> problem,
        std::shared_ptr<GridGeometry> gridGeometry,
        std::shared_ptr<const GridVariables> gridVariables,
        SolutionVector& sol
    ) : ParentType()
    , problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , adaptionMap_(gridGeometry->gridView().grid(), 0)
    {}

    void store(const Grid& grid) override
    {
        adaptionMap_.resize();
        auto fvGeometry = localView(*gridGeometry_);
        // finest level first so child momenta are summed into fathers before we visit them
        for (auto level = grid.maxLevel(); level >= 0; --level)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                auto& stored = adaptionMap_[element];
                if (element.isLeaf())
                {
                    stored.u = ElementSolution(element, sol_, *gridGeometry_);
                    fvGeometry.bindElement(element);
                    // conserve momentum on the VERTEX dofs only; the element-local pq1bubble
                    // centroid dof is left to interpolation (collapsing it to the element mean
                    // destroys inf-sup stability) and is recomputed by the next solve.
                    const auto numCorners = static_cast<std::size_t>(element.subEntities(dim));
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        if (scv.localDofIndex() >= numCorners) continue;
                        const auto vol = Extrusion::volume(fvGeometry, scv);
                        for (int d = 0; d < PrimaryVariables::size(); ++d)
                            stored.associatedMomentum[d] += stored.u[scv.localDofIndex()][d] * vol;
                    }
                }
                else
                    storeNonLeafElementSolution_(element, stored.u, grid.maxLevel());

                // carry the momentum up so a coarsened father holds its children's momentum
                if (element.level() > 0)
                {
                    auto& fatherValues = adaptionMap_[element.father()];
                    if (&stored != &fatherValues)
                        fatherValues.associatedMomentum += stored.associatedMomentum;
                }
            }
        }
    }

    void reconstruct(const Grid& grid) override
    {
        gridGeometry_->update(grid.leafGridView());
        adaptionMap_.resize();
        sol_.resize(gridGeometry_->numDofs());

        // Momentum-conserving transfer: redistribute the stored element momentum int_E u dV
        // via a lumped L2 projection so the global momentum sum_i V_i u_i is preserved across
        // the remesh (the velocity otherwise is only interpolated, so a remesh injects a small
        // spurious momentum change). NOTE default OFF: tested at N=5 it conserves momentum but
        // the lumped redistribution also smooths the SURVIVING-element vertex velocities, which
        // damps the rise-velocity peak (0.248 -> 0.218) and worsens max|v| roughness (the
        // smoothed vertices mismatch the interpolated, un-smoothed pq1bubble centroid dof). A
        // net win needs an identity-on-surviving consistent projection (not lumped) - left as
        // future work. The residual V_c roughness is the inherent mesh-dependence of the
        // converged solution at each remesh, not a fixable transfer error.
        static const bool conserve = getParamFromGroup<bool>(
            problem_->paramGroup(), "Adaptive.ConserveMomentum", false);
        std::vector<PrimaryVariables> momCoeff, momentum;
        if (conserve)
        {
            momCoeff.assign(gridGeometry_->numDofs(), PrimaryVariables(0.0));
            momentum.assign(gridGeometry_->numDofs(), PrimaryVariables(0.0));
        }

        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView().grid().leafGridView(), Dune::Partitions::interior))
        {
            fvGeometry.bindElement(element);

            // the element that holds the stored source data (self if surviving, else ancestor).
            // NOTE: iterate localDofs(), not scvs() -- for a plain box/pq1bubble scheme these are
            // the same index set (scvs() already includes the element-local bubble dof there), but
            // for a hybrid PQ2 scheme scvs() is vertex-only and would silently skip the P2 edge
            // dofs, leaving them at their (zero) resize default on every remesh.
            auto source = element;
            if (!element.isNew())
            {
                const auto& stored = adaptionMap_[element];
                for (const auto& localDof : localDofs(fvGeometry))
                    sol_[localDof.dofIndex()] = stored.u[localDof.index()];
            }
            else
            {
                if (!element.hasFather())
                    DUNE_THROW(Dune::InvalidStateException,
                        "New element does not have a father; cannot reconstruct solution.");

                source = element.father();
                while (source.isNew() && source.level() > 0)
                    source = source.father();

                const auto& stored = adaptionMap_[source];
                const auto sourceGeometry = source.geometry();
                for (const auto& localDof : localDofs(fvGeometry))
                    sol_[localDof.dofIndex()] = evalSolution(
                        source, sourceGeometry, *gridGeometry_, stored.u, ipData(fvGeometry, localDof).global()
                    );
            }

            if (conserve)
            {
                const auto& storedSrc = adaptionMap_[source];
                const auto sourceVolume = volume(source.geometry(), Extrusion{});
                const auto numCorners = static_cast<std::size_t>(element.subEntities(dim));
                for (const auto& scv : scvs(fvGeometry))
                {
                    if (scv.localDofIndex() >= numCorners) continue; // skip bubble dof (keep interpolated)
                    const auto vol = Extrusion::volume(fvGeometry, scv);
                    const auto dofIdx = scv.dofIndex();
                    for (int d = 0; d < PrimaryVariables::size(); ++d)
                    {
                        momCoeff[dofIdx][d] += vol;
                        momentum[dofIdx][d] += vol / sourceVolume * storedSrc.associatedMomentum[d];
                    }
                }
            }
        }

        if (conserve)
            for (std::size_t i = 0; i < gridGeometry_->numDofs(); ++i)
                for (int d = 0; d < PrimaryVariables::size(); ++d)
                    if (momCoeff[i][d] > 0.0)
                        sol_[i][d] = momentum[i][d] / momCoeff[i][d];

        adaptionMap_.resize(typename PersistentContainer::Value());
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill(typename PersistentContainer::Value());
    }

private:
    // Construct an element solution for a non-leaf element when the discretization
    // has element-local DOFs (e.g. the pq1bubble centroid bubble).
    void storeNonLeafElementSolution_(const Element& element, ElementSolution& u, int maxLevel)
    {
        // Initialize size via any leaf descendant (sibling/child elements share
        // the same element topology, hence the same number of local DOFs).
        for (auto it = element.hbegin(maxLevel); it != element.hend(maxLevel); ++it)
        {
            if (it->isLeaf())
            {
                u = ElementSolution(*it, sol_, *gridGeometry_);
                break;
            }
        }

        const auto& lfe = gridGeometry_->feCache().get(element.type());
        const auto& localCoeff = lfe.localCoefficients();
        const auto& localBasis = lfe.localBasis();

        // (1) Vertex DOFs: persistent nodal values.
        for (int i = 0; i < localCoeff.size(); ++i)
            if (localCoeff.localKey(i).codim() == dim)
            {
                const auto vIdx = gridGeometry_->dofMapper().subIndex(
                    element, localCoeff.localKey(i).subEntity(), dim);
                u[i] = sol_[vIdx];
            }

        // (2) Non-vertex DOFs (the pq1bubble centroid, or P2's edge-midpoint DOFs): solve the
        // LOCAL L2 projection problem for this whole dof block JOINTLY, given the (fixed)
        // vertex DOFs. We seek the coefficients u_a minimising the L2 error
        //     || u_src - (u_vert + Σ_a u_a B_a) ||^2_{L2(E)}
        // i.e. the small dense system  M u = r  with
        //     M_ab = ∫_E B_a B_b dV                  (assembled directly over E; independent
        //                                              of the leaf children below)
        //     r_a  = ∫_E B_a (u_src - u_vert) dV      (u_src lives on the leaf children, so
        //                                              integrated piecewise over them, which
        //                                              tile E and still hold the fine solution
        //                                              at store time)
        // For exactly one non-vertex dof (pq1bubble) M is 1x1 and this reduces to the previous
        // per-dof divide. For P2's THREE coupled edge dofs the off-diagonal M_ab (a != b) terms
        // are NOT negligible -- solving each dof independently (as if the others were absent,
        // i.e. using only the diagonal) is a real cross-coupling error. (Zeroing the coefficients
        // instead would collapse the field on coarsening -> jumpy cell-to-cell values + spurious
        // storage-term kick next solve.)
        static const bool nodalSeed = getParamFromGroup<bool>(
            problem_->paramGroup(), "Adaptive.BubbleNodalSeed", true);
        using ctype = typename Element::Geometry::ctype;
        using ShapeVal = typename std::decay_t<decltype(localBasis)>::Traits::RangeType;
        const auto fgeo = element.geometry();

        std::vector<int> nonVertexIdx;
        for (int i = 0; i < localCoeff.size(); ++i)
            if (localCoeff.localKey(i).codim() != dim)
                nonVertexIdx.push_back(i);

        if (nonVertexIdx.empty())
            return;

        if (!nodalSeed) // old behaviour for A/B
        {
            for (const auto i : nonVertexIdx)
                u[i] = PrimaryVariables(0.0);
            return;
        }

        const auto numNV = nonVertexIdx.size();
        Dune::DynamicMatrix<ctype> M(numNV, numNV, 0.0);
        std::vector<PrimaryVariables> r(numNV, PrimaryVariables(0.0));

        // Mass matrix restricted to the non-vertex-dof block, assembled directly over the
        // father element E (does not depend on the leaf children).
        std::vector<ShapeVal> shpF;
        const auto& massQuad = Dune::QuadratureRules<ctype, dim>::rule(element.type(), 2*localBasis.order());
        for (const auto& qp : massQuad)
        {
            localBasis.evaluateFunction(qp.position(), shpF);
            const auto w = qp.weight() * fgeo.integrationElement(qp.position());
            for (std::size_t a = 0; a < numNV; ++a)
                for (std::size_t c = 0; c < numNV; ++c)
                    M[a][c] += w * shpF[nonVertexIdx[a]][0] * shpF[nonVertexIdx[c]][0];
        }

        // RHS: u_src (fine children field) minus the (fixed) vertex part, projected onto each
        // non-vertex basis function B_a.
        std::vector<ShapeVal> shp;
        const auto& rhsQuad = Dune::QuadratureRules<ctype, dim>::rule(element.type(), 4);
        for (auto it = element.hbegin(maxLevel); it != element.hend(maxLevel); ++it)
        {
            if (!it->isLeaf()) continue;
            const auto cgeo = it->geometry();
            const auto csol = ElementSolution(*it, sol_, *gridGeometry_);
            for (const auto& qp : rhsQuad)
            {
                const auto global = cgeo.global(qp.position());
                localBasis.evaluateFunction(fgeo.local(global), shp);

                PrimaryVariables uVert(0.0); // father vertex part Σ_v u_v Ñ_v at this point
                for (int k = 0; k < localCoeff.size(); ++k)
                    if (localCoeff.localKey(k).codim() == dim)
                        uVert.axpy(shp[k][0], u[k]);

                const auto uSrc = evalSolution(*it, cgeo, *gridGeometry_, csol, global);
                const auto w = qp.weight() * cgeo.integrationElement(qp.position());
                auto diff = uSrc; diff -= uVert;

                for (std::size_t a = 0; a < numNV; ++a)
                    r[a].axpy(w * shp[nonVertexIdx[a]][0], diff);
            }
        }

        // Solve the small dense system M u = r, once per PrimaryVariables component.
        Dune::DynamicVector<ctype> rhsComp(numNV), solComp(numNV);
        for (int d = 0; d < PrimaryVariables::size(); ++d)
        {
            for (std::size_t a = 0; a < numNV; ++a)
                rhsComp[a] = r[a][d];

            auto Mcopy = M; // solve() factorizes in place
            Mcopy.solve(solComp, rhsComp);

            for (std::size_t a = 0; a < numNV; ++a)
                u[nonVertexIdx[a]][d] = solComp[a];
        }
    }

    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<GridGeometry> gridGeometry_;
    std::shared_ptr<const GridVariables> gridVariables_;
    SolutionVector& sol_;
    PersistentContainer adaptionMap_;
};


/*!
 * \brief Composite data transfer that drives both sub-domain transfers.
 */
template<class Grid, class GDT1, class GDT2>
class TwoDomainOneGridDataTransfer
: public GridDataTransfer<Grid>
{
    using ParentType = GridDataTransfer<Grid>;
public:
    TwoDomainOneGridDataTransfer(std::shared_ptr<GDT1> gdt1,
                                 std::shared_ptr<GDT2> gdt2)
    : ParentType()
    , gdt1_(gdt1)
    , gdt2_(gdt2)
    {}

    void store(const Grid& grid) override
    {
        gdt1_->store(grid);
        gdt2_->store(grid);
    }

    void reconstruct(const Grid& grid) override
    {
        gdt1_->reconstruct(grid);
        gdt2_->reconstruct(grid);
    }

private:
    std::shared_ptr<GDT1> gdt1_;
    std::shared_ptr<GDT2> gdt2_;
};

} // end namespace Dumux

#endif
