// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Coupling manager for the hydrogel actuator (gripper) model.
 */
#ifndef DUMUX_GRIPPER_MODEL_COUPLINGMANAGER_HH
#define DUMUX_GRIPPER_MODEL_COUPLINGMANAGER_HH

#include <cmath>
#include <memory>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \brief Coupling manager for the gripper poromechanical problem.
 *
 * Both subdomains share the same grid.
 * - Subdomain 0 (momentum, PQ1Bubble) → subdomain 1 (pressure, Box) coupling:
 *   stencil = vertex DOF indices of the element in the pressure grid geometry.
 * - Subdomain 1 (pressure, Box) → subdomain 0 (momentum, PQ1Bubble) coupling:
 *   stencil = all DOF indices of the element in the momentum grid geometry
 *             (vertices + bubble DOFs).
 *
 * Key services:
 * - `pressureAtPoint(element, ip)`: evaluates \f$p\f$ at a global position
 *   inside the element using Box P1 shape functions and the stored pressure solution.
 *   Used by the momentum local residual.
 * - `deformationGradientAtPoint(element, ip)`: evaluates \f$\mathbf{F}\f$ at a
 *   global position using PQ1Bubble shape functions and the stored displacement solution.
 *   Used by the pressure local residual.
 * - `prevJacobian(dofIdx)`: returns \f$J^n\f$ stored per Box DOF.
 *   Updated via `savePrevTimeStepJacobian()` from main.cc after each converged step.
 */
template<class MDTraits>
class GripperCouplingManager : public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using Indices
        = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>::VolumeVariables::Indices;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    static constexpr std::size_t momentumIdx_ = 0;
    static constexpr std::size_t massIdx_ = 1;

    using MomFVGeometry  = typename GridGeometry<momentumIdx_>::LocalView;
    using PresFVGeometry = typename GridGeometry<massIdx_>::LocalView;

    using GlobalPosition = typename Element<0>::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView<0>::dimension, GridView<0>::dimension>;
    static constexpr int dim = GridView<0>::dimension;

    static constexpr auto pressureIdx = Indices<massIdx_>::pressureIdx;

public:
    static constexpr auto momentumIdx = Dune::index_constant<momentumIdx_>();
    static constexpr auto massIdx = Dune::index_constant<massIdx_>();

    using CouplingStencilType = std::vector<std::size_t>;

    GripperCouplingManager(std::shared_ptr<GridGeometry<momentumIdx_>> momentumGG,
                           std::shared_ptr<GridGeometry<massIdx_>> massGG)
    : ParentType()
    {
        // Initialise prevJ_ to 1 (reference configuration J = det I = 1).
        prevJ_.assign(massGG->numDofs(), Scalar(1.0));

        buildCouplingStencils_(*momentumGG, *massGG);
    }

    /*!
     * \brief Initialize the coupling manager.
     *
     * Stores grid geometry pointers and pre-builds coupling stencils for both directions.
     */
    void init(std::shared_ptr<Problem<momentumIdx_>> momProblem,
              std::shared_ptr<Problem<massIdx_>> presProblem,
              const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        this->setSubProblems(std::make_tuple(momProblem, presProblem));
    }

    template<std::size_t i, std::size_t j>
    const CouplingStencilType& couplingStencil(Dune::index_constant<i>,
                                               const Element<j>& element,
                                               Dune::index_constant<j>) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");
        if constexpr (i == momentumIdx_ && j == massIdx_)
        {
            const auto eIdx = this->problem(massIdx).gridGeometry().elementMapper().index(element);
            return stencilsMomToPres_[eIdx];
        }
        else
        {
            const auto eIdx = this->problem(momentumIdx).gridGeometry().elementMapper().index(element);
            return stencilsPresToMom_[eIdx];
        }
    }

    /*!
     * \brief Evaluate the fluid pore pressure at a global position inside an element.
     *
     * Called from the momentum local residual.
     */
    Scalar pressureAtPoint(const typename GridGeometry<momentumIdx_>::LocalView& fvGeometry,
                           const GlobalPosition& ip) const
    {
        const auto& gg = this->problem(massIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(massIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            ip
        )[pressureIdx];
    }

    /*!
     * \brief Evaluate the solid bulk pressure \f$p_s\f$ at a global position.
     *
     * Called from the momentum local residual. Reads index 1 from the mass solution
     * (P1 discretisation, same basis as fluid pressure).
     */
    Scalar solidBulkPressureAtPoint(const typename GridGeometry<momentumIdx_>::LocalView& fvGeometry,
                                    const GlobalPosition& ip) const
    {
        const auto& gg = this->problem(massIdx).gridGeometry();
        const auto elemSol = elementSolution(fvGeometry.element(), curSol(massIdx), gg);
        return evalSolution(
            fvGeometry.element(),
            fvGeometry.element().geometry(),
            gg, elemSol,
            ip
        )[Indices<massIdx_>::solidBulkPressureIdx];
    }

    /*!
     * \brief Evaluate the deformation gradient \f$\mathbf{F}\f$ at a global position
     * inside an element using the PQ1Bubble basis.
     *
     * Called from the pressure local residual.
     */
    Tensor deformationGradientAtPoint(const typename GridGeometry<massIdx_>::LocalView& fvGeometry,
                                      const GlobalPosition& ip) const
    {
        const auto& gg = this->problem(momentumIdx).gridGeometry();
        const auto fvGeomMom = localView(gg).bindElement(fvGeometry.element());
        return deformationGradientAtPoint_(fvGeometry.element(), ip, fvGeomMom);
    }

    /*!
     * \brief Return the stored \f$J^n\f$ (det F from previous converged time step)
     * at a given pressure-subdomain DOF index.
     */
    Scalar prevJacobian(std::size_t dofIdx) const
    {
        return prevJ_[dofIdx];
    }

    /*!
     * \brief Compute and store \f$J^n\f$ at every Box DOF from the current momentum
     * solution. Call this from main.cc after each converged time step.
     */
    void savePrevTimeStepJacobian()
    {
        const auto& massGG = this->problem(massIdx).gridGeometry();
        const auto& momGG = this->problem(momentumIdx).gridGeometry();

        prevJ_.assign(massGG.numDofs(), Scalar(1.0));
        std::vector<bool> visited(massGG.numDofs(), false);

        auto fvGeomMom  = localView(momGG);
        auto fvGeomPres = localView(massGG);

        for (const auto& element : elements(massGG.gridView()))
        {
            fvGeomMom.bindElement(element);
            fvGeomPres.bindElement(element);

            for (const auto& scv : scvs(fvGeomPres))
            {
                if (visited[scv.dofIndex()]) continue;
                visited[scv.dofIndex()] = true;
                prevJ_[scv.dofIndex()] =
                    deformationGradientAtPoint_(element, scv.dofPosition(), fvGeomMom)
                    .determinant();
            }
        }
    }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    auto& curSol(Dune::index_constant<i> domainIdx)
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    const auto& curSol(Dune::index_constant<i> domainIdx) const
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    { elementSets_ = computeColoring(this->problem(momentumIdx).gridGeometry()).sets; }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainId the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (elementSets_.empty())
            DUNE_THROW(Dune::InvalidStateException,
                "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix or modifying grid variable caches
        // each color can be assembled using multiple threads
        const auto& grid = this->problem(domainId).gridGeometry().gridView().grid();
        for (const auto& elements : elementSets_)
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t n)
            {
                const auto element = grid.entity(elements[n]);
                assembleElement(element);
            });
        }
    }

private:
    Tensor deformationGradientAtPoint_(const Element<massIdx_>& element,
                                       const GlobalPosition& ip,
                                       const MomFVGeometry& fvGeomMom) const
    {
        const auto& geom = element.geometry();
        const auto ipLocal = geom.local(ip);

        const auto& localBasis = fvGeomMom.feLocalBasis();
        using ShapeJacobian = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
        std::vector<ShapeJacobian> shapeJac;
        localBasis.evaluateJacobian(ipLocal, shapeJac);

        const auto jacInvT = geom.jacobianInverseTransposed(ipLocal);

        Tensor F(0.0);
        for (const auto& localDof : localDofs(fvGeomMom))
        {
            GlobalPosition gradN;
            jacInvT.mv(shapeJac[localDof.index()][0], gradN);
            for (int dir = 0; dir < dim; ++dir)
                F[dir].axpy(
                    this->curSol(momentumIdx)[localDof.dofIndex()][dir],
                    gradN);
        }
        for (int dir = 0; dir < dim; ++dir)
            F[dir][dir] += 1.0;
        return F;
    }

    void buildCouplingStencils_(const GridGeometry<momentumIdx_>& momGG, const GridGeometry<massIdx_>& massGG)
    {
        const std::size_t numElements = massGG.gridView().size(0);
        stencilsMomToPres_.clear();
        stencilsPresToMom_.clear();
        stencilsMomToPres_.resize(numElements);
        stencilsPresToMom_.resize(numElements);

        auto fvGeomMom = localView(momGG);
        auto fvGeomPres = localView(massGG);

        for (const auto& element : elements(massGG.gridView()))
        {
            fvGeomMom.bindElement(element);
            fvGeomPres.bindElement(element);

            const std::size_t eIdx = massGG.elementMapper().index(element);

            // Momentum → Pressure stencil: pressure vertex DOFs of this element.
            auto& momToPres = stencilsMomToPres_[eIdx];
            for (const auto& scv : scvs(fvGeomPres))
                momToPres.push_back(scv.dofIndex());

            // Pressure → Momentum stencil: all momentum DOFs of this element (vertices + bubbles).
            auto& presToMom = stencilsPresToMom_[eIdx];
            for (const auto& scv : scvs(fvGeomMom))
                presToMom.push_back(scv.dofIndex());
        }
    }

    //! coloring for multithreaded assembly
    std::deque<std::vector<ElementSeed<momentumIdx>>> elementSets_;

    std::vector<CouplingStencilType> stencilsMomToPres_; //!< indexed by element
    std::vector<CouplingStencilType> stencilsPresToMom_; //!< indexed by element

    std::vector<Scalar> prevJ_; //!< J^n stored at each pressure (Box) DOF
};

// we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<GripperCouplingManager<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif // DUMUX_GRIPPER_MODEL_HH
