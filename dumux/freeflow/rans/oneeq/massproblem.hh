// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::RANSMassOneEqProblem
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_PROBLEM_HH
#define DUMUX_RANS_ONEEQ_MASS_PROBLEM_HH

#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the one-equation (Spalart-Allmaras)
 *        RANS turbulence closure: the Spalart-Allmaras production/destruction/
 *        cross-diffusion source term (through the standard problem.source() extension
 *        point), and the once-per-time-step (lagged, see below) update of the ν̃ gradient
 *        needed for the cross-diffusion term.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::RANSMomentumProblem
 *         or, in practice, Dumux::OneEqMomentumProblem, subclass) - explicit template
 *         parameter rather than property-system plumbing, since the mass sub-model has no
 *         other way to learn its momentum sibling's type. Set post-construction via
 *         setMomentumProblem(), mirroring how e.g. Problem::setTimeLoop() is wired in
 *         main.cc elsewhere in this codebase.
 *
 * wallDistance()/vorticityTensorScalarProduct() are *not* recomputed here - they are frozen
 * (lagged) values computed once per time step by the momentum domain's
 * Dumux::RANSMomentumProblem and simply forwarded, read-only, through the stored
 * momentum-problem pointer (both grids share the same underlying element indexing). This is
 * deliberately *not* routed through a CouplingManager: since these values are frozen for the
 * whole time step, there is nothing for Newton to differentiate against.
 *
 * The gradient ∇ν̃ is reconstructed like releases/3.10 did: a central difference over the two
 * neighboring elements along each axis (structured, axis-aligned meshes), and, at Dirichlet
 * boundaries, a difference between the boundary value and the next interior element (see
 * updateDynamicWallProperties() below).
 */
template<class TypeTag, class MomentumProblem>
class RANSMassOneEqProblem : public NavierStokesMassProblem<TypeTag>
{
    using ParentType = NavierStokesMassProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    static constexpr int dim = GridView::dimension;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    template<class... Args>
    RANSMassOneEqProblem(Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    // ν̃'s gradient is read (OneEqMassVolumeVariables::update()) as soon as the grid variables
    // are first initialized in main.cc, i.e. before updateDynamicWallProperties() is ever
    // called (that only happens once per time step, inside the time loop) - so this must be
    // pre-sized here, to zero, rather than left empty until the first update.
    , viscosityTildeGradient_(this->gridGeometry().elementMapper().size(), DimVector(0.0))
    {}

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar wallDistance(std::size_t eIdx) const
    { return momentumProblem_->wallDistance(eIdx); }

    Scalar vorticityTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->vorticityTensorScalarProduct(eIdx); }

    DimVector storedViscosityTildeGradient(std::size_t eIdx) const
    { return viscosityTildeGradient_[eIdx]; }

    /*!
     * \brief Updates the (lagged) ν̃ gradient used by the cross-diffusion source term, from
     *        the given (typically: the last converged) mass solution. Called once per time
     *        step from main.cc, before the next assembly - not re-evaluated within Newton,
     *        matching how releases/3.10 treated this term (see class documentation above).
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& sol)
    {
        const auto& gridGeometry = this->gridGeometry();
        const auto numElements = gridGeometry.elementMapper().size();
        viscosityTildeGradient_.assign(numElements, DimVector(0.0));

        std::vector<Scalar> viscosityTilde(numElements);
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            const auto& insideScv = *scvs(fvGeometry).begin();
            viscosityTilde[insideScv.elementIndex()] = sol[insideScv.dofIndex()][Indices::viscosityTildeIdx];
        }

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bind(element);

            // central difference over the two neighbors along each axis (a neighbor that does not
            // exist is the element itself, which gives a one-sided difference at boundaries)
            for (int axisIdx = 0; axisIdx < dim; ++axisIdx)
            {
                const auto neighborIdx0 = momentumProblem_->neighborIndex(eIdx, axisIdx, 0);
                const auto neighborIdx1 = momentumProblem_->neighborIndex(eIdx, axisIdx, 1);
                const auto distance = momentumProblem_->cellCenter(neighborIdx1)[axisIdx]
                                     - momentumProblem_->cellCenter(neighborIdx0)[axisIdx];
                if (std::abs(distance) < 1e-8)
                    continue;

                viscosityTildeGradient_[eIdx][axisIdx] = (viscosityTilde[neighborIdx1] - viscosityTilde[neighborIdx0])/distance;
            }

            // Dirichlet boundaries (inlet and, as a weak condition, the walls where nu-tilde is zero):
            // difference the boundary value against the next interior element instead
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary())
                    continue;

                Scalar boundaryValue;
                if (this->asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::viscosityTildeEqIdx))
                    boundaryValue = this->asImp_().dirichlet(element, scvf)[Indices::viscosityTildeIdx];
                else if (this->asImp_().isOnWallAtPos(scvf.center()))
                    boundaryValue = 0.0;
                else
                    continue;

                const auto& normal = scvf.unitOuterNormal();
                const int axisIdx = std::abs(normal[0]) > 0.5 ? 0 : 1;
                auto neighborIdx = momentumProblem_->neighborIndex(eIdx, axisIdx, 0);
                if (scvf.center()[axisIdx] < momentumProblem_->cellCenter(eIdx)[axisIdx])
                    neighborIdx = momentumProblem_->neighborIndex(eIdx, axisIdx, 1);

                viscosityTildeGradient_[eIdx][axisIdx] = (viscosityTilde[neighborIdx] - boundaryValue)
                    /(momentumProblem_->cellCenter(neighborIdx)[axisIdx] - scvf.center()[axisIdx]);
            }
        }
    }

    /*!
     * \brief The Spalart-Allmaras production/destruction/cross-diffusion source term,
     *        ported from releases/3.10:dumux/freeflow/rans/oneeq/staggered/localresidual.hh's
     *        computeSourceForCellCenter(). Added through the standard problem.source(...)
     *        extension point (Dumux::DiscretizationDefaultLocalOperator's default
     *        computeSource() already forwards here generically, so no localResidual
     *        override is needed for this term, see localresidual.hh).
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];

        source[Indices::viscosityTildeEqIdx] += volVars.cb1()*volVars.stressTensorScalarProductTilde()
                                                * volVars.viscosityTilde()*volVars.density();

        source[Indices::viscosityTildeEqIdx] -= volVars.cw1()*volVars.fw()
                                                * volVars.viscosityTilde()*volVars.viscosityTilde()
                                                / (volVars.wallDistance()*volVars.wallDistance())*volVars.density();

        source[Indices::viscosityTildeEqIdx] += volVars.cb2()/volVars.sigma()
                                                * volVars.viscosityTildeGradientSquared()*volVars.density();

        return source;
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::shared_ptr<const MomentumProblem> momentumProblem_;
    std::vector<DimVector> viscosityTildeGradient_;
};

} // end namespace Dumux

#endif
