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
 * The gradient ∇ν̃ is obtained via a standard Green-Gauss reconstruction over CCTpfa's own
 * face/neighbor connectivity - works on general grids, no flat-wall-bounded assumption needed
 * (see updateDynamicWallProperties() below).
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
        viscosityTildeGradient_.assign(gridGeometry.elementMapper().size(), DimVector(0.0));

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            const auto& insideScv = *scvs(fvGeometry).begin();
            const auto eIdx = insideScv.elementIndex();
            const Scalar viscosityTildeInside = sol[insideScv.dofIndex()][Indices::viscosityTildeIdx];

            DimVector gradient(0.0);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                Scalar viscosityTildeFace;
                if (scvf.boundary())
                {
                    if (this->asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::viscosityTildeEqIdx))
                        viscosityTildeFace = this->asImp_().dirichlet(element, scvf)[Indices::viscosityTildeIdx];
                    else
                        viscosityTildeFace = viscosityTildeInside; // zero-gradient extrapolation
                }
                else
                {
                    const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                    const Scalar viscosityTildeOutside = sol[outsideScv.dofIndex()][Indices::viscosityTildeIdx];
                    viscosityTildeFace = 0.5*(viscosityTildeInside + viscosityTildeOutside);
                }

                gradient.axpy((viscosityTildeFace - viscosityTildeInside)*scvf.area(), scvf.unitOuterNormal());
            }
            gradient /= insideScv.volume();
            viscosityTildeGradient_[eIdx] = gradient;
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
