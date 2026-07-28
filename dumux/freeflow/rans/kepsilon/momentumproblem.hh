// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KEpsilonMomentumProblem
 */
#ifndef DUMUX_RANS_KEPSILON_MOMENTUM_PROBLEM_HH
#define DUMUX_RANS_KEPSILON_MOMENTUM_PROBLEM_HH

#include <cmath>
#include <functional>

#include <dumux/common/properties.hh>
#include <dumux/freeflow/rans/common/ransmomentumproblem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the momentum-domain side of the two-equation, high-Reynolds-number
 *        wall-function k-epsilon RANS turbulence closure to a face-centered-staggered
 *        Navier-Stokes momentum problem: the wall-distance/velocity-gradient/stress-tensor
 *        bookkeeping inherited from Dumux::RANSMomentumProblem, forwarding the turbulent
 *        viscosity and turbulent kinetic energy (both computed on the *mass* sub-model, see
 *        dumux/freeflow/rans/kepsilon/volumevariables.hh) into effectiveViscosity() and the
 *        normal-stress correction respectively (unchanged from Dumux::KOmegaMomentumProblem),
 *        plus the log-law momentum wall function applied *only* at the matching-point cell of
 *        each wall-normal column (see dumux/freeflow/rans/kepsilon/massproblem.hh for what
 *        "matching point" means).
 *
 * The wall function is implemented as a plain Neumann momentum flux at the matching-point wall
 * face (the same mechanism already used for the outlet's fixedPressureMomentumFlux, just gated
 * to matching-point cells instead) rather than a dumux/freeflow/navierstokes/momentum/
 * slipcondition.hh SlipConditions policy - see whatisimplemented.md/proposedimplementation.md
 * for why: the slip-velocity mechanism there solves for an *equivalent slip velocity* that
 * reproduces a desired flux through the standard finite-difference viscous-flux formula, a
 * nontrivial reformulation of the deleted releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/
 * problem.hh's wallFunction(), which instead returns the shear-stress flux directly.
 *
 * The matching-point classification and uStarNominal/yPlusNominal values are lagged (computed
 * once per time step on the mass domain, see Dumux::KEpsilonMassProblem) - reading them via the
 * coupling manager (as effectiveViscosity()/turbulentKineticEnergy() above legitimately do,
 * since those are meant to be live/differentiated) triggered a segfault here: boundaryTypes()
 * is evaluated very early, as part of building the momentum grid's own flux-variables cache,
 * from a multithreaded pass that isn't compatible with the coupling manager's single-context
 * binding cache. Since this data is lagged anyway (no Newton differentiation needed - the same
 * reasoning already used for wallDistance()/stressTensorScalarProduct() etc.), it is instead
 * forwarded directly via plain, type-erased callbacks (setWallFunctionQueries(), wired once in
 * main.cc after both sub-problems exist) - this also sidesteps what would otherwise be a
 * circular template dependency (a concrete Dumux::KEpsilonMassProblem<TypeTag, MomentumProblem>
 * template parameter here would need MomentumProblem's own concrete type, which in turn would
 * need the mass problem's type, and so on).
 */
template<class TypeTag>
class KEpsilonMomentumProblem : public RANSMomentumProblem<TypeTag>
{
    using ParentType = RANSMomentumProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using BoundaryTypes = typename ParentType::BoundaryTypes;

    using ParentType::ParentType;

    //! Type-erased read-only queries into the mass domain's (lagged) k-epsilon wall-function
    //! bookkeeping - see class documentation above for why this isn't routed through the
    //! coupling manager. Must be called once, after both sub-problems have been constructed.
    void setWallFunctionQueries(std::function<bool(std::size_t)> isMatchingPoint,
                               std::function<Scalar(std::size_t)> uStarNominal,
                               std::function<Scalar(std::size_t)> yPlusNominal)
    {
        isMatchingPointQuery_ = std::move(isMatchingPoint);
        uStarNominalQuery_ = std::move(uStarNominal);
        yPlusNominalQuery_ = std::move(yPlusNominal);
    }

    //! Adds the (current-iterate) turbulent viscosity, read from the mass domain's k/epsilon
    //! through the coupling manager, to the molecular one.
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scvf) + this->couplingManager().turbulentViscosity(element, fvGeometry, scvf); }

    //! \copydoc effectiveViscosity(const Element&,const FVElementGeometry&,const SubControlVolumeFace&) const
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scv) + this->couplingManager().turbulentViscosity(element, fvGeometry, scv); }

    //! The turbulent kinetic energy k, read live from the mass domain through the coupling
    //! manager - used by the isotropic 2/dim*rho*k momentum normal-stress correction.
    Scalar turbulentKineticEnergy(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf) const
    { return this->couplingManager().turbulentKineticEnergy(element, fvGeometry, scvf); }

    /*!
     * \brief At the matching-point cell's own wall (lateral, boundary) face, switch the
     *        velocity boundary condition from the usual Dirichlet no-slip to Neumann, so the
     *        log-law wall-shear-stress flux (wallFunctionMomentumFlux() below) can be applied
     *        there instead - exactly mirroring how releases/3.10's useWallFunction() gated its
     *        own wallFunction() hook.
     */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (scvf.isLateral() && scvf.boundary() && this->isOnWall(scvf) && isMatchingPointQuery_)
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            if (isMatchingPointQuery_(eIdx))
            {
                BoundaryTypes values;
                values.setAllNeumann();
                return values;
            }
        }

        return ParentType::boundaryTypes(element, scvf);
    }

    /*!
     * \brief The log-law wall-shear-stress momentum flux, ported from the deleted
     *        releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/problem.hh's
     *        tangentialMomentumWallFunction()/wallFunction(). Only ever called at the
     *        matching-point cell's own wall face (see boundaryTypes() above) - the most derived
     *        test problem's neumann() must dispatch here for wall faces, falling back to its
     *        own (e.g. outlet) Neumann treatment elsewhere.
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    Scalar wallFunctionMomentumFlux(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFluxVariablesCache& elemFluxVarsCache,
                                    const SubControlVolumeFace& scvf) const
    {
        using std::log;
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto eIdx = scv.elementIndex();
        const auto& volVars = elemVolVars[scv];

        const auto uStarNominal = uStarNominalQuery_(eIdx);
        const auto yPlusNominal = yPlusNominalQuery_(eIdx);
        const auto velocityNominal = uStarNominal*(1.0/this->karmanConstant()*log(yPlusNominal) + 5.0);
        const auto tangentialMomentumWallFunction = uStarNominal*uStarNominal*volVars.velocity()/velocityNominal;

        return tangentialMomentumWallFunction*this->storedDensity(eIdx);
    }

private:
    std::function<bool(std::size_t)> isMatchingPointQuery_;
    std::function<Scalar(std::size_t)> uStarNominalQuery_;
    std::function<Scalar(std::size_t)> yPlusNominalQuery_;
};

} // end namespace Dumux

#endif
