// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumFluxCVFE
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FLUXVARIABLES_HH

#include <concepts>

#include <dune/common/fmatrix.hh>
#include <dune/istl/scaledidmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/ipdata_.hh>

#include <dumux/discretization/cvfe/localdof.hh>

#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Context for computing fluxes
 *
 * \tparam Problem the problem type to solve
 * \tparam ElementDiscretization the element geometry type
 * \tparam ElementVariables the element variables type
 * \tparam ElementFluxVariablesCache the element flux variables cache type
 */
template<class Problem,
         class ElementDiscretization,
         class ElementVariables,
         class ElementFluxVariablesCache>
class NavierStokesMomentumFluxContext
{
    using Element = typename ElementDiscretization::Element;
    using SubControlVolumeFace = typename ElementDiscretization::SubControlVolumeFace;
public:

    //! Initialize the flux variables storing some temporary pointers
    NavierStokesMomentumFluxContext(
        const Problem& problem,
        const ElementDiscretization& elemDisc,
        const ElementVariables& elemVars,
        const ElementFluxVariablesCache& elemFluxVarsCache,
        const SubControlVolumeFace& scvf
    )
    : problem_(problem)
    , elemDisc_(elemDisc)
    , elemVars_(elemVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    , scvf_(scvf)
    {}

    const Problem& problem() const
    { return problem_; }

    const Element& element() const
    { return elemDisc_.element(); }

    const SubControlVolumeFace& scvFace() const
    { return scvf_; }

    const ElementDiscretization& elemDisc() const
    { return elemDisc_; }

    const ElementVariables& elemVars() const
    { return elemVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

private:
    const Problem& problem_;
    const ElementDiscretization& elemDisc_;
    const ElementVariables& elemVars_;
    const ElementFluxVariablesCache& elemFluxVarsCache_;
    const SubControlVolumeFace& scvf_;
};

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using control-volume finite element schemes
 */
template<class GridGeometry, class NumEqVector>
class NavierStokesMomentumFluxCVFE
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename NumEqVector::value_type;

    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    static_assert(NumEqVector::dimension == dimWorld, "Wrong dimension of velocity vector");

public:
    /*!
     * \brief Returns the advective momentum flux
     */
    template<class Context>
    NumEqVector advectiveMomentumFlux(const Context& context) const
    {
        if (!context.problem().enableInertiaTerms())
            return NumEqVector(0.0);

        const auto& elemDisc = context.elemDisc();
        const auto& elemVars = context.elemVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();

        // interpolate velocity at scvf
        NumEqVector v(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            v.axpy(shapeValues[localDof.index()][0], elemVars[localDof.index()].velocity());

        // get density from the problem
        const Scalar density = context.problem().density(context.element(), context.elemDisc(), fluxVarCache.ipData());

        const auto vn = v*scvf.unitOuterNormal();
        const auto& insideVars = elemVars[elemDisc.scv(scvf.insideScvIdx())];
        const auto& outsideVars = elemVars[elemDisc.scv(scvf.outsideScvIdx())];
        const auto upwindVelocity = vn > 0 ? insideVars.velocity() : outsideVars.velocity();
        const auto downwindVelocity = vn > 0 ? outsideVars.velocity() : insideVars.velocity();
        static const auto upwindWeight = getParamFromGroup<Scalar>(context.problem().paramGroup(), "Flux.UpwindWeight");
        const auto advectiveTermIntegrand = density*vn * (upwindWeight * upwindVelocity + (1.0-upwindWeight)*downwindVelocity);

        return advectiveTermIntegrand * Extrusion::area(elemDisc, scvf) * insideVars.extrusionFactor();
    }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    template<class Context>
    NumEqVector diffusiveMomentumFlux(const Context& context) const
    {
        const auto& element = context.element();
        const auto& elemDisc = context.elemDisc();
        const auto& elemVars = context.elemVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        // interpolate velocity gradient at scvf
        Tensor gradV(0.0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof];
            for (int dir = 0; dir < dim; ++dir)
                gradV[dir].axpy(vars.velocity(dir), fluxVarCache.gradN(localDof.index()));
        }

        // get viscosity from the problem
        const auto mu = context.problem().effectiveViscosity(element, elemDisc, fluxVarCache.ipData());

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        // compute -mu*gradV*n*dA
        NumEqVector diffusiveFlux = enableUnsymmetrizedVelocityGradient ?
                mv(gradV, scvf.unitOuterNormal())
                : mv(gradV + getTransposed(gradV),scvf.unitOuterNormal());

        diffusiveFlux *= -mu;

        static const bool enableDilatationTerm = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
            diffusiveFlux += 2.0/3.0 * mu * trace(gradV) * scvf.unitOuterNormal();

        diffusiveFlux *= Extrusion::area(elemDisc, scvf) * elemVars[elemDisc.scv(scvf.insideScvIdx())].extrusionFactor();
        return diffusiveFlux;
    }

    template<class Context>
    NumEqVector pressureContribution(const Context& context) const
    {
        const auto& element = context.element();
        const auto& elemDisc = context.elemDisc();
        const auto& elemVars = context.elemVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        // The pressure force needs to take the extruded scvf area into account
        const auto pressure = context.problem().pressure(element, elemDisc, fluxVarCache.ipData());

        // The pressure contribution calculated above might have a much larger numerical value compared to the viscous or inertial forces.
        // This may lead to numerical inaccuracies due to loss of significance (cancellation) for the final residual value.
        // In the end, we are only interested in a pressure gradient between the two relevant faces so we can
        // subtract a constant reference value from the actual pressure contribution.
        const auto referencePressure = context.problem().referencePressure();

        NumEqVector pn(scvf.unitOuterNormal());
        pn *= (pressure-referencePressure)*Extrusion::area(elemDisc, scvf)*elemVars[elemDisc.scv(scvf.insideScvIdx())].extrusionFactor();

        return pn;
    }
};

} // end namespace Dumux

namespace Dumux::Experimental {

/*!
 * \ingroup NavierStokesModel
 * \brief Momentum flux of the Navier-Stokes equations for control-volume finite element schemes,
 *        evaluated at an interpolation point
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 */
template<class Scalar, class GridDiscretization>
class CVFENavierStokesMomentumFluxAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using PressureTensor = Dune::ScaledIdentityMatrix<Scalar, dimWorld>;

    /*!
     * \brief Returns the velocity at an interpolation point
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipData The interpolation point data
     */
    template<class ElementVariables, Concept::IpData IpData>
    static GlobalPosition velocity(const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpData& ipData)
    { return velocity(elemDisc, elemVars, cache(elemVars, ipData)); }

    /*!
     * \brief Returns the velocity at an interpolation point
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipCache The interpolation point cache
     */
    template<class ElementVariables, class IpCache>
    static GlobalPosition velocity(const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpCache& ipCache)
    {
        const auto& shapeValues = ipCache.shapeValues();

        GlobalPosition velocity(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            velocity.axpy(shapeValues[localDof.index()][0], elemVars[localDof.index()].velocity());

        return velocity;
    }

    /*!
     * \brief Returns the velocity gradient at an interpolation point
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipData The interpolation point data
     */
    template<class ElementVariables, Concept::IpData IpData>
    static Tensor velocityGradient(const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpData& ipData)
    { return velocityGradient(elemDisc, elemVars, cache(elemVars, ipData)); }

    /*!
     * \brief Returns the velocity gradient at an interpolation point
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipCache The interpolation point cache
     */
    template<class ElementVariables, class IpCache>
    static Tensor velocityGradient(const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpCache& ipCache)
    {
        Tensor gradV(0.0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof.index()];
            for (int dir = 0; dir < dimWorld; ++dir)
                gradV[dir].axpy(vars.velocity(dir), ipCache.gradN(localDof.index()));
        }

        return gradV;
    }

    /*!
     * \brief Returns the diffusive part of the flux tensor at an interpolation point
     *
     * \param problem The problem
     * \param element The element
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipData The interpolation point data
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData>
    static Tensor diffusiveFluxTerm(const Problem& problem,
                                    const Element& element,
                                    const ElementDiscretization& elemDisc,
                                    const ElementVariables& elemVars,
                                    const IpData& ipData)
    { return diffusiveFluxTerm(problem, element, elemDisc, elemVars, ipData, cache(elemVars, ipData)); }

    /*!
     * \brief Returns the diffusive part of the flux tensor at an interpolation point
     *
     * \param problem The problem
     * \param element The element
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipData The interpolation point data
     * \param ipCache The interpolation point cache
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class IpCache>
    static Tensor diffusiveFluxTerm(const Problem& problem,
                                    const Element& element,
                                    const ElementDiscretization& elemDisc,
                                    const ElementVariables& elemVars,
                                    const IpData& ipData,
                                    const IpCache& ipCache)
    {
        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        static const bool enableDilatationTerm
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableDilatationTerm", false);

        const auto gradV = velocityGradient(elemDisc, elemVars, ipCache);
        const Scalar mu = problem.effectiveViscosity(element, elemDisc, ipData);

        Tensor fluxTerm = enableUnsymmetrizedVelocityGradient ? gradV : gradV + getTransposed(gradV);
        fluxTerm *= -mu;

        if (enableDilatationTerm)
        {
            const Scalar dilatation = 2.0/3.0 * mu * trace(gradV);
            for (int dir = 0; dir < dimWorld; ++dir)
                fluxTerm[dir][dir] += dilatation;
        }

        return fluxTerm;
    }

    /*!
     * \brief Returns the pressure part of the flux tensor at an interpolation point
     *
     * \param problem The problem
     * \param element The element
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     */
    template<class Problem, Concept::IpData IpData>
    static PressureTensor pressureFluxTerm(const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const IpData& ipData)
    { return PressureTensor(problem.pressure(element, elemDisc, ipData) - problem.referencePressure()); }

    /*!
     * \brief Returns the advective part of the flux tensor at an interpolation point
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param ipData The interpolation point data
     * \param transportedMomentum The transported momentum
     * \param transportingVelocity The transporting velocity
     */
    template<class ElementVariables, Concept::IpData IpData,
             class TransportedMomentum, class TransportingVelocity>
        requires std::invocable<TransportedMomentum, const ElementDiscretization&, const ElementVariables&, const IpData&>
              && std::invocable<TransportingVelocity, const ElementDiscretization&, const ElementVariables&, const IpData&>
    static Tensor advectiveFluxTerm(const ElementDiscretization& elemDisc,
                                    const ElementVariables& elemVars,
                                    const IpData& ipData,
                                    const TransportedMomentum& transportedMomentum,
                                    const TransportingVelocity& transportingVelocity)
    {
        const GlobalPosition momentum = transportedMomentum(elemDisc, elemVars, ipData);
        const GlobalPosition velocity = transportingVelocity(elemDisc, elemVars, ipData);

        Tensor fluxTerm(0.0);
        for (int i = 0; i < dimWorld; ++i)
            fluxTerm[i].axpy(momentum[i], velocity);

        return fluxTerm;
    }
};

} // end namespace Dumux::Experimental

#endif
