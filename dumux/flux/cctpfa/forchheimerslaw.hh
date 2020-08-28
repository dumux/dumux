// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CCTpfaFlux
 * \brief Forchheimers's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FORCHHEIMERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FORCHHEIMERS_LAW_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/cctpfa/darcyslaw.hh>

namespace Dumux {

// forward declarations
template<class TypeTag, DiscretizationMethod discMethod>
class ForchheimersLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Forchheimer's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Forchheimer's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 * \tparam isNetwork whether we are computing on a network grid embedded in a higher world dimension
 */
template<class Scalar, class GridGeometry, bool isNetwork>
class CCTpfaForchheimersLaw;

/*!
 * \ingroup CCTpfaFlux
 * \brief Forchheimer's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Forchheimer's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 */
template <class TypeTag>
class ForchheimersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>
: public CCTpfaForchheimersLaw<GetPropType<TypeTag, Properties::Scalar>,
                         GetPropType<TypeTag, Properties::GridGeometry>,
                         (GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension < GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld)>
{};

/*!
 * \ingroup CCTpfaFlux
 * \brief Class that fills the cache corresponding to tpfa Forchheimer's Law
 */
template<class GridGeometry>
class TpfaForchheimersLawCacheFiller
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! Function to fill a TpfaForchheimersLawCache of a given scvf
    //! This interface has to be met by any advection-related cache filler class
    //! TODO: Probably get cache type out of the filler
    template<class FluxVariablesCache, class Problem, class ElementVolumeVariables, class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {
        scvfFluxVarsCache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf);
    }
};

/*!
 * \ingroup CCTpfaFlux
 * \brief The cache corresponding to tpfa Forchheimer's Law
 */
template<class AdvectionType, class GridGeometry>
class TpfaForchheimersLawCache
{
    using Scalar = typename AdvectionType::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using Filler = TpfaForchheimersLawCacheFiller<GridGeometry>;

    template<class Problem, class ElementVolumeVariables>
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
        harmonicMeanSqrtK_ = AdvectionType::calculateHarmonicMeanSqrtPermeability(problem, elemVolVars, scvf);
    }

    const Scalar& advectionTij() const
    { return tij_; }

    const DimWorldMatrix& harmonicMeanSqrtPermeability() const
    { return harmonicMeanSqrtK_; }

private:
    Scalar tij_;
    DimWorldMatrix harmonicMeanSqrtK_;
};

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaForchheimersLaw grids where dim=dimWorld
 */
template<class ScalarType, class GridGeometry>
class CCTpfaForchheimersLaw<ScalarType, GridGeometry, /*isNetwork*/ false>
{
    using ThisType = CCTpfaForchheimersLaw<ScalarType, GridGeometry, /*isNetwork*/ false>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimWorldVector = Dune::FieldVector<ScalarType, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<ScalarType, dimWorld, dimWorld>;

    using DarcysLaw = CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ false>;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaForchheimersLawCache<ThisType, GridGeometry>;

    /*! \brief Compute the advective flux of a phase across
    *          the given sub-control volume face uing the Forchheimer equation.
    *
    *          The flux is given in N*m, and can be converted
    *          into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
    *          for the mobility or the product of density and mobility, respectively.
    *
    * see e.g. Nield & Bejan: Convection in Porous Media \cite nield2006
    *
    * The relative passability \f$ \eta_r\f$ is the "Forchheimer-equivalent" to the relative
    * permeability \f$ k_r\f$.
    * We use the same function as for \f$ k_r \f$ (VG, BC, linear) other authors use a simple
    * power law e.g.: \f$\eta_{rw} = S_w^3\f$
    *
    * Some rearrangements have been made to the original Forchheimer relation:
    *
    * \f[ \mathbf{v_\alpha} + c_F \sqrt{\mathbf{K}} \frac{\rho_\alpha}{\mu_\alpha }
    *     |\mathbf{v_\alpha}| \mathbf{v_\alpha}
    *     + \frac{k_{r \alpha}}{\mu_\alpha} \mathbf{K} \nabla \left(p_\alpha
    *     + \rho_\alpha g z \right)=  0
    * \f]
    *
    * This already includes the assumption \f$ k_r(S_w) = \eta_r(S_w)\f$:
    * - \f$\eta_{rw} = S_w^x\f$ looks very similar to e.g. Van Genuchten relative permeabilities
    * - Fichot, et al. (2006), Nuclear Engineering and Design, state that several authors claim
    *   that \f$ k_r(S_w), \eta_r(S_w)\f$ can be chosen equal
    * - It leads to the equation not degenerating for the case of \f$S_w=1\f$, because I do not
    *   need to multiply with two different functions, and therefore there are terms not being
    *   zero.
    * - If this assumption is not to be made: Some regularization needs to be introduced ensuring
    *   that not all terms become zero for \f$S_w=1\f$.
    *
    * This non-linear equations is solved for \f$\mathbf{v_\alpha}\f$ using Newton's method
    * and an analytical derivative w.r.t. \f$\mathbf{v_\alpha}\f$.
    *
    * The gradient of the Forchheimer relations looks as follows (mind that \f$\sqrt{\mathbf{K}}\f$
    * is a tensor):
    *
    * \f[  f\left(\mathbf{v_\alpha}\right) =
    * \left(
    * \begin{array}{ccc}
    * 1 & 0 &0 \\
    * 0 & 1 &0 \\
    * 0 & 0 &1 \\
    * \end{array}
    * \right)
    * +
    * c_F \frac{\rho_\alpha}{\mu_\alpha} |\mathbf{v}_\alpha| \sqrt{\mathbf{K}}
    * +
    * c_F \frac{\rho_\alpha}{\mu_\alpha}\frac{1}{|\mathbf{v}_\alpha|} \sqrt{\mathbf{K}}
    * \left(
    * \begin{array}{ccc}
    * v_x^2 & v_xv_y & v_xv_z \\
    * v_yv_x & v_{y}^2 & v_yv_z \\
    * v_zv_x & v_zv_y &v_{z}^2 \\
    * \end{array}
    * \right)
    *  \f]
    *
    * \note We restrict the use of Forchheimer's law to diagonal permeability tensors so far. This might be changed to
    * general tensors using eigenvalue decomposition to get \f$\sqrt{\mathbf{K}}\f$
    */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto velocity = velocity_(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);
        Scalar flux = velocity * scvf.unitOuterNormal();
        flux *= Extrusion::area(scvf);
        return flux;
    }

    //! The flux variables cache has to be bound to an element prior to flux calculations
    //! During the binding, the transmissibility will be computed and stored using the method below.
    template<class Problem, class ElementVolumeVariables>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        return DarcysLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    /*! \brief Returns the harmonic mean of \f$\sqrt{K_0}\f$ and \f$\sqrt{K_1}\f$.
     *
     * This is a specialization for scalar-valued permeabilities which returns a tensor with identical diagonal entries.
     */
    template<class Problem, class ElementVolumeVariables,
             bool scalarPerm = std::is_same<typename Problem::SpatialParams::PermeabilityType, Scalar>::value,
             std::enable_if_t<scalarPerm, int> = 0>
    static DimWorldMatrix calculateHarmonicMeanSqrtPermeability(const Problem& problem,
                                                                const ElementVolumeVariables& elemVolVars,
                                                                const SubControlVolumeFace& scvf)
    {
        using std::sqrt;
        Scalar harmonicMeanSqrtK(0.0);

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const Scalar Ki = getPermeability_(problem, insideVolVars, scvf.ipGlobal());
        const Scalar sqrtKi = sqrt(Ki);

        if (!scvf.boundary())
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar Kj = getPermeability_(problem, outsideVolVars, scvf.ipGlobal());
            const Scalar sqrtKj = sqrt(Kj);
            harmonicMeanSqrtK = problem.spatialParams().harmonicMean(sqrtKi, sqrtKj, scvf.unitOuterNormal());
        }
        else
            harmonicMeanSqrtK = sqrtKi;

        // Convert to tensor
        DimWorldMatrix result(0.0);
        for (int i = 0; i < dimWorld; ++i)
            result[i][i] = harmonicMeanSqrtK;

        return result;
    }

    /*! \brief Returns the harmonic mean of \f$\sqrt{\mathbf{K_0}}\f$ and \f$\sqrt{\mathbf{K_1}}\f$.
     *
     * This is a specialization for tensor-valued permeabilities.
     */
    template<class Problem, class ElementVolumeVariables,
             bool scalarPerm = std::is_same<typename Problem::SpatialParams::PermeabilityType, Scalar>::value,
             std::enable_if_t<!scalarPerm, int> = 0>
    static DimWorldMatrix calculateHarmonicMeanSqrtPermeability(const Problem& problem,
                                                                const ElementVolumeVariables& elemVolVars,
                                                                const SubControlVolumeFace& scvf)
    {
        using std::sqrt;
        DimWorldMatrix sqrtKi(0.0);
        DimWorldMatrix sqrtKj(0.0);
        DimWorldMatrix harmonicMeanSqrtK(0.0);

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto& Ki = getPermeability_(problem, insideVolVars, scvf.ipGlobal());
        // Make sure the permeability matrix does not have off-diagonal entries
        // TODO implement method using eigenvalues and eigenvectors for general tensors
        if (!isDiagonal_(Ki))
            DUNE_THROW(Dune::NotImplemented, "Only diagonal permeability tensors are supported.");

        for (int i = 0; i < dim; ++i)
            sqrtKi[i][i] = sqrt(Ki[i][i]);

        if (!scvf.boundary())
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto& Kj = getPermeability_(problem, outsideVolVars, scvf.ipGlobal());
            // Make sure the permeability matrix does not have off-diagonal entries
            if (!isDiagonal_(Kj))
                DUNE_THROW(Dune::NotImplemented, "Only diagonal permeability tensors are supported.");

            for (int i = 0; i < dim; ++i)
                sqrtKj[i][i] = sqrt(Kj[i][i]);

            harmonicMeanSqrtK = problem.spatialParams().harmonicMean(sqrtKi, sqrtKj, scvf.unitOuterNormal());
        }
        else
            harmonicMeanSqrtK = sqrtKi;

        return harmonicMeanSqrtK;
    }

private:

    //! Computes the face velocity based on the Forchheimer equation
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static DimWorldVector velocity_(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    int phaseIdx,
                                    const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // Get the volume flux based on Darcy's law. The value returned by this method needs to be multiplied with the
        // mobility (upwinding).
        Scalar darcyFlux = DarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);
        auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.mobility(phaseIdx); };
        const Scalar darcyUpwindMobility = doUpwind_(scvf, elemVolVars, upwindTerm, !std::signbit(darcyFlux) /*insideIsUpstream*/);
        darcyFlux *= darcyUpwindMobility;

        // Convert to Darcy velocity
        DimWorldVector darcyVelocity = scvf.unitOuterNormal();
        darcyVelocity *= darcyFlux / Extrusion::area(scvf);

        // Get the harmonic mean of the square root of permeability
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& sqrtK = fluxVarsCache.harmonicMeanSqrtPermeability();

        // Obtain the Forchheimer coefficient from the spatial parameters
        const Scalar forchCoeff = problem.spatialParams().forchCoeff(scvf);

        // Initial guess of velocity: Darcy relation
        DimWorldVector velocity = darcyVelocity;

        DimWorldVector deltaV(0.0);           // the change in velocity between Newton iterations
        DimWorldVector residual(10e10);  // the residual (function value that is to be minimized)
        DimWorldMatrix  gradF(0.0);            // slope of equation that is to be solved

        // Search by means of the Newton method for a root of Forchheimer equation
        static const Scalar epsilon = getParamFromGroup<Scalar>(problem.paramGroup(), "Forchheimer.NewtonTolerance", 1e-12);
        static const std::size_t maxNumIter = getParamFromGroup<std::size_t>(problem.paramGroup(), "Forchheimer.MaxIterations", 30);
        for (int k = 0; residual.two_norm() > epsilon ; ++k)
        {
            if (k >= maxNumIter)
                DUNE_THROW(NumericalProblem, "could not determine forchheimer velocity within "
                                             << k << " iterations");

            // calculate current value of Forchheimer equation
            forchheimerResidual_(residual,
                                 forchCoeff,
                                 sqrtK,
                                 darcyVelocity,
                                 velocity,
                                 elemVolVars,
                                 scvf,
                                 phaseIdx);

            // newton's method: slope (gradF) and current value (residual) of function is known,
            forchheimerDerivative_(gradF,
                                   forchCoeff,
                                   sqrtK,
                                   velocity,
                                   elemVolVars,
                                   scvf,
                                   phaseIdx);

            // solve for change in velocity ("x-Axis intercept")
            gradF.solve(deltaV, residual);
            velocity -= deltaV;
        }

        // The fluxvariables expect a value on which an upwinding of the mobility is performed.
        // We therefore divide by the mobility here.
        const Scalar forchheimerUpwindMobility = doUpwind_(scvf, elemVolVars, upwindTerm,
                                                           !std::signbit(velocity * scvf.unitOuterNormal()) /*insideIsUpstream*/);

        // Do not divide by zero. If the mobility is zero, the resulting flux with upwinding will be zero anyway.
        if (forchheimerUpwindMobility > 0.0)
            velocity /= forchheimerUpwindMobility;

        return velocity;
    }

     //! calculate the residual of the Forchheimer equation
     template<class ElementVolumeVariables>
     static void forchheimerResidual_(DimWorldVector& residual,
                                      const Scalar forchCoeff,
                                      const DimWorldMatrix& sqrtK,
                                      const DimWorldVector& darcyVelocity,
                                      const DimWorldVector& oldVelocity,
                                      const ElementVolumeVariables& elemVolVars,
                                      const SubControlVolumeFace& scvf,
                                      const int phaseIdx)
     {
         residual = oldVelocity;

         // residual += k_r/mu  K grad p
         residual -= darcyVelocity;

         // residual += c_F rho / mu abs(v) sqrt(K) v
         auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.density(phaseIdx) / volVars.viscosity(phaseIdx); };
         bool insideIsUpstream = !std::signbit(oldVelocity * scvf.unitOuterNormal());
         const Scalar rhoOverMu = doUpwind_(scvf, elemVolVars, upwindTerm, insideIsUpstream);
         sqrtK.usmv(forchCoeff * rhoOverMu * oldVelocity.two_norm(), oldVelocity, residual);
     }

     //! The analytical derivative of the the Forcheimer equation's residual
    template<class ElementVolumeVariables>
    static void forchheimerDerivative_(DimWorldMatrix& derivative,
                                       const Scalar forchCoeff,
                                       const DimWorldMatrix& sqrtK,
                                       const DimWorldVector& velocity,
                                       const ElementVolumeVariables& elemVolVars,
                                       const SubControlVolumeFace& scvf,
                                       const int phaseIdx)
    {


        // Initialize because for low velocities, we add and do not set the entries.
        derivative = 0.0;

        auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.density(phaseIdx) / volVars.viscosity(phaseIdx); };
        bool insideIsUpstream = !std::signbit(velocity * scvf.unitOuterNormal());

        const Scalar absV = velocity.two_norm() ;
        // This part of the derivative is only calculated if the velocity is sufficiently small:
        // do not divide by zero!
        // This should not be very bad because the derivative is only intended to give an
        // approximation of the gradient of the
        // function that goes into the Newton scheme.
        // This is important in the case of a one-phase region in two-phase flow. The non-existing
        // phase has a velocity of zero (kr=0).
        // f = sqrtK * vTimesV (this is  a matrix)
        // f *= forchCoeff density / |velocity| / viscosity (multiply all entries with scalars)
        const Scalar rhoOverMu = doUpwind_(scvf, elemVolVars, upwindTerm, insideIsUpstream);
        if (absV > 1e-20)
        {
            for (int i = 0; i < dim; i++)
            {
                for (int k = 0; k <= i; k++)
                {
                    derivative[i][k] = sqrtK[i][i] * velocity[i] * velocity[k] * forchCoeff
                    / absV  * rhoOverMu;

                    if (k < i)
                        derivative[k][i] = derivative[i][k];
                }
            }
        }

        // add on the main diagonal:
        // 1 + sqrtK_i forchCoeff density |v| / viscosity
        for (int i = 0; i < dim; i++)
            derivative[i][i] += (1.0 + (sqrtK[i][i]*forchCoeff * absV * rhoOverMu));
    }

     template<class ElementVolumeVariables, class UpwindTermFunction>
     static Scalar doUpwind_(const SubControlVolumeFace& scvf,
                             const ElementVolumeVariables& elemVolVars,
                             const UpwindTermFunction& upwindTerm,
                             const bool insideIsUpstream)
     {
         static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

         const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
         const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

         if (!insideIsUpstream) // if sign of flux is negative
             return (upwindWeight*upwindTerm(outsideVolVars)
                          + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
         else
             return (upwindWeight*upwindTerm(insideVolVars)
                          + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
     }

    /*!
     * \brief Check whether all off-diagonal entries of a tensor are zero.
     *
     * \param K the tensor that is to be checked.
     *
     * \return True if all off-diagonals are zero.
     *
    */
    static bool isDiagonal_(const DimWorldMatrix & K)
    {
        using std::abs;
        for (int i = 0; i < dim; i++) {
            for (int k = 0; k < dim; k++) {
                if ((i != k) && (abs(K[i][k]) >= 1e-25)) {
                  return false;
                }
            }
        }
        return true;
    }

    template<class Problem, class VolumeVariables,
             std::enable_if_t<!Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return volVars.permeability(); }

    template<class Problem, class VolumeVariables,
             std::enable_if_t<Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return problem.spatialParams().permeabilityAtPos(scvfIpGlobal); }
};

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaForchheimersLaw grids where dim<dimWorld
 */
template<class ScalarType, class GridGeometry>
class CCTpfaForchheimersLaw<ScalarType, GridGeometry, /*isNetwork*/ true>
{
    static_assert(AlwaysFalse<ScalarType>::value, "Forchheimer not implemented for network grids");
};

} // end namespace Dumux

#endif
