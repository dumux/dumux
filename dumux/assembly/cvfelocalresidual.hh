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
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Calculates the element-wise residual for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_CVFE_LOCAL_RESIDUAL_HH

#include <dune/common/std/type_traits.hh>
#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/fvlocalresidual.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux::Detail {

template<class Imp, class P, class G, class S, class V>
using TimeInfoInterfaceCVFEDetector = decltype(
    std::declval<Imp>().computeStorage(
        std::declval<P>(), std::declval<G>(), std::declval<S>(), std::declval<V>(), true
    )
);

template<class Imp, class P, class G, class S, class V>
constexpr inline bool hasTimeInfoInterfaceCVFE()
{ return Dune::Std::is_detected<TimeInfoInterfaceCVFEDetector, Imp, P, G, S, V>::value; }

template<class Imp>
using SCVFIsOverlappingDetector = decltype(
    std::declval<Imp>().isOverlapping()
);

template<class Imp>
constexpr inline bool hasScvfIsOverlapping()
{ return Dune::Std::is_detected<SCVFIsOverlappingDetector, Imp>::value; }

} // end namespace Dumux::Detail


namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief The element-wise residual for control-volume finite element schemes
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class CVFELocalResidual : public FVLocalResidual<TypeTag>
{
    using ParentType = FVLocalResidual<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;
public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate flux residuals for one sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const
    {
        if constexpr (Detail::hasScvfIsOverlapping<SubControlVolumeFace>())
        {
            if (!scvf.isOverlapping())
            {
                const auto flux = evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
                if (!scvf.boundary())
                {
                    const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                    residual[insideScv.localDofIndex()] += flux;
                    residual[outsideScv.localDofIndex()] -= flux;
                }
                else
                {
                    const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                    residual[insideScv.localDofIndex()] += flux;
                }
            }
        }
        else
        {
            const auto flux = evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
            if (!scvf.boundary())
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                residual[insideScv.localDofIndex()] += flux;
                residual[outsideScv.localDofIndex()] -= flux;
            }
            else
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                residual[insideScv.localDofIndex()] += flux;
            }
        }
    }

    //! evaluate flux residuals for one sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementBoundaryTypes& elemBcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

        // inner faces
        if (!scvf.boundary())
            flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        // boundary faces
        else
        {
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& bcTypes = elemBcTypes.get(fvGeometry, scv);

            // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
            // For Dirichlet there is no addition to the residual here but they
            // are enforced strongly by replacing the residual entry afterwards.
            if (bcTypes.hasNeumann())
            {
                auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                neumannFluxes *= Extrusion::area(fvGeometry, scvf)*elemVolVars[scv].extrusionFactor();

                // only add fluxes to equations for which Neumann is set
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        flux[eqIdx] += neumannFluxes[eqIdx];
            }
        }

        return flux;
    }

    using ParentType::evalStorage;
    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // Compute storage with the model specific storage residual
        // This addresses issues #792/#940 in ad-hoc way by additionally providing crude time level information (previous or current)
        // to the low-level interfaces if this is supported by the LocalResidual implementation
        NumEqVector prevStorage = computeStorageImpl_(problem, fvGeometry, scv, prevVolVars, /*previous time level?*/true);
        NumEqVector storage = computeStorageImpl_(problem, fvGeometry, scv, curVolVars, /*previous time level?*/false);

        prevStorage *= prevVolVars.extrusionFactor();
        storage *= curVolVars.extrusionFactor();

        storage -= prevStorage;
        storage *= Extrusion::volume(fvGeometry, scv);
        storage /= this->timeLoop().timeStepSize();

        residual[scv.localDofIndex()] += storage;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        auto source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);

        if constexpr (Detail::hasScvfIsOverlapping<SubControlVolume>())
        {
            if (scv.isOverlapping())
            {
                const auto geometry = element.geometry();
                const auto& localBasis = fvGeometry.feLocalBasis();

                const auto h = diameter(geometry);
                auto velocity = elemVolVars[scv].velocity();
                static const auto stabParam = getParam<Scalar>("Problem.StabilizationParameter", 1.0);
                velocity *= 1.0/(stabParam*h*h);
                source -= velocity;

                // compute the gradN in the center of the element
                const auto ipLocal = geometry.local(scv.center());
                const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
                std::vector<ShapeJacobian> shapeJacobian;
                localBasis.evaluateJacobian(ipLocal, shapeJacobian);
                std::vector<GlobalPosition> gradN;
                gradN.resize(fvGeometry.numScv());
                for (const auto& scv : scvs(fvGeometry))
                    jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN[scv.indexInElement()]);

                Dune::FieldVector<Scalar, GridView::dimensionworld> gradP(0.0);
                for (const auto& scv : scvs(fvGeometry))
                    if (!scv.isOverlapping())
                        gradP.axpy(problem.couplingManager().pressure(element, fvGeometry, scv), gradN[scv.indexInElement()]);

                source -= gradP;
            }
        }

        return source;
    }

private:
    NumEqVector computeStorageImpl_(const Problem& problem,
                                    const FVElementGeometry& fvGeometry,
                                    const SubControlVolume& scv,
                                    const VolumeVariables& volVars,
                                    [[maybe_unused]] bool isPreviousTimeStep) const
    {
        if constexpr (Detail::hasTimeInfoInterfaceCVFE<Implementation, Problem, FVElementGeometry, SubControlVolume, VolumeVariables>())
            return this->asImp().computeStorage(problem, fvGeometry, scv, volVars, isPreviousTimeStep);
        else
            return this->asImp().computeStorage(problem, scv, volVars);
    }
};

} // end namespace Dumux

#endif
