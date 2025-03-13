// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Local residual for the hyperelastic model
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_LOCAL_RESIDUAL_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail::Hyperelastic {
// helper struct detecting if the user-defined problem class has a solidDensity method
template <typename T, typename ...Ts>
using SolidDensityDetector = decltype(std::declval<T>().solidDensity(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool hasSolidDensity()
{ return Dune::Std::is_detected<SolidDensityDetector, T, Args...>::value; }

// helper struct detecting if the user-defined problem class has a acceleration method
template <typename T, typename ...Ts>
using AccelerationDetector = decltype(std::declval<T>().acceleration(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool hasAcceleration()
{ return Dune::Std::is_detected<AccelerationDetector, T, Args...>::value; }
} // end namespace Detail::Hyperelastic
#endif

/*!
 * \ingroup Hyperelastic
 * \brief Local residual for the hyperelastic model
 */
template<class TypeTag>
class HyperelasticLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;
public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    using ParentType::evalStorage;

    /*!
     * \brief Evaluate the storage contribution to the residual: rho*d^2u/dt^2
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

        if constexpr (Detail::Hyperelastic::hasSolidDensity<Problem, const Element&, const SubControlVolume&>()
                      && Detail::Hyperelastic::hasAcceleration<Problem, const Element&, const SubControlVolume&, Scalar, decltype(curVolVars.displacement())>())
        {
            const auto& curVolVars = curElemVolVars[scv];
            NumEqVector storage = problem.solidDensity(element, scv)*problem.acceleration(
                element, scv, this->timeLoop().timeStepSize(), curVolVars.displacement()
            );
            storage *= curVolVars.extrusionFactor()*Extrusion::volume(fvGeometry, scv);
            residual[scv.localDofIndex()] += storage;
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Problem class must implement solidDensity and acceleration"
                " methods to evaluate storage term"
            );
    }

    /*!
     * \brief Evaluate the fluxes over a face of a sub control volume
     * flux: -div(P) where P is the 1st Piola-Kirchoff stress tensor (stress in reference configuration)
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // the deformation gradient F = grad(d) + I
        const auto F = [&]()
        {
            const auto& fluxVarCache = elemFluxVarsCache[scvf];
            Dune::FieldMatrix<Scalar, dimWorld, dimWorld> F(0.0);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                for (int dir = 0; dir < dimWorld; ++dir)
                    F[dir].axpy(volVars.displacement(dir), fluxVarCache.gradN(scv.indexInElement()));
            }

            for (int dir = 0; dir < dimWorld; ++dir)
                F[dir][dir] += 1;

            return F;
        }();

        // numerical flux ==> -P*n*dA
        NumEqVector flux(0.0);

        // problem implements the constitutive law P=P(F)
        const auto P = problem.firstPiolaKirchhoffStressTensor(F);
        P.mv(scvf.unitOuterNormal(), flux);
        flux *= -scvf.area();
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        // loading and body forces are implemented in the problem
        return ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);
    }
};

} // end namespace Dumux

#endif
