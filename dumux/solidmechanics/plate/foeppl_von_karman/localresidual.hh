// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FoepplVonKarmanPlate
 * \brief Local residuals for the Föppl-von Kármán model
 */
#ifndef DUMUX_FOEPPL_VON_KARMAN_PLATE_LOCAL_RESIDUAL_HH
#define DUMUX_FOEPPL_VON_KARMAN_PLATE_LOCAL_RESIDUAL_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <type_traits>

#include <dumux/solidmechanics/plate/kirchhoff_love/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief In-plane stress resultant \f$ \mathbf{N} \f$ for the Föppl-von Kármán membrane strain
 *
 * \f[
 *   \boldsymbol{\varepsilon} = \tfrac{1}{2}\left(\nabla\mathbf{u} + (\nabla\mathbf{u})^T
 *                                                + \nabla w \otimes \nabla w\right),\quad
 *   \mathbf{N} = C\left\{(1-\nu)\boldsymbol{\varepsilon}
 *                        + \nu\operatorname{tr}(\boldsymbol{\varepsilon})\,\mathbf{I}\right\}
 * \f]
 *
 * An eigenstrain \f$ \boldsymbol{\varepsilon}_g \f$, the incompatible strain of growth or
 * thermal expansion, is subtracted from the membrane strain before the stress is formed.
 *
 * \param gradDisplacement in-plane displacement gradient, `gradDisplacement[i][j]` \f$ = \partial_j u_i \f$
 * \param gradDeformation gradient of the vertical deformation \f$ \nabla w \f$
 * \param membraneStiffness the extensional stiffness \f$ C = Et/(1-\nu^2) \f$
 * \param poissonRatio the Poisson ratio \f$ \nu \f$
 * \param eigenstrain the eigenstrain \f$ \boldsymbol{\varepsilon}_g \f$
 */
template<class GradDisplacement, class GradDeformation, class Scalar, class Eigenstrain>
auto vonKarmanStressResultant(const GradDisplacement& gradDisplacement,
                              const GradDeformation& gradDeformation,
                              Scalar membraneStiffness,
                              Scalar poissonRatio,
                              const Eigenstrain& eigenstrain)
{
    static constexpr int dim = GradDeformation::dimension;
    Dune::FieldMatrix<Scalar, dim, dim> strain(0.0);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            strain[i][j] = 0.5*(gradDisplacement[i][j] + gradDisplacement[j][i]
                                + gradDeformation[i]*gradDeformation[j]) - eigenstrain[i][j];

    Scalar trace(0.0);
    for (int i = 0; i < dim; ++i)
        trace += strain[i][i];

    auto N = strain;
    N *= membraneStiffness*(1.0 - poissonRatio);
    for (int i = 0; i < dim; ++i)
        N[i][i] += membraneStiffness*poissonRatio*trace;

    return N;
}

//! The stress resultant without eigenstrain
template<class GradDisplacement, class GradDeformation, class Scalar>
auto vonKarmanStressResultant(const GradDisplacement& gradDisplacement,
                              const GradDeformation& gradDeformation,
                              Scalar membraneStiffness,
                              Scalar poissonRatio)
{
    static constexpr int dim = GradDeformation::dimension;
    return vonKarmanStressResultant(gradDisplacement, gradDeformation, membraneStiffness, poissonRatio,
                                    Dune::FieldMatrix<Scalar, dim, dim>(0.0));
}

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief The eigenstrain of the in-plane problem at a position, zero unless the problem provides one
 */
template<class Problem, class GlobalPosition>
auto plateEigenstrain(const Problem& problem, const GlobalPosition& globalPos)
{
    using Scalar = std::decay_t<decltype(problem.membraneStiffness(globalPos))>;
    static constexpr int dim = GlobalPosition::dimension;
    if constexpr (requires { problem.eigenstrain(globalPos); })
        return Dune::FieldMatrix<Scalar, dim, dim>(problem.eigenstrain(globalPos));
    else
        return Dune::FieldMatrix<Scalar, dim, dim>(0.0);
}

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Local residual for the Föppl-von Kármán model (deformation and potentials)
 *
 * Identical to the Kirchhoff-Love deformation residual except that the membrane
 * traction \f$ \mathbf{N}\nabla w \f$ is added to the flux of the shear gradient
 * potential equation, which turns the linear plate into the Föppl-von Kármán plate.
 *
 * The transverse balance carries the inertia \f$ \rho h\,\ddot w \f$ of the plate if the
 * problem provides `massPerArea(globalPos)` and `acceleration(element, scv, dt, priVars)`,
 * the latter typically evaluated by a Newmark scheme the problem holds.
 */
template<class TypeTag>
class FoepplVonKarmanPlateLocalResidualDeformation
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
     * \brief Evaluate the storage contribution to the residual: the pseudo-time damping and,
     *        if the problem supplies a mass per area and an acceleration, the inertia of the
     *        transverse balance
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        ParentType::evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scv);

        const auto& volVars = curElemVolVars[scv];
        if constexpr (requires { problem.massPerArea(scv.dofPosition());
                                 problem.acceleration(element, scv, Scalar{}, volVars.priVars()); })
        {
            // the same sign as the damping, see computeStorage
            const auto acceleration = problem.acceleration(
                element, scv, this->timeLoop().timeStepSize(), volVars.priVars()
            );
            residual[scv.localDofIndex()][Indices::shearGradPotentialEqIdx]
                -= problem.massPerArea(scv.dofPosition())*acceleration[Indices::verticalDeformationIdx]
                   *volVars.extrusionFactor()*Extrusion::volume(fvGeometry, scv);
        }
    }

    /*!
     * \brief Evaluate the rate of change of all conserved quantities
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // Damping, for relaxing towards an equilibrium rather than solving for one
        // directly. Integrating \f$ c\,\dot{w} = -\partial E/\partial w \f$ in pseudo-time
        // takes the plate downhill in energy, so it can only come to rest at a stable
        // equilibrium; a Newton solve, which looks for stationary points, is equally
        // happy to sit on an unstable one.
        //
        // The damping belongs in the shear gradient potential equation, which is the
        // transverse force balance; the deformation equation is the kinematic constraint
        // \f$ \boldsymbol{\theta} = \nabla w \f$ and carries no force. It enters with a
        // MINUS sign: the residual of that equation is \f$ -\partial E/\partial w \f$,
        // because \f$ \nabla\cdot\nabla\varphi = -D\nabla^4 w \f$ and the load enters as
        // \f$ F = -p \f$.
        NumEqVector storage(0.0);
        storage[Indices::shearGradPotentialEqIdx]
            = -problem.damping(scv.dofPosition())*volVars.verticalDeformation();
        return storage;
    }

    /*!
     * \brief Evaluate the fluxes over a face of a sub control volume
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if (scvf.boundary())
            DUNE_THROW(Dune::InvalidStateException, "Calling computeFlux for boundary scvf");

        const auto& cm = problem.couplingManager();

        // rotate with J = [0 1, -1 0]
        const auto tangent = [&](){
            auto tangent = scvf.unitOuterNormal();
            std::swap(tangent[0], tangent[1]);
            tangent[1] = -tangent[1];
            return tangent;
        }();

        NumEqVector flux(0.0);
        for (const auto& qp : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& fluxVarCache = elemFluxVarsCache[qp.ipData()];
            Dune::FieldVector<Scalar, dimWorld> gradShearGradPotential(0.0);
            Dune::FieldVector<Scalar, dimWorld> gradDeformation(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto& volVars = elemVolVars[localDof];
                const auto& gradN = fluxVarCache.gradN(localDof.index());
                gradShearGradPotential.axpy(volVars.shearGradPotential(), gradN);
                gradDeformation.axpy(volVars.verticalDeformation(), gradN);
            }

            const auto rotation = cm.rotation(fvGeometry, qp.ipData());

            const auto& globalPos = qp.ipData().global();
            const auto N = vonKarmanStressResultant(
                cm.inPlaneDisplacementGradient(fvGeometry, qp.ipData()), gradDeformation,
                problem.membraneStiffness(globalPos), problem.poissonRatio(globalPos),
                plateEigenstrain(problem, globalPos)
            );

            Dune::FieldVector<Scalar, dimWorld> membraneTraction(0.0);
            N.mv(gradDeformation, membraneTraction);

            flux[Indices::shearGradPotentialEqIdx]
                += vtmv(scvf.unitOuterNormal(), 1.0, gradShearGradPotential+membraneTraction)*qp.weight();
            flux[Indices::deformationEqIdx]
                -= vtmv(scvf.unitOuterNormal(), 1.0, (gradDeformation-rotation))*qp.weight();
            flux[Indices::shearCurlPotentialEqIdx] -= vtmv(tangent, 1.0, rotation)*qp.weight();
        }
        return flux;
    }
};

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Local residual for the Föppl-von Kármán model (rotations)
 */
template<class TypeTag>
using FoepplVonKarmanPlateLocalResidualRotation = KirchhoffLovePlateLocalResidualRotation<TypeTag>;

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Local residual for the Föppl-von Kármán model (in-plane displacements)
 *
 * Implements \f$ -\nabla\cdot\mathbf{N} = \mathbf{f} \f$, where the stress resultant
 * \f$ \mathbf{N} \f$ carries the quadratic von Kármán contribution of the vertical
 * deformation and \f$ \mathbf{f} \f$ is an in-plane body force per unit area.
 */
template<class TypeTag>
class FoepplVonKarmanPlateLocalResidualInPlane
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
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;
public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    /*!
     * \brief Evaluate the rate of change of all conserved quantities
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // As in the deformation residual, but the in-plane force balance already has the
        // residual equal to \f$ +\partial E/\partial\mathbf{u} \f$, so the sign is the plain one.
        NumEqVector storage(0.0);
        const auto damping = problem.damping(scv.dofPosition());
        for (int dir = 0; dir < dimWorld; ++dir)
            storage[dir] = damping*volVars.displacement(dir);
        return storage;
    }

    /*!
     * \brief Evaluate the fluxes over a face of a sub control volume
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if (scvf.boundary())
            DUNE_THROW(Dune::InvalidStateException, "Calling computeFlux for boundary scvf");

        NumEqVector flux(0.0);
        for (const auto& qp : CVFE::quadratureRule(fvGeometry, scvf))
        {
            NumEqVector contribution(0.0);
            stressResultant(problem, fvGeometry, elemVolVars, qp.ipData(), elemFluxVarsCache)
                .mv(scvf.unitOuterNormal(), contribution);
            flux.axpy(-qp.weight(), contribution);
        }
        return flux;
    }

    /*!
     * \brief The in-plane stress resultant at a point of an element
     */
    template<class IpData>
    auto stressResultant(const Problem& problem,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const IpData& ipData,
                         const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[ipData];
        Dune::FieldMatrix<Scalar, dimWorld, dimWorld> gradDisplacement(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof];
            const auto& gradN = fluxVarCache.gradN(localDof.index());
            for (int dir = 0; dir < dimWorld; ++dir)
                gradDisplacement[dir].axpy(volVars.displacement(dir), gradN);
        }

        const auto& globalPos = ipData.global();
        return vonKarmanStressResultant(
            gradDisplacement, problem.couplingManager().verticalDeformationGradient(fvGeometry, ipData),
            problem.membraneStiffness(globalPos), problem.poissonRatio(globalPos),
            plateEigenstrain(problem, globalPos)
        );
    }
};

} // end namespace Dumux

#endif
