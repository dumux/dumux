// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MindlinReissnerPlate
 * \brief Local residual for the Mindlin-Reissner model
 */
#ifndef DUMUX_MINDLIN_REISSNER_PLATE_LOCAL_RESIDUAL_HH
#define DUMUX_MINDLIN_REISSNER_PLATE_LOCAL_RESIDUAL_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux {

/*!
 * \ingroup MindlinReissnerPlate
 * \brief Local residual for the Mindlin-Reissner model (deformation and potentials)
 */
template<class TypeTag>
class MindlinReissnerPlateLocalResidualDeformation
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
        // so far we only implement the equilibrium model
        return NumEqVector(0.0);
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

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradShearGradPotential(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradShearCurlPotential(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradDeformation(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof];
            const auto& gradN = fluxVarCache.gradN(localDof.index());
            gradShearGradPotential.axpy(volVars.shearGradPotential(), gradN);
            gradShearCurlPotential.axpy(volVars.shearCurlPotential(), gradN);
            gradDeformation.axpy(volVars.verticalDeformation(), gradN);
        }

        const auto rotation = problem.couplingManager().rotation(fvGeometry, scvf);

        // rotate with J = [0 1, -1 0]
        const auto tangent = [&](){
            auto tangent = scvf.unitOuterNormal();
            std::swap(tangent[0], tangent[1]);
            tangent[1] = -tangent[1];
            return tangent;
        }();

        const auto thickness = problem.thickness();
        const auto shearCorrectionFactor = problem.shearCorrectionFactor();
        const auto shearModulus = problem.shearModulus();
        std::swap(gradShearCurlPotential[0], gradShearCurlPotential[1]);
        gradShearCurlPotential[1] = -gradShearCurlPotential[1];
        gradShearCurlPotential /= shearModulus*thickness*shearCorrectionFactor;

        NumEqVector flux(0.0);
        flux[Indices::shearGradPotentialEqIdx] = 1.0*vtmv(scvf.unitOuterNormal(), 1.0, gradShearGradPotential)*scvf.area();
        gradShearGradPotential /= shearModulus*thickness*shearCorrectionFactor;

        flux[Indices::deformationEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), 1.0, gradDeformation-rotation-gradShearGradPotential)*scvf.area();
        flux[Indices::shearCurlPotentialEqIdx] = -1.0*vtmv(tangent, 1.0, rotation+gradShearCurlPotential)*scvf.area();
        return flux;
    }
};

/*!
 * \ingroup MindlinReissnerPlate
 * \brief Local residual for the Mindlin-Reissner model (rotations)
 */
template<class TypeTag>
class MindlinReissnerPlateLocalResidualRotation
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
        // so far we only implement the equilibrium equation
        return NumEqVector(0.0);
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

        const auto& cm =  problem.couplingManager();
        const auto vars = cm.deformationAndPotentials(fvGeometry, scvf);
        const auto psi = vars[cm.shearCurlPotentialIdx()];
        const auto phi = vars[cm.shearGradPotentialIdx()];

        Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Q(0.0);
        Q[0][0] = phi; Q[0][1] = psi;
        Q[1][0] = -psi;
        Q[1][1] = phi;

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldMatrix<Scalar, dimWorld, dimWorld> gradRotations(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof];
            const auto& gradN = fluxVarCache.gradN(localDof.index());
            for (int dir = 0; dir < dimWorld; ++dir)
                gradRotations[dir].axpy(volVars.rotation(dir), gradN);
        }

        const auto& globalPos = scvf.ipGlobal();
        const auto poissonRatio = problem.poissonRatio(globalPos);
        const auto D = problem.D(globalPos);
        Dune::FieldMatrix<Scalar, dimWorld, dimWorld> M(0.0);
        M[0][0] = -D*(gradRotations[0][0] + poissonRatio*gradRotations[1][1]);
        M[1][1] = -D*(gradRotations[1][1] + poissonRatio*gradRotations[0][0]);
        M[0][1] = -D*(1.0-poissonRatio)*0.5*(gradRotations[0][1] + gradRotations[1][0]);
        M[1][0] = M[0][1];
        M -= Q;

        NumEqVector flux(0.0);
        M.mv(scvf.unitOuterNormal(), flux);
        flux *= -scvf.area();
        return flux;
    }
};

} // end namespace Dumux

#endif
