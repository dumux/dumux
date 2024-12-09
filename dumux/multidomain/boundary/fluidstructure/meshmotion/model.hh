// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_MULTIDOMAIN_FLUIDSTRUCTURE_MESH_MOTION_MODEL_HH
#define DUMUX_MULTIDOMAIN_FLUIDSTRUCTURE_MESH_MOTION_MODEL_HH

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/method.hh>

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/geomechanics/fvproblem.hh>
#include <dumux/geomechanics/elastic/model.hh>
#include <dumux/geomechanics/elastic/fvspatialparams.hh>

namespace Dumux::Properties::TTag {
struct BiharmonicMeshMotionModel {};
struct LinearElasticMeshMotionModel { using InheritsFrom = std::tuple<Elastic>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux {

/*!
 * \ingroup FluidStructure
 * \brief Assembling the biharmonic operator in mixed form:
 *   -∇⋅(∇x_i) = φ_i
 *   -∇⋅(∇φ_i) = 0
 * for each coordinate direction i.
 */
template<class TypeTag>
class BiharmonicMeshMotionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

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

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    { return NumEqVector(0.0); }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        // Compute ∇x_i, ∇φ_i
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        std::array<Dune::FieldVector<Scalar, dimWorld>, dimWorld> gradCoordinates = {};
        std::array<Dune::FieldVector<Scalar, dimWorld>, dimWorld> gradLaplacian = {};
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            for (int i = 0; i < dimWorld; ++i)
            {
                gradCoordinates[i].axpy(
                    volVars.priVar(i),
                    fluxVarCache.gradN(scv.indexInElement())
                );

                gradLaplacian[i].axpy(
                    volVars.priVar(i + dimWorld),
                    fluxVarCache.gradN(scv.indexInElement())
                );
            }
        }

        const auto alpha = problem.spatialParams().alpha(scvf.ipGlobal());
        const auto beta = problem.spatialParams().beta(scvf.ipGlobal());

        // Compute Fx_i = -∇x_i⋅n dA and Fw_i = -∇φ_i⋅n dA
        NumEqVector flux;
        for (int i = 0; i < dimWorld; ++i)
        {
            flux[i] = -alpha*(gradCoordinates[i]*scvf.unitOuterNormal())*scvf.area();
            flux[i + dimWorld] = -beta*(gradLaplacian[i]*scvf.unitOuterNormal())*scvf.area();
        }
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];
        for (int i = 0; i < dimWorld; ++i)
            source[i] = volVars.priVar(i + dimWorld);
        return source;
    }

};

/*!
 * \ingroup FluidStructure
 * \brief Problem definition for the biharmonic problem for mesh motion
 */
template<class TypeTag>
class BiharmonicMeshMotionProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dimWorld = GridView::dimensionworld;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
public:
    BiharmonicMeshMotionProblem(
        std::shared_ptr<const GridGeometry> gridGeometry,
        std::shared_ptr<const CouplingManager> couplingManager
    )
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        if (couplingManager_->isCoupledDof(scv.dofIndex(), CouplingManager::meshMotionIdx))
        {
            const auto d = couplingManager_->displacement(element, scv);
            PrimaryVariables values(0.0);
            for (int i = 0; i < dimWorld; ++i)
                values[i] = d[i];
            return values;
        }

        return PrimaryVariables(0.0);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setDirichlet(0, 2);
        values.setDirichlet(1, 3);
        values.setNeumann(0);
        values.setNeumann(1);
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

private:
    static constexpr Scalar eps_ = 3e-6;
    std::shared_ptr<const CouplingManager> couplingManager_;
};

template<class Scalar, class GridGeometry>
class BiharmonicSpatialParams
: public FVElasticSpatialParams<
    GridGeometry, Scalar, BiharmonicSpatialParams<Scalar, GridGeometry>
>
{
    using ThisType = BiharmonicSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    BiharmonicSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , alpha_(getParam<Scalar>("MeshMotion.BiharmonicAlpha", 1.0))
    , beta_(getParam<Scalar>("MeshMotion.BiharmonicBeta", 1.0))
    {}

    Scalar alpha(const GlobalPosition& globalPos) const
    { return alpha_; }

    Scalar beta(const GlobalPosition& globalPos) const
    { return beta_; }

private:
    Scalar alpha_, beta_;
};

/*!
 * \ingroup FluidStructure
 * \brief Problem definition for the linear elastic problem for mesh motion
 */
template<class TypeTag>
class LinearElasticMeshMotionProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
public:
    LinearElasticMeshMotionProblem(
        std::shared_ptr<const GridGeometry> gridGeometry,
        std::shared_ptr<const CouplingManager> couplingManager
    )
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        if (couplingManager_->isCoupledDof(scv.dofIndex(), CouplingManager::meshMotionIdx))
            return couplingManager_->displacement(element, scv);
        return PrimaryVariables(0.0);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

private:
    static constexpr Scalar eps_ = 3e-6;
    std::shared_ptr<const CouplingManager> couplingManager_;
};

template<class Scalar, class GridGeometry>
class LinearElasticSpatialParams
: public FVElasticSpatialParams<
    GridGeometry, Scalar, LinearElasticSpatialParams<Scalar, GridGeometry>
>
{
    using ThisType = LinearElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    LinearElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        lameParams_.setLambda(getParam<Scalar>("MeshMotion.LinearElasticLambda", 1.0));
        lameParams_.setMu(getParam<Scalar>("MeshMotion.LinearElasticMu", 1.0));
    }

    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    GlobalPosition gravity(const GlobalPosition& globalPos) const
    { return GlobalPosition(0.0); }

private:
    LameParams lameParams_;
};

} // end namespace Dumux

namespace Dumux::Properties {

//////////////////////////////////
// BiharmonicMeshMotionModel
//////////////////////////////////
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BiharmonicMeshMotionModel>
{ using type = BiharmonicMeshMotionModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::BiharmonicMeshMotionModel>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BiharmonicMeshMotionModel>
{
    static constexpr int dimWorld = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld;
    struct type
    {
        struct Indices
        {
            static constexpr int coordIdx(int dimIdx){ return dimIdx; };
            static constexpr int phiIdx(int dimIdx){ return dimWorld + dimIdx; };
        };

        static constexpr int numEq() { return 2*dimWorld; }
    };
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::BiharmonicMeshMotionModel>
{
    struct type
    {
        template <class OutputModule>
        static void initOutputModule(OutputModule& out) {}
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::BiharmonicMeshMotionModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::BiharmonicMeshMotionModel>
{
    struct Traits
    {
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };
    using type = BasicVolumeVariables<Traits>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::BiharmonicMeshMotionModel>
{ using type = BiharmonicMeshMotionProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BiharmonicMeshMotionModel>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = BiharmonicSpatialParams<Scalar, GridGeometry>;
};

//////////////////////////////////
// LinearElasticMeshMotionModel
//////////////////////////////////
template<class TypeTag>
struct Problem<TypeTag, TTag::LinearElasticMeshMotionModel>
{ using type = LinearElasticMeshMotionProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LinearElasticMeshMotionModel>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = LinearElasticSpatialParams<Scalar, GridGeometry>;
};

} // end namespace Dumux::Properties

#endif
