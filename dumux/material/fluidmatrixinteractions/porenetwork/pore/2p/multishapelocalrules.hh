// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Implementation of capillary pressure curves for multiple pore body geometries
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_HH
#define DUMUX_PNM_2P_LOCAL_RULES_HH

#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include "localrulesforplatonicbody.hh"

namespace Dumux::PoreNetwork::FluidMatrix {
/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief LocalRulesTraits for implementation of capillary pressure curves for multiple pore body geometries
 */
template<class ScalarT>
struct LocalRulesTraits
{
    using Tetrahedron = TwoPLocalRulesPlatonicBodyDefault<Pore::Shape::tetrahedron, ScalarT>;
    using Cube = TwoPLocalRulesPlatonicBodyDefault<Pore::Shape::cube, ScalarT>;
    using Octahedron = TwoPLocalRulesPlatonicBodyDefault<Pore::Shape::octahedron, ScalarT>;
    using Icosahedron = TwoPLocalRulesPlatonicBodyDefault<Pore::Shape::icosahedron, ScalarT>;
    using Dodecahedron = TwoPLocalRulesPlatonicBodyDefault<Pore::Shape::dodecahedron, ScalarT>;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Implementation of capillary pressure curves for multiple pore body geometries
 */
template<class ScalarT>
class MultiShapeTwoPLocalRules : public Dumux::FluidMatrix::Adapter<MultiShapeTwoPLocalRules<ScalarT>, Dumux::FluidMatrix::PcKrSw>
{
public:
    using Scalar = ScalarT;
    using LocalRules = LocalRulesTraits<Scalar>;

    struct BasicParams
    {
        BasicParams() {}

        template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
        BasicParams(const SpatialParams& spatialParams,
                    const Element& element,
                    const SubControlVolume& scv,
                    const ElemSol& elemSol)
        {
            shape_ = spatialParams.gridGeometry().poreGeometry(scv.dofIndex());

            switch (shape_)
            {
                case Pore::Shape::tetrahedron:
                case Pore::Shape::cube:
                case Pore::Shape::octahedron:
                case Pore::Shape::icosahedron:
                case Pore::Shape::dodecahedron:
                    platonicBodyParams_ = std::make_unique<PlatonicBodyParams<Scalar>>(spatialParams, element, scv, elemSol);
                    break;
                default:
                    DUNE_THROW(Dune::NotImplemented, "Invalid shape");
            }
        }

         // we use cube specialization here, but all platonic bodies have the same params
        PlatonicBodyParams<Scalar> platonicBodyParams() const
        { return *platonicBodyParams_; }

        Pore::Shape poreShape() const
        { return shape_; }

        void setParams(const PlatonicBodyParams<Scalar>& platonicBodyParams)
        {
            shape_ = platonicBodyParams.poreShape();
            platonicBodyParams_ = std::make_unique<PlatonicBodyParams<Scalar>>(platonicBodyParams);
        }

    private:
        std::unique_ptr<PlatonicBodyParams<Scalar>> platonicBodyParams_;
        Pore::Shape shape_;
    };

    static constexpr bool supportsMultipleGeometries()
    { return true; }

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    MultiShapeTwoPLocalRules(const BasicParams& baseParams,
                             //const RegularizationParams& regParams = {}, TODO
                             const std::string& paramGroup = "")
    {
        shape_ = baseParams.poreShape();
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                localRulesForTetrahedron_ = std::make_shared<typename LocalRules::Tetrahedron>(baseParams.platonicBodyParams(),
                                                                                               typename LocalRules::Tetrahedron::RegularizationParams{},
                                                                                               paramGroup);
                break;
            case Pore::Shape::cube:
                localRulesForCube_ = std::make_shared<typename LocalRules::Cube>(baseParams.platonicBodyParams(),
                                                                                 typename LocalRules::Cube::RegularizationParams{},
                                                                                 paramGroup);
                break;
            case Pore::Shape::octahedron:
                localRulesForOctahedron_ = std::make_shared<typename LocalRules::Octahedron>(baseParams.platonicBodyParams(),
                                                                                             typename LocalRules::Octahedron::RegularizationParams{},
                                                                                             paramGroup);
                break;
            case Pore::Shape::icosahedron:
                localRulesForIcosahedron_ = std::make_shared<typename LocalRules::Icosahedron>(baseParams.platonicBodyParams(),
                                                                                               typename LocalRules::Icosahedron::RegularizationParams{},
                                                                                               paramGroup);
                break;
            case Pore::Shape::dodecahedron:
                localRulesForDodecahedron_ = std::make_shared<typename LocalRules::Dodecahedron>(baseParams.platonicBodyParams(),
                                                                                                 typename LocalRules::Dodecahedron::RegularizationParams{},
                                                                                                 paramGroup);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }


    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void updateParams(const SpatialParams& spatialParams,
                      const Element& element,
                      const SubControlVolume& scv,
                      const ElemSol& elemSol)
    {
        shape_ = spatialParams.gridGeometry().poreGeometry(scv.dofIndex());
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                return localRulesForTetrahedron_->updateParams(spatialParams, element, scv, elemSol);
            case Pore::Shape::cube:
                return localRulesForCube_->updateParams(spatialParams, element, scv, elemSol);
            case Pore::Shape::octahedron:
                return localRulesForOctahedron_->updateParams(spatialParams, element, scv, elemSol);
            case Pore::Shape::icosahedron:
                return localRulesForIcosahedron_->updateParams(spatialParams, element, scv, elemSol);
            case Pore::Shape::dodecahedron:
                return localRulesForDodecahedron_->updateParams(spatialParams, element, scv, elemSol);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    Scalar pc(const Scalar sw) const
    {
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                return localRulesForTetrahedron_->pc(sw);
            case Pore::Shape::cube:
                return localRulesForCube_->pc(sw);
            case Pore::Shape::octahedron:
                return localRulesForOctahedron_->pc(sw);
            case Pore::Shape::icosahedron:
                return localRulesForIcosahedron_->pc(sw);
            case Pore::Shape::dodecahedron:
                return localRulesForDodecahedron_->pc(sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The saturation-capilllary-pressure curve
     */
    Scalar sw(const Scalar pc) const
    {
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                return localRulesForTetrahedron_->sw(pc);
            case Pore::Shape::cube:
                return localRulesForCube_->sw(pc);
            case Pore::Shape::octahedron:
                return localRulesForOctahedron_->sw(pc);
            case Pore::Shape::icosahedron:
                return localRulesForIcosahedron_->sw(pc);
            case Pore::Shape::dodecahedron:
                return localRulesForDodecahedron_->sw(pc);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    Scalar dpc_dsw(const Scalar sw) const
    {
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                return localRulesForTetrahedron_->dpc_dsw(sw);
            case Pore::Shape::cube:
                return localRulesForCube_->dpc_dsw(sw);
            case Pore::Shape::octahedron:
                return localRulesForOctahedron_->dpc_dsw(sw);
            case Pore::Shape::icosahedron:
                return localRulesForIcosahedron_->dpc_dsw(sw);
            case Pore::Shape::dodecahedron:
                return localRulesForDodecahedron_->dpc_dsw(sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    Scalar dsw_dpc(const Scalar pc) const
    {
        switch (shape_)
        {
            case Pore::Shape::tetrahedron:
                return localRulesForTetrahedron_->dsw_dpc(pc);
            case Pore::Shape::cube:
                return localRulesForCube_->dsw_dpc(pc);
            case Pore::Shape::octahedron:
                return localRulesForOctahedron_->dsw_dpc(pc);
            case Pore::Shape::icosahedron:
                return localRulesForIcosahedron_->dsw_dpc(pc);
            case Pore::Shape::dodecahedron:
                return localRulesForDodecahedron_->dsw_dpc(pc);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     * \note This is only for compatibility. Will not be used.
     */
    Scalar krw(const Scalar sw) const
    { return 1.0; }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     * \note This is only for compatibility. Will not be used.
     */
    Scalar dkrw_dsw(const Scalar sw) const
    { return 0; }

    /*!
     * \brief The relative permeability for the non-wetting phase
     * \note This is only for compatibility. Will not be used.
     */
    Scalar krn(const Scalar sw) const
    { return 1.0; }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     * \note This is only for compatibility. Will not be used.
     */
    Scalar dkrn_dsw(const Scalar sw) const
    { return 0.0; }

private:
    std::shared_ptr<typename LocalRules::Tetrahedron> localRulesForTetrahedron_;
    std::shared_ptr<typename LocalRules::Cube> localRulesForCube_;
    std::shared_ptr<typename LocalRules::Octahedron> localRulesForOctahedron_;
    std::shared_ptr<typename LocalRules::Icosahedron> localRulesForIcosahedron_;
    std::shared_ptr<typename LocalRules::Dodecahedron> localRulesForDodecahedron_;

    // TODO add more shapes here

    Pore::Shape shape_;
};

} // end namespace Dumux::FluidMatrix

#endif
