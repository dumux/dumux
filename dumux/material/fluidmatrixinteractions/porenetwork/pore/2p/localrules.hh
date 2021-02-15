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
 *
 * \brief Implementation of capillary pressure curves for multiple pore body geometries
 *
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_HH
#define DUMUX_PNM_2P_LOCAL_RULES_HH

#include <dumux/porenetworkflow/common/poreproperties.hh>
#include "baselocalrules.hh"
#include "localrulesforcube.hh"

namespace Dumux::FluidMatrix {

template<class ScalarT, class LocalRulesForCube = TwoPLocalRulesCubeJoekarNiasarDefault<ScalarT>>
class MultiShapeTwoPLocalRules : public Adapter<MultiShapeTwoPLocalRules<ScalarT, LocalRulesForCube>, PcKrSw>
{
public:
    using Scalar = ScalarT;

    struct BasicParams
    {
        template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
        BasicParams(const SpatialParams& spatialParams,
                    const Element& element,
                    const SubControlVolume& scv,
                    const ElemSol& elemSol)
        {
            shape_ = spatialParams.gridGeometry().poreGeometry(scv.dofIndex());

            switch (shape_)
            {
                case Pore::Shape::cube:
                    cubeParams_ = std::make_unique<typename LocalRulesForCube::BasicParams>(spatialParams, element, scv, elemSol);
                    break;
                default:
                    DUNE_THROW(Dune::NotImplemented, "Invalid shape");
            }
        }

        typename LocalRulesForCube::BasicParams cubeParams() const
        { return *cubeParams_; }

        Pore::Shape poreShape() const
        { return shape_; }

    private:
        std::unique_ptr<typename LocalRulesForCube::BasicParams> cubeParams_;
        Pore::Shape shape_;
    };

    static constexpr bool supportsMultipleGeometries()
    { return true; }

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    static BasicParams makeParams(const SpatialParams& spatialParams,
                                  const Element& element,
                                  const SubControlVolume& scv,
                                  const ElemSol& elemSol)
    {
        return BasicParams(spatialParams, element, scv, elemSol);
    }

    MultiShapeTwoPLocalRules(const BasicParams& baseParams,
                             //const RegularizationParams& regParams = {}, TODO
                             const std::string& paramGroup = "")
    {
        shape_ = baseParams.poreShape();
        switch (shape_)
        {
            case Pore::Shape::cube:
                localRulesForCube_ = std::make_shared<LocalRulesForCube>(baseParams.cubeParams());
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
            case Pore::Shape::cube:
                return localRulesForCube_->updateParams(spatialParams, element, scv, elemSol);
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
            case Pore::Shape::cube:
                return localRulesForCube_->pc(sw);
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
            case Pore::Shape::cube:
                return localRulesForCube_->sw(pc);
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
    std::shared_ptr<LocalRulesForCube> localRulesForCube_;
    // TODO add more shapes here

    Pore::Shape shape_;
};

}

#endif
