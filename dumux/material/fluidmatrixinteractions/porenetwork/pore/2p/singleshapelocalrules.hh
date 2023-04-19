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
 * \brief Base classes for standard pore-local pc-Sw curves.
 */
#ifndef DUMUX_PNM_2P_SINGLE_SHAPE_LOCAL_RULES_HH
#define DUMUX_PNM_2P_SINGLE_SHAPE_LOCAL_RULES_HH

#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include <dumux/material/fluidmatrixinteractions/2p/noregularization.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::PoreNetwork::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Base class for all standard pore-local pc-Sw curves.
 */
template<class ScalarType,
         class BaseLaw,
         class Regularization = Dumux::FluidMatrix::NoRegularization>
class SingleShapeTwoPLocalRules : public Dumux::FluidMatrix::Adapter<SingleShapeTwoPLocalRules<ScalarType, BaseLaw, Regularization>, Dumux::FluidMatrix::PcKrSw>
{
public:

    using Scalar = ScalarType;

    using BasicParams = typename BaseLaw::template Params<Scalar>;
    using RegularizationParams = typename Regularization::template Params<Scalar>;

    static constexpr bool supportsMultipleGeometries()
    { return false; }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void updateParams(const SpatialParams& spatialParams,
                      const Element& element,
                      const SubControlVolume& scv,
                      const ElemSol& elemSol)
    {
        basicParams_.update(spatialParams, element, scv, elemSol);
        regularization_.updateParams(spatialParams, element, scv, elemSol);
    }

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    /*!
     * \brief Return whether this law is regularized
     */
    static constexpr bool isRegularized()
    { return !std::is_same<Regularization, Dumux::FluidMatrix::NoRegularization>::value; }

    /*!
     * \brief Construct from parameter structs
     * \note More efficient constructor but you need to ensure all parameters are initialized
     */
    SingleShapeTwoPLocalRules(const BasicParams& baseParams = {},
                              const RegularizationParams& regParams = {},
                              const std::string& paramGroup = "")
    : basicParams_(baseParams)

    {
        if constexpr (isRegularized())
            regularization_.init(this, regParams, paramGroup);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pc(const Scalar sw) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.pc(sw);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::pc(sw, basicParams_);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dpc_dsw(const Scalar sw) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dpc_dsw(sw);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::dpc_dsw(sw, basicParams_);
    }

    /*!
     * \brief The saturation-capillary-pressure curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar sw(const Scalar pc) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.sw(pc);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::sw(pc, basicParams_);
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dsw_dpc(const Scalar pc) const
    {
        if constexpr (enableRegularization)
        {
            const auto regularized = regularization_.dsw_dpc(pc);
            if (regularized)
                return regularized.value();
        }

        return BaseLaw::dsw_dpc(pc, basicParams_);
    }

    /*!
     * \brief The relative permeability for the wetting phase
     * \note This is only for compatibility. Will not be used.
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krw(const Scalar sw) const
    { return 1.0; }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     * \note This is only for compatibility. Will not be used.
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dkrw_dsw(const Scalar sw) const
    { return 0; }

    /*!
     * \brief The relative permeability for the non-wetting phase
     * \note This is only for compatibility. Will not be used.
     */
    template<bool enableRegularization = isRegularized()>
    Scalar krn(const Scalar sw) const
    { return 1.0; }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     * \note This is only for compatibility. Will not be used.
     */
    template<bool enableRegularization = isRegularized()>
    Scalar dkrn_dsw(const Scalar sw) const
    { return 0.0; }

    /*!
     * \brief Equality comparison with another instance
     */
    bool operator== (const SingleShapeTwoPLocalRules& o) const
    {
        return basicParams_ == o.basicParams_
               && regularization_ == o.regularization_;
    }

    /*!
     * \brief Return the base law's parameters
     */
    const BasicParams& basicParams() const
    { return basicParams_; }


private:
    BasicParams basicParams_;
    Regularization regularization_;
};

} // end namespace Dumux::PoreNetwork::FluidMatrix

#endif
