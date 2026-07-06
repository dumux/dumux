// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Co-located momentum + Cahn-Hilliard CVFE model (Taylor-Hood P2 velocity / P2 c / P2 mu;
 *        pressure lives in a separate P1 subdomain). numEq = dim + 2.
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_MODEL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/localresidual.hh>
#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/volumevariables.hh>
#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/indices.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the co-located momentum + Cahn-Hilliard CVFE model.
 * \tparam dimension The spatial dimension
 */
template<int dimension>
struct NavierStokesMomentumCHNSCVFEModelTraits
{
    static constexpr int dim() { return dimension; }

    //! dim momentum-balance equations + phase-field transport + chemical-potential definition
    static constexpr int numEq() { return dim() + 2; }

    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return 1; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }

    using Indices = NavierStokesMomentumCHNSCVFEIndices<dim()>;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the co-located momentum + Cahn-Hilliard volume variables.
 */
template<class PV, class MT>
struct NavierStokesMomentumCHNSCVFEVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
//! The type tag for the co-located momentum + Cahn-Hilliard CVFE model
struct NavierStokesMomentumCHNSCVFE { using InheritsFrom = std::tuple<FreeFlow>; };
} // end namespace TTag

//! model traits (numEq = dim + 2)
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMomentumCHNSCVFE>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr auto dim = GridView::dimension;
public:
    using type = NavierStokesMomentumCHNSCVFEModelTraits<dim>;
};

//! the local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMomentumCHNSCVFE>
{ using type = NavierStokesMomentumCHNSCVFELocalResidual<TypeTag>; };

//! the volume variables
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMomentumCHNSCVFE>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = NavierStokesMomentumCHNSCVFEVolumeVariablesTraits<PV, MT>;
public:
    using type = NavierStokesMomentumCHNSCVFEVolumeVariables<Traits>;
};

//! default: not coupled (the test overrides with a pressure<->momentum coupling manager)
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::NavierStokesMomentumCHNSCVFE>
{
    struct EmptyCouplingManager {};
    using type = EmptyCouplingManager;
};

} // end namespace Dumux::Properties

#endif
