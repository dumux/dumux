// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MembranePlate
 * \brief Local residual for the membrane plate model
 */
#ifndef DUMUX_MEMBRANE_PLATE_LOCAL_RESIDUAL_HH
#define DUMUX_MEMBRANE_PLATE_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux {

/*!
 * \ingroup MembranePlate
 * \brief Local residual for the membrane plate model
 *
 * Implements \f$ \nabla\cdot(T\,\nabla w) = F \f$ as a finite-volume flux:
 * the flux over a face is \f$ T\,(\nabla w \cdot \mathbf{n})\,|\Gamma| \f$.
 */
template<class TypeTag>
class MembranePlateLocalResidual
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

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // so far we only implement the equilibrium equation
        return NumEqVector(0.0);
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Scalar gradWn = 0.0;
        for (const auto& scv : scvs(fvGeometry))
            gradWn += elemVolVars[scv].deformation()
                      * (fluxVarCache.gradN(scv.indexInElement()) * scvf.unitOuterNormal());

        const auto tension = problem.spatialParams().tension(scvf.ipGlobal());

        NumEqVector flux(0.0);
        flux[Indices::deformationEqIdx] = tension * gradWn * scvf.area();
        return flux;
    }
};

} // end namespace Dumux

#endif
