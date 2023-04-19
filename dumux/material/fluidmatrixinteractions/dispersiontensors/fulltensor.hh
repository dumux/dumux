// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::FullDispersionTensor
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_FULLTENSOR_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_FULLTENSOR_HH

#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Full dispersion tensor
 */
template<class TypeTag>
class FullDispersionTensor
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static const int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:

    template <class ElementFluxVariablesCache>
    static DimWorldMatrix compositionalDispersionTensor(const Problem& problem,
                                                        const SubControlVolumeFace& scvf,
                                                        const FVElementGeometry&,
                                                        const ElementVolumeVariables&,
                                                        const ElementFluxVariablesCache&,
                                                        const int phaseIdx,
                                                        const int compIdx)
    { return problem.spatialParams().dispersionTensor(scvf.center(), phaseIdx, compIdx); }

    template <class ElementFluxVariablesCache>
    static DimWorldMatrix thermalDispersionTensor(const Problem& problem,
                                                  const SubControlVolumeFace& scvf,
                                                  const FVElementGeometry&,
                                                  const ElementVolumeVariables&,
                                                  const ElementFluxVariablesCache&,
                                                  const int phaseIdx)
    { return problem.spatialParams().dispersionTensor(scvf.center(), phaseIdx); }

};

} // end namespace Dumux

#endif
