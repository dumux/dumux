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
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::FullDispersionTensor
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_FULLTENSOR_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_FULLTENSOR_HH

namespace Dumux {

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
    static DimWorldMatrix compositionalDispersionTensor([[maybe_unused]] const Problem& problem,
                                                        [[maybe_unused]] const SubControlVolumeFace& scvf,
                                                        [[maybe_unused]] const FVElementGeometry& fvGeometry,
                                                        [[maybe_unused]] const ElementVolumeVariables& elemVolVars,
                                                        [[maybe_unused]] const ElementFluxVariablesCache& elemFluxVarsCache,
                                                        [[maybe_unused]] const int phaseIdx,
                                                        [[maybe_unused]] const int compIdx)
    { return problem.spatialParams().dispersionTensor(scvf.center(), phaseIdx, compIdx); }

};

} // end namespace Dumux

#endif
