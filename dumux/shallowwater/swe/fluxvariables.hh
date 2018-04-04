// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup SweModel
 * \copydoc Dumux::SweFluxVariables
 */
#ifndef DUMUX_SWE_FLUXVARIABLES_HH
#define DUMUX_SWE_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablesbase.hh>

namespace Dumux
{
/*!
 * \ingroup SweModel
 * \brief The flux variables class for the shallow water model.
 *
 */
template<class TypeTag> //LEO TypeTag ist der Template parameter es ginge auch template<typename TypeTag>
class SweFluxVariables : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

//    using NumericalFluxType = typename GET_PROP_TYPE(TypeTag, NumericalFluxType);
//    using TurbulenceModelType = typename GET_PROP_TYPE(TypeTag, TurbulenceModelType);
    static constexpr bool enableNumericalFlux = GET_PROP_VALUE(TypeTag, EnableNumericalFlux);
    static constexpr bool enableTurbulenceModel = GET_PROP_VALUE(TypeTag, EnableTurbulenceModel);

public:

    //! The constructor
    SweFluxVariables()
    {
    }

    NumEqVector numericalFlux() const
    {

        NumEqVector fluxVector(0.0);
        if (enableNumericalFlux)
        {

//             return NumericalFluxType::flux(this->problem(),
//                                             this->element(),
//                                             this->fvGeometry(),
//                                             this->elemVolVars(),
//                                             this->scvFace(),
//                                             this->elemFluxVarsCache());
        }
        else
        {
            return fluxVector;
        }
    }

    NumEqVector turbulenceFlux() const
    {
        NumEqVector fluxVector(0.0);
        if (enableTurbulenceModel)
        {
//             return TurbulenceModelType::flux(this->problem(),
//                                             this->element(),
//                                             this->fvGeometry(),
//                                             this->elemVolVars(),
//                                             this->scvFace(),
//                                             this->elemFluxVarsCache());
        }
        else
        {
            return fluxVector;
        }
    }
};

} // end namespace

#endif
