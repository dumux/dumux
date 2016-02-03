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
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 */

#ifndef DUMUX_2CNI_STOKES_2P2CNI_PROBLEM_HH
#define DUMUX_2CNI_STOKES_2P2CNI_PROBLEM_HH

#include <dumux/freeflow/stokesncni/properties.hh>
#include <dumux/multidomain/2cstokes2p2c/problem.hh>
#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 */
template <class TypeTag>
class TwoCNIStokesTwoPTwoCNIProblem : public TwoCStokesTwoPTwoCProblem<TypeTag>
{
    typedef TwoCStokesTwoPTwoCProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

public:
    //! The constructor
    template<class GridView>
    TwoCNIStokesTwoPTwoCNIProblem(TimeManager &timeManager,
                                  GridView gridView)
    : ParentType(timeManager, gridView)
    { }

    /*!
     * \brief Returns the temperature gradient through the boundary layer
     *
     * \param cParams a parameter container
     * \param scvIdx1 The local index of the sub-control volume of the Stokes domain
     */
    template<typename CParams>
    Scalar evalBoundaryLayerTemperatureGradient(CParams cParams, const int scvIdx) const
    {
        const Scalar temperatureOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
        Scalar normalTemperatureGrad = cParams.elemVolVarsCur1[scvIdx].temperature()
                                       - temperatureOut;
        return normalTemperatureGrad
               / asImp_().evalBoundaryLayerModel(cParams, scvIdx).thermalBoundaryLayerThickness();
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Dumux

#endif // DUMUX_2CNI_STOKES_2P2CNI_PROBLEM_HH
