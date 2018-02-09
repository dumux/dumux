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
 * \ingroup ZeroEqModel
 *
 * \copydoc Dumux::ZEROEQVolumeVariables
 */
#ifndef DUMUX_ZEROEQ_VOLUME_VARIABLES_HH
#define DUMUX_ZEROEQ_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class ZeroEqVolumeVariablesImplementation;

/*!
 * \ingroup Reynolds-Averaged NavierStokesModel
 * \brief Volume variables for the single-phase Reynolds-Averaged Navier-Stokes model.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using ZeroEqVolumeVariables = ZeroEqVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup Reynolds-Averaged NavierStokesModel
 * \brief Volume variables for the isothermal single-phase Reynolds-Averaged Navier-Stokes model.
 */
template <class TypeTag>
class ZeroEqVolumeVariablesImplementation<TypeTag, false>
: public RANSVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = RANSVolumeVariablesImplementation<TypeTag, false>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    static const int defaultPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:
    ZeroEqVolumeVariablesImplementation()
    {
        eddyViscosityModel_ = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup),
                                                     "FreeFlow.EddyViscosityModel");
    }

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
    };


    /*!
     * \brief Calculate the eddy viscosity
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void calculateEddyViscosity(const ElementSolutionVector &elemSol,
                                const Problem &problem,
                                const Element &element,
                                const SubControlVolume& scv)
    {
        Scalar kinematicEddyViscosity = 0.0;
        if (eddyViscosityModel_ == Indices::prandtl)
        {
            Scalar mixingLength = 0.41 * asImp_().wallDistance_;
            Scalar velGrad = asImp_().velocityGradients()[0][1];
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        asImp_().setDynamicEddyViscosity(kinematicEddyViscosity * asImp_().density(defaultPhaseIdx));
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    int eddyViscosityModel_;
};

/*!
 * \ingroup ZeroEqModel
 * \brief Volume variables for the non-isothermal single-phase Reynolds-Averaged Navier-Stokes model.
 */
template <class TypeTag>
class ZeroEqVolumeVariablesImplementation<TypeTag, true>
: public RANSVolumeVariablesImplementation<TypeTag, true>
{ };
}

#endif
