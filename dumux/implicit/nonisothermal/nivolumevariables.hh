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
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase two-component
 *        model.
 */
#ifndef DUMUX_NI_VOLUME_VARIABLES_HH
#define DUMUX_NI_VOLUME_VARIABLES_HH

#include "niproperties.hh"

namespace Dumux
{

/*!
 * \ingroup NIModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase two-component
 *        model.
 */
template <class TypeTag>
class NIVolumeVariables : public GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables)
{
    //! \cond false
    typedef typename GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables) ParentType;
    typedef typename ParentType::FluidState FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    //! \endcond

public:
    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(const int phaseIdx) const
    { return this->fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(const int phaseIdx) const
    { return this->fluidState_.enthalpy(phaseIdx); };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(kg K)]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidHeatCapacity_; };

    Scalar heatCapacity() const
    DUNE_DEPRECATED_MSG("use solidHeatCapacity()*solidDensity()*(1 - porosity()) instead")
    { return solidHeatCapacity()*solidDensity()*(1 - this->porosity()); }

    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    { return solidDensity_; };

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of a fluid phase in
     *        the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(this->fluidState_, phaseIdx); };

    Scalar thermalConductivityFluid(const int phaseIdx) const
    DUNE_DEPRECATED_MSG("use fluidThermalConductivity() instead")
    { return fluidThermalConductivity(phaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of the solid phase in
     *        the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return solidThermalConductivity_; };

    Scalar thermalConductivitySolid() const
    DUNE_DEPRECATED_MSG("use solidThermalConductivity() instead")
    { return solidThermalConductivity(); }

protected:
    // The methods below get called by the parent class. Since they
    // are protected, we are friends with our parent.
#if __clang__ || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    friend typename GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables);
#else
    friend class GET_PROP_TYPE(TypeTag, IsothermalVolumeVariables);
#endif

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               const int scvIdx)
    {
        return priVars[temperatureIdx];
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            const int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param sol The solution primary variables
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    {
        solidHeatCapacity_ = problem.spatialParams().solidHeatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidHeatCapacity_);

        solidDensity_ = problem.spatialParams().solidDensity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidDensity_);

        solidThermalConductivity_
            = problem.spatialParams().solidThermalConductivity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidThermalConductivity_);
    };

    Scalar solidHeatCapacity_;
    Scalar solidDensity_;
    Scalar solidThermalConductivity_;
};

} // end namespace

#endif
