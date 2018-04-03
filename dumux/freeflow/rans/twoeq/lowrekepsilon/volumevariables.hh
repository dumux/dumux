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
 * \ingroup LowReKEpsilonModel
 *
 * \copydoc Dumux::LowReKEpsilonVolumeVariables
 */
#ifndef DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH
#define DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

#include "models.hh"

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class LowReKEpsilonVolumeVariablesImplementation;

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Volume variables for the single-phase 0-Eq. model.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using LowReKEpsilonVolumeVariables = LowReKEpsilonVolumeVariablesImplementation<TypeTag, GET_PROP_TYPE(TypeTag, ModelTraits)::enableEnergyBalance()>;

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Volume variables for the isothermal single-phase 0-Eq. model.
 */
template <class TypeTag>
class LowReKEpsilonVolumeVariablesImplementation<TypeTag, false>
: virtual public RANSVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = RANSVolumeVariablesImplementation<TypeTag, false>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        updateRANSProperties(elemSol, problem, element, scv);
    }

    /*!
     * \brief Update all turbulent quantities for a given control volume
     *
     * Wall and roughness related quantities are stored. Eddy viscosity is set.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void updateRANSProperties(const ElementSolution &elemSol,
                              const Problem &problem,
                              const Element &element,
                              const SubControlVolume& scv)
    {
        ParentType::updateRANSProperties(elemSol, problem, element, scv);
        lowReKEpsilonModel_ = problem.lowReKEpsilonModel();
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipationTilde_ = elemSol[0][Indices::dissipationIdx];
        storedDissipationTilde_ = problem.storedDissipationTilde_[ParentType::elementID()];
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy_[ParentType::elementID()];
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct_[ParentType::elementID()];
        dofPosition_ = scv.dofPosition();
        if (getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.UseStoredEddyViscosity", true))
            ParentType::setDynamicEddyViscosity(problem.storedKinematicEddyViscosity_[ParentType::elementID()] * ParentType::density());
        else
            calculateEddyViscosity();
    }

    /*!
     * \brief Calculate and set the dynamic eddy viscosity.
     */
    void calculateEddyViscosity()
    {
        Scalar kinematicEddyViscosity
            = cMu() * fMu() * turbulentKineticEnergy() * turbulentKineticEnergy()
              / std::max(dissipationTilde(), getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.MinDissipation", -1));
        ParentType::setDynamicEddyViscosity(kinematicEddyViscosity * ParentType::density());
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar turbulentKineticEnergy() const
    {
        return turbulentKineticEnergy_;
    }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar dissipationTilde() const
    {
        return dissipationTilde_;
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar storedTurbulentKineticEnergy() const
    {
        return storedTurbulentKineticEnergy_;
    }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar storedDissipationTilde() const
    {
        return storedDissipationTilde_;
    }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    {
        return stressTensorScalarProduct_;
    }

    //! \brief Returns the \$f Re_\textrm{T} \$f value
    const Scalar reT() const
    {
        return turbulentKineticEnergy() * turbulentKineticEnergy()
               / ParentType::kinematicViscosity() / dissipationTilde();
    }

    //! \brief Returns the \$f Re_\textrm{y} \$f value
    const Scalar reY() const
    {
        using std::sqrt;
        return sqrt(turbulentKineticEnergy()) * ParentType::wallDistance()
               / ParentType::kinematicViscosity();
    }

    //! \brief Returns the \$f C_\mu \$f constant
    const Scalar cMu() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 0.09;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 0.09;
        else // use default model by Chien
            return 0.09;
    }

    //! \brief Returns the \$f \sigma_\textrm{k} \$f constant
    const Scalar sigmaK() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.0;
        else // use default model by Chien
            return 1.0;
    }

    //! \brief Returns the \$f \sigma_\varepsilon \$f constant
    const Scalar sigmaEpsilon() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.3;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.3;
        else // use default model by Chien
            return 1.3;
    }

    //! \brief Returns the \$f C_{1\tilde{\varepsilon}}  \$f constant
    const Scalar cOneEpsilon() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.44;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.44;
        else // use default model by Chien
            return 1.35;
    }

    //! \brief Returns the \$f C_{2\tilde{\varepsilon}} \$f constant
    const Scalar cTwoEpsilon() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.92;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.92;
        else // use default model by Chien
            return 1.8;
    }

    //! \brief Returns the \$f D \$f value
    const Scalar dValue() const
    {
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 0.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 0.0;
        else // use default model by Chien
            return 2.0 * ParentType::kinematicViscosity() * turbulentKineticEnergy()
                   / ParentType::wallDistance() / ParentType::wallDistance();
    }

    //! \brief Returns the \$f f_\mu \$f value
    const Scalar fMu() const
    {
        using std::exp;
        using std::pow;
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return pow((1.0 - exp(-0.0165 * reY())), 2.0)
                   * (1.0 + 20.5 / reT());
        else // use default model by Chien
            return 1.0 - exp(-0.0115 * ParentType::yPlus());
    }

    //! \brief Returns the \$f f_1 \$f value
    const Scalar fOne() const
    {
        using std::pow;
        if (lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.0 + pow((0.05 / fMu()), 3.0);
        else // use default model by Chien
            return 1.0;
    }

    //! \brief Returns the \$f f_2 \$f value
    const Scalar fTwo() const
    {
        using std::exp;
        if(lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 1.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 1.0 - exp(-1.0 * (reT() * reT()));
        else // use default model by Chien
            return 1.0 - 0.22 * exp(-1.0 * (reT() * reT() / 6.0 / 6.0));
    }

    //! \brief Returns the \$f E \$f value
    const Scalar eValue() const
    {
        using std::exp;
        if(lowReKEpsilonModel_ == LowReKEpsilonModels::standardHighReKEpsilon)
            return 0.0;
        else if (lowReKEpsilonModel_ == LowReKEpsilonModels::lamBremhorst)
            return 0.0;
        else // use default model by Chien
            return -2.0 * ParentType::kinematicViscosity() * dissipationTilde()
                   / ParentType::wallDistance() / ParentType::wallDistance()
                   * exp(-0.5 * ParentType::yPlus());
    }

protected:
    Dune::FieldVector<Scalar, 2> dofPosition_;
    int lowReKEpsilonModel_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipationTilde_;
    Scalar storedTurbulentKineticEnergy_;
    Scalar storedDissipationTilde_;
    Scalar stressTensorScalarProduct_;
};

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Volume variables for the non-isothermal single-phase 0-Eq. model.
 */
template <class TypeTag>
class LowReKEpsilonVolumeVariablesImplementation<TypeTag, true>
: virtual public LowReKEpsilonVolumeVariablesImplementation<TypeTag, false>,
  virtual public RANSVolumeVariablesImplementation<TypeTag, true>
{
    using ParentTypeNonIsothermal = RANSVolumeVariablesImplementation<TypeTag, true>;
    using ParentTypeIsothermal = LowReKEpsilonVolumeVariablesImplementation<TypeTag, false>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentTypeNonIsothermal::update(elemSol, problem, element, scv);
        updateRANSProperties(elemSol, problem, element, scv);
    }

    /*!
     * \brief Update all turbulent quantities for a given control volume
     *
     * Wall and roughness related quantities are stored. Eddy viscosity is set.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void updateRANSProperties(const ElementSolution &elemSol,
                              const Problem &problem,
                              const Element &element,
                              const SubControlVolume& scv)
    {
        ParentTypeIsothermal::updateRANSProperties(elemSol, problem, element, scv);
        ParentTypeNonIsothermal::calculateEddyThermalConductivity();
    }
};
}

#endif
