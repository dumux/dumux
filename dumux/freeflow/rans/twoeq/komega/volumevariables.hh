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
 * \ingroup KOmegaEqModel
 *
 * \copydoc Dumux::KOmegaVolumeVariables
 */
#ifndef DUMUX_KOMEGA_VOLUME_VARIABLES_HH
#define DUMUX_KOMEGA_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

#include "models.hh"

namespace Dumux
{
<<<<<<< HEAD

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class KOmegaVolumeVariablesImplementation;

/*!
 * \ingroup KOmegaEqModel
 * \brief Volume variables for the single-phase k-omega 2-Eq model
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using KOmegaVolumeVariables = KOmegaVolumeVariablesImplementation<TypeTag, GET_PROP_TYPE(TypeTag, ModelTraits)::enableEnergyBalance()>;
=======
>>>>>>> f6e896abb... First set of K-omega classes

/*!
 * \ingroup KOmegaEqModel
 * \brief Volume variables for the isothermal single-phase k-omega 2-Eq model.
 */
template <class TypeTag>
class KOmegaVolumeVariablesImplementation<TypeTag, false>
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
     * \param problem The object specifying the problem which ought to be simulated
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
        kOmegaModel_ = problem.kOmegaModel();
        betaOmega_ = problem.betaOmega();
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];
        storedDissipation_ = problem.storedDissipation_[ParentType::elementID()];
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy_[ParentType::elementID()];
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct_[ParentType::elementID()];
        dofPosition_ = scv.dofPosition();
        if (problem.useStoredEddyViscosity_)
            ParentType::setKinematicEddyViscosity(problem.storedKinematicEddyViscosity_[ParentType::elementID()]);
        else
            calculateEddyViscosity();
    }

    /*!
     * \brief Calculate and set the dynamic eddy viscosity.
     */
    void calculateEddyViscosity()
    {
        using std::sqrt;
        using std::max;
        Scalar kinematicEddyViscosity = 0.0;
        if(kOmegaModel_ == KOmegaModels::wilcox88)
            kinematicEddyViscosity = turbulentKineticEnergy() / dissipation();
        else if (kOmegaModel_ == KOmegaModels::wilcox08)
          {
            Scalar limitiedDissipation = (7.0 / 8.0) * std::sqrt( 2.0 * stressTensorScalarProduct() * stressTensorScalarProduct() / betaK() );
            Scalar Wbar = std::max(dissipation(), limitiedDissipation);
            kinematicEddyViscosity = turbulentKineticEnergy() / Wbar;
          }
        else
              DUNE_THROW(Dune::NotImplemented, "This model is not implemented.");
        ParentType::setKinematicEddyViscosity(kinematicEddyViscosity);
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
    Scalar dissipation() const
    {
        return dissipation_;
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
    Scalar storedDissipation() const
    {
        return storedDissipation_;
    }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    {
        return stressTensorScalarProduct_;
    }

    /**
     *  \brief Here is a list of the background coefficients that are called in the Komega Model
     */

    //! \brief Returns the \$f \alpha \$f value
    const Scalar alpha() const
    {
        if(kOmegaModel_ == KOmegaModels::wilcox88)
            return 0.55555555556;
        else if (kOmegaModel_ == KOmegaModels::wilcox08)
            return 0.520;
        else
            DUNE_THROW(Dune::NotImplemented, "This Model is not implemented.");
    }

    //! \brief Returns the \$f \sigma_k \$f constant
    const Scalar sigmaK() const
    {
        if(kOmegaModel_ == KOmegaModels::wilcox88)
            return 0.50;
        else if (kOmegaModel_ == KOmegaModels::wilcox08)
            return 0.60;
        else
            DUNE_THROW(Dune::NotImplemented, "This Model is not implemented.");
    }

    //! \brief Returns the \$f \sigma_{\omega} \$f constant
    const Scalar sigmaOmega() const
    {
        if(kOmegaModel_ == KOmegaModels::wilcox88)
            return 0.50;
        else if (kOmegaModel_ == KOmegaModels::wilcox08)
            return 0.50;
        else
            DUNE_THROW(Dune::NotImplemented, "This Model is not implemented.");
    }

    //! \brief Returns the \$f \beta_k \$f constant
    const Scalar betaK() const
    {
        if(kOmegaModel_ == KOmegaModels::wilcox88)
            return 0.09;
        else if (kOmegaModel_ == KOmegaModels::wilcox08)
            return 0.09;
        else
            DUNE_THROW(Dune::NotImplemented, "This Model is not implemented.");
    }
public:
    Scalar betaOmega_;

protected:
    Dune::FieldVector<Scalar,2> dofPosition_;
    int kOmegaModel_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    Scalar storedTurbulentKineticEnergy_;
    Scalar storedDissipation_;
    Scalar stressTensorScalarProduct_;
};

}

#endif
