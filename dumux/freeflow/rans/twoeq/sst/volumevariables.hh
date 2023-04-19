// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SSTModel
 *
 * \copydoc Dumux::SSTVolumeVariables
 */
#ifndef DUMUX_SST_VOLUME_VARIABLES_HH
#define DUMUX_SST_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup SSTModel
 * \brief Volume variables for the isothermal single-phase SST 2-Eq model.
 */
template <class Traits, class NSVolumeVariables>
class SSTVolumeVariables
:  public RANSVolumeVariables<Traits, NSVolumeVariables>
{
    using RANSParentType = RANSVolumeVariables<Traits, NSVolumeVariables>;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using DimVector = Dune::FieldVector<Scalar, Traits::ModelTraits::dim()>;

public:
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        RANSParentType::updateNavierStokesVolVars(elemSol, problem, element, scv);
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
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void updateRANSProperties(const ElementSolution &elemSol,
                              const Problem &problem,
                              const Element &element,
                              const SubControlVolume& scv)
    {
        RANSParentType::updateRANSProperties(elemSol, problem, element, scv);

        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];

        storedDissipation_ = problem.storedDissipation(RANSParentType::elementIdx());
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy(RANSParentType::elementIdx());

        storedDissipationGradient_ = problem.storedDissipationGradient(RANSParentType::elementIdx());
        storedTurbulentKineticEnergyGradient_ = problem.storedTurbulentKineticEnergyGradient(RANSParentType::elementIdx());

        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(RANSParentType::elementIdx());
        vorticityTensorScalarProduct_ = problem.vorticityTensorScalarProduct(RANSParentType::elementIdx());
        wallDistance_ = problem.wallDistance(RANSParentType::elementIdx());
        kinematicViscosity_ = problem.kinematicViscosity(RANSParentType::elementIdx());

        if (problem.useStoredEddyViscosity())
            RANSParentType::setDynamicEddyViscosity_(problem.storedDynamicEddyViscosity(RANSParentType::elementIdx()));
        else
            RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity(problem));

        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$ for the SST-model.
     */
    template<class Problem>
    Scalar calculateEddyViscosity(const Problem& problem)
    {
        using std::sqrt;
        using std::max;

        if(problem.sstModelVersion() == SSTModel::BSL)
            return turbulentKineticEnergy() / dissipation() * RANSParentType::density();
        else if(problem.sstModelVersion() == SSTModel::SST)
        {
            const Scalar dividend = a1SST() * turbulentKineticEnergy();
            const Scalar possibleMax1 = a1SST()*dissipation();
            const Scalar possibleMax2 = absoluteValueVorticity() * F2();
            const Scalar divisor = std::max(possibleMax1,possibleMax2);
            return dividend / divisor * RANSParentType::density();
        }
        else
            DUNE_THROW(Dune::NotImplemented, "\nThis SST Model is not implemented.\n");
    }


    //! \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
    Scalar turbulentKineticEnergy() const
    { return turbulentKineticEnergy_; }

    //! \brief Returns an effective dissipation \f$ m^2/s^3 \f$
    Scalar dissipation() const
    { return dissipation_; }

    //! \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
    Scalar storedTurbulentKineticEnergy() const
    { return storedTurbulentKineticEnergy_; }

    //! \brief Returns an effective dissipation \f$ m^2/s^3 \f$
    Scalar storedDissipation() const
    { return storedDissipation_; }

    //! \brief Returns the gradient of the turbulent kinetic energy \f$ m^2/s^2 \f$
    DimVector storedTurbulentKineticEnergyGradient() const
    { return storedTurbulentKineticEnergyGradient_; }

    //! \brief Returns the gradient of the effective dissipation \f$ m^2/s^3 \f$
    DimVector storedDissipationGradient() const
    { return storedDissipationGradient_; }

    //! \brief Returns the scalar product of the stress tensor
    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    //! \brief Returns the scalar product of the vorticity tensor
    Scalar vorticityTensorScalarProduct() const
    { return vorticityTensorScalarProduct_; }

    //! \brief Returns the kinematic viscosity
    Scalar kinematicViscosity() const
    { return kinematicViscosity_; }

    Scalar wallDistance() const
    { return wallDistance_; }

    //! \brief Returns the \f$ \beta_{\omega} \f$ constant
    const Scalar betaOmega() const
    { return 0.0708; }

    //! \brief Returns the \f$ \alpha \f$ value
    const Scalar alpha() const
    { return 0.52; }

    //! \brief Returns the \f$ \sigma_k \f$ constant
    const Scalar sigmaK() const
    { return 0.6; }

    //! \brief Returns the \f$ \sigma_{\omega} \f$ constant
    const Scalar sigmaOmega() const
    { return 0.5; }

    //! \brief Returns the \f$ \beta_k \f$ constant
    const Scalar betaK() const
    { return 0.09; }

    //! \brief Returns the \f$ \sigma_{k1,BSL} \f$ constant
    const Scalar sigmaK1BSL() const
    { return 0.5; }

    //! \brief Returns the \f$ \sigma_{k1,SST} \f$ constant
    const Scalar sigmaK1SST() const
    { return 0.85; }

    //! \brief Returns the \f$ \sigma_{k2} \f$ constant
    const Scalar sigmaK2() const
    { return 1.0; }

    //! \brief Returns the \f$ \sigma_{\omega1,BSL} \f$ constant
    const Scalar sigmaOmega1BSL() const
    { return 0.5; }

    //! \brief Returns the \f$ \sigma_{\omega1,SST} \f$ constant
    const Scalar sigmaOmega1SST() const
    { return 0.5; }

    //! \brief Returns the \f$ \sigma_{\omega2} \f$ constant
    const Scalar sigmaOmega2() const
    { return 0.856; }

    //! \brief Returns the \f$ \beta_{1,BSL} \f$ constant
    const Scalar beta1BSL() const
    { return 0.0750; }

    //! \brief Returns the \f$ \beta_{1,SST} \f$ constant
    const Scalar beta1SST() const
    { return 0.0750; }

    //! \brief Returns the \f$ \beta_{2} \f$ constant
    const Scalar beta2() const
    { return 0.0820; }

    //! \brief Returns the \f$ \beta^{*}_{1,BSL} \f$ constant
    const Scalar betaStar1BSL() const
    { return 0.09; }

    //! \brief Returns the \f$ \beta^{*}_{1,SST} \f$ constant
    const Scalar betaStar1SST() const
    { return 0.09; }

    //! \brief Returns the \f$ \beta^{*}_{2} \f$ constant
    const Scalar betaStar2() const
    { return 0.09; }

    //! \brief Returns the \f$ \kappa_{1,BSL} \f$ constant
    const Scalar kappa1BSL() const
    { return 0.41; }

    //! \brief Returns the \f$ \kappa_{1,SST} \f$ constant
    const Scalar kappa1SST() const
    { return 0.41; }

    //! \brief Returns the \f$ \kappa_{2} \f$ constant
    const Scalar kappa2() const
    { return 0.41; }

    //! \brief Returns the \f$ \gamma_{1,BSL} \f$ constant
    const Scalar gamma1BSL() const
    { return (beta1BSL()/betaStar1BSL()) - ((sigmaOmega1BSL()*kappa1BSL()*kappa1BSL())/(std::sqrt(betaStar1BSL()))); }

    //! \brief Returns the \f$ \gamma_{1,SST} \f$ constant
    const Scalar gamma1SST() const
    { return (beta1SST()/betaStar1SST()) - ((sigmaOmega1SST()*kappa1SST()*kappa1SST())/(std::sqrt(betaStar1SST()))); }

    //! \brief Returns the \f$ \gamma_{2} \f$ constant
    const Scalar gamma2() const
    { return (beta2()/betaStar2()) - ((sigmaOmega2()*kappa2()*kappa2())/(std::sqrt(betaStar2()))); }

    //! \brief Returns the \f$ a_{1,SST} \f$ constant
    const Scalar a1SST() const
    { return 0.31; }

    //! \brief Returns the absolute value of the vorticity\f$ \Omega \f$
    Scalar absoluteValueVorticity() const
    { return std::sqrt(2.0 * vorticityTensorScalarProduct()); }

    //! \brief Returns the transformation function \f$ F_{1} \f$  for the constants of the BSL- and SST-model
    Scalar F1() const
    {
        Scalar gradientProduct = 0.0;
        for (unsigned int i = 0; i < Traits::ModelTraits::dim(); ++i){
            gradientProduct += storedTurbulentKineticEnergyGradient()[i] * storedDissipationGradient()[i];}

        Scalar positiveCrossDiffusion = 2.0 * RANSParentType::density() * sigmaOmega2() / storedDissipation() * gradientProduct;
        positiveCrossDiffusion = std::max(positiveCrossDiffusion, 1e-20);

        const Scalar possibleMin2 = (4.0 * RANSParentType::density() * sigmaOmega2() * storedTurbulentKineticEnergy())
                                  / (positiveCrossDiffusion * wallDistance() * wallDistance());

        const Scalar possibleMax1 = (std::sqrt(storedTurbulentKineticEnergy())) / (0.09 * storedDissipation() * wallDistance());
        const Scalar possibleMax2 = (500.00 * kinematicViscosity()) / (wallDistance() * wallDistance() * storedDissipation());
        const Scalar possibleMin1 = std::max(possibleMax1, possibleMax2);

        const Scalar argument = std::min(possibleMin1,possibleMin2);

        return tanh(argument * argument * argument * argument);
    }


    //! \brief Returns the transformation function \f$ F_{2} \f$ for the eddy viscosity of the SST-model
    Scalar F2() const
    {
        const Scalar possibleMax1 = 2.0 * (std::sqrt(storedTurbulentKineticEnergy())) / (0.09 * storedDissipation() * wallDistance());
        const Scalar possibleMax2 = (500.0 * kinematicViscosity()) / (wallDistance() * wallDistance() * storedDissipation());
        const Scalar argument = std::max(possibleMax1, possibleMax2);

        return tanh(argument * argument);
    }

    //------------------------Constants for SST--------------------------------
    //! \brief Returns the \f$ \sigma_{k,SST} \f$ constant for the SST-model
    const Scalar sigmaKSST() const
    { return F1()*sigmaK1SST() + (1-F1())*sigmaK2(); }

    //! \brief Returns the \f$ \sigma_{\omega,SST} \f$ constant for the SST-model
    const Scalar sigmaOmegaSST() const
    { return F1()*sigmaOmega1SST() + (1-F1())*sigmaOmega2(); }

    //! \brief Returns the \f$ \beta_{SST} \f$ constant for the SST-model
    const Scalar betaSST() const
    { return F1()*beta1SST() + (1-F1())*beta2(); }

    //! \brief Returns the \f$ \beta^{*}_{SST} \f$ constant for the SST-model. \f$ \beta^{*} \f$ is the same for all models.
    const Scalar betaStarSST() const
    { return F1()*betaStar1SST() + (1-F1())*betaStar2(); }

    //! \brief Returns the \f$ \kappa_{SST} \f$ constant for the SST-model. \f$ \kappa \f$ is the same for all models.
    const Scalar kappaSST() const
    { return F1()*kappa1SST() + (1-F1())*kappa2(); }

    //! \brief Returns the \f$ \gamma_{SST} \f$ constant for the SST-model
    const Scalar gammaSST() const
    { return F1()*gamma1SST() + (1-F1())*gamma2(); }



    //------------------------Constants for BSL--------------------------------
    //! \brief Returns the \f$ \sigma_{k,BSL} \f$ constant for the BSL-model
    const Scalar sigmaKBSL() const
    { return F1()*sigmaK1BSL() + (1-F1())*sigmaK2(); }

    //! \brief Returns the \f$ \sigma_{\omega,BSL} \f$ constant for the BSL-model
    const Scalar sigmaOmegaBSL() const
    { return F1()*sigmaOmega1BSL() + (1-F1())*sigmaOmega2(); }

    //! \brief Returns the \f$ \beta_{BSL} \f$ constant for the BSL-model
    const Scalar betaBSL() const
    { return F1()*beta1BSL() + (1-F1())*beta2(); }

    //! \brief Returns the \f$ \beta^{*}_{BSL} \f$ constant for the BSL-model. \f$ \beta^{*} \f$ is the same for all models.
    const Scalar betaStarBSL() const
    { return F1()*betaStar1BSL() + (1-F1())*betaStar2(); }

    //! \brief Returns the \f$ \kappa_{BSL} \f$ constant for the BSL-model. \f$ \kappa \f$ is the same for all models.
    const Scalar kappaBSL() const
    { return F1()*kappa1BSL() + (1-F1())*kappa2(); }

    //! \brief Returns the \f$ \gamma_{BSL} \f$ constant for the BSL-model
    const Scalar gammaBSL() const
    { return F1()*gamma1BSL() + (1-F1())*gamma2(); }


protected:
    Scalar betaOmega_ = 0.0;

    Scalar dissipation_ = 0.0;
    Scalar storedDissipation_ = 0.0;
    DimVector storedDissipationGradient_ = DimVector(0.0);
    Scalar turbulentKineticEnergy_ = 0.0;
    Scalar storedTurbulentKineticEnergy_ = 0.0;
    DimVector storedTurbulentKineticEnergyGradient_ = DimVector(0.0);

    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar vorticityTensorScalarProduct_ = 0.0;
    Scalar wallDistance_ = 0.0;
    Scalar kinematicViscosity_ = 0.0;
};
}

#endif
