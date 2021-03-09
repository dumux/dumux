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
 * \ingroup OneEqModel
 *
 * \copydoc Dumux::OneEqVolumeVariables
 */
#ifndef DUMUX_ONEEQ_VOLUME_VARIABLES_HH
#define DUMUX_ONEEQ_VOLUME_VARIABLES_HH

#include <dune/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup OneEqModel
 * \brief Volume variables for the isothermal single-phase one-equation
 *        turbulence model by Spalart-Allmaras
 */
template <class Traits, class NSVolumeVariables>
class OneEqVolumeVariables
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
     * \param problem The object specifying the problem which ought to
     *                be simulated
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
        viscosityTilde_ = elemSol[0][Indices::viscosityTildeIdx];
        storedViscosityTilde_ = problem.storedViscosityTilde(RANSParentType::elementIdx());
        storedViscosityTildeGradient_ = problem.storedViscosityTildeGradient(RANSParentType::elementIdx());
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(RANSParentType::elementIdx());
        vorticityTensorScalarProduct_ = problem.vorticityTensorScalarProduct(RANSParentType::elementIdx());
        if (problem.useStoredEddyViscosity())
            RANSParentType::setDynamicEddyViscosity_(problem.storedDynamicEddyViscosity(RANSParentType::elementIdx()));
        else
            RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity());
        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$.
     */
    Scalar calculateEddyViscosity()
    { return viscosityTilde() * fv1() *  RANSParentType::density(); }

    /*!
     * \brief Returns the viscosity parameter \f$ m^2/s \f$
     */
    Scalar viscosityTilde() const
    { return viscosityTilde_; }

    /*!
     * \brief Returns the viscosity parameter from the last iteration \f$ m^2/s \f$
     */
    Scalar storedViscosityTilde() const
    { return storedViscosityTilde_; }

    /*!
     * \brief Returns the gradient of the viscosity parameter
     */
    DimVector storedViscosityTildeGradient() const
    { return storedViscosityTildeGradient_; }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    /*!
     * \brief Returns damping function for the eddy viscosity
     */
    Scalar fv1() const
    {
        return viscosityRatio() * viscosityRatio() * viscosityRatio()
               / (viscosityRatio() * viscosityRatio() * viscosityRatio()
                  + cv1() * cv1() * cv1());
    }

    //! \brief Returns a model function
    Scalar fv2() const
    { return 1.0 - viscosityRatio() / (1.0 + viscosityRatio() * fv1()); }

    //! \brief Returns a model function
    Scalar ft2() const
    {
        using std::exp;
        // the trip correction term is dropped according to Versteeg2009 and Wilcox2006
        // return ct3() * exp(-ct4() * viscosityRatio()  * viscosityRatio());
        return 0.0;
    }

    //! \brief Returns a model function
    Scalar fW() const
    {
        using std::pow;
        using Dune::power;
        return g() * pow(1.0 + power(cw3(), 6) / (power(g(), 6) + power(cw3(), 6)), 1.0/6.0);
    }

    //! \brief Returns a model function
    Scalar g() const
    {
        using Dune::power;
        return r() + cw2() * (power(r(), 6) - r());
    }

    //! \brief Returns a model function
    Scalar r() const
    {
        using std::min;
        return min(10.0,
                   viscosityTilde() / stressTensorScalarProductTilde()
                   / RANSParentType::karmanConstant() / RANSParentType::karmanConstant()
                   / RANSParentType::wallDistance() / RANSParentType::wallDistance());
    }

    //! \brief Returns the ratio of the kinematic viscosity and the viscosity parameter
    Scalar viscosityRatio() const
    { return viscosityTilde() / RANSParentType::kinematicViscosity(); }

    /*
     * ! \brief Returns a modified version of the stress tensor scalar product
     *
     * According to <a href="https://turbmodels.larc.nasa.gov/spalart.html">NASA</a>
     * this term should never be zero and different limiters might be used.
     * The Implementation uses the one proposed in:
     * Allmaras, S. R., Johnson, F. T., and Spalart, P. R.,
     * "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model," ICCFD7-1902
     */
    Scalar stressTensorScalarProductTilde() const
    {
        // original form
        // return vorticityMagnitude()
        //        + viscosityTilde() * fv2()
        //          / RANSParentType::karmanConstant() / RANSParentType::karmanConstant()
        //          / RANSParentType::wallDistance() / RANSParentType::wallDistance();

        // limiter form, literature source see above
        Scalar sBar = viscosityTilde() * fv2()
                      / RANSParentType::karmanConstant() / RANSParentType::karmanConstant()
                      / RANSParentType::wallDistance() / RANSParentType::wallDistance();
        return sBar < -c2() * vorticityMagnitude()
               ? vorticityMagnitude()
                 + (vorticityMagnitude() * (c2() * c2() * vorticityMagnitude() + c3() * sBar))
                   / ((c3() - 2.0 * c2()) * vorticityMagnitude() - sBar)
               : vorticityMagnitude() + sBar;
    }

    //! \brief Returns the magnitude of the vorticity
    Scalar vorticityMagnitude() const
    {
        using std::sqrt;
        return sqrt(2.0 * vorticityTensorScalarProduct_);
    }

    //! \brief Returns a model constant
    Scalar c2() const
    { return 0.7; }

    //! \brief Returns a model constant
    Scalar c3() const
    { return 0.9; }

    //! \brief Returns a model constant
    Scalar sigma() const
    { return 2.0/3.0; }

    //! \brief Returns a model constant
    Scalar cb1() const
    { return 0.1355; }

    //! \brief Returns a model constant
    Scalar cb2() const
    { return 0.622; }

    //! \brief Returns a model constant
    Scalar cv1() const
    { return 7.1; }

    //! \brief Returns a model constant
    Scalar ct3() const
    { return 1.2; }

    //! \brief Returns a model constant
    Scalar ct4() const
    { return 0.5; }

    //! \brief Returns a model constant
    Scalar cw1() const
    {
        return cb1() / RANSParentType::karmanConstant() / RANSParentType::karmanConstant()
               + (1.0 + cb2()) / sigma();
    }

    //! \brief Returns a model constant
    Scalar cw2() const
    { return 0.3; }

    //! \brief Returns a model constant
    Scalar cw3() const
    { return 2.0; }

protected:
    Scalar viscosityTilde_ = 0.0;
    Scalar storedViscosityTilde_ = 0.0;
    DimVector storedViscosityTildeGradient_ = DimVector(0.0);
    Scalar stressTensorScalarProduct_;
    Scalar vorticityTensorScalarProduct_;
};

} // end namespace Dumux

#endif
