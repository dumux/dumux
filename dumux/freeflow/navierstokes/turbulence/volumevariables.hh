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
 * \ingroup RANSModel
 *
 * \copydoc Dumux::RANSVolumeVariables
 */
#ifndef DUMUX_RANS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_VOLUME_VARIABLES_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class Traits, class NSVolumeVariables, int dim>
class RANSVolumeVariables
: public NavierStokesMassOnePVolumeVariables<Traits>
{
    using NavierStokesParentType = NavierStokesMassOnePVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using VelocityGradientMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();
public:
    static constexpr bool usesTKETurbulenceModel = Traits::ModelTraits::numTurbulenceEqs() > 0;

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
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        NavierStokesParentType::update(elemSol, problem, element, scv);
    }

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void updateRANSVariables(const ElementSolution& elemSol,
                             const Problem& problem,
                             const Element& element,
                             const SubControlVolume& scv)
    {
        wallDistance_ = problem.wallDistance(element);
    }

    /*!
     * \brief Calculates the eddy diffusivity \f$\mathrm{[m^2/s]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent Schmidt number
     */
    void calculateEddyDiffusivity()
    { eddyDiffusivity_ = kinematicEddyViscosity() / turbulentSchmidtNumber(); }

    /*!
     * \brief Calculates the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent Prandtl number
     */
    void calculateEddyThermalConductivity()
    {
        if constexpr (enableEnergyBalance)
        {
            eddyThermalConductivity_ = kinematicEddyViscosity() * NavierStokesParentType::density()
                                     * NavierStokesParentType::heatCapacity() / turbulentPrandtlNumber();
        }
        else
            eddyThermalConductivity_ = 0.0;
    }

    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    Scalar effectiveViscosity() const
    { return NavierStokesParentType::viscosity() + dynamicEddyViscosity(); }

    Scalar kinematicEddyViscosity() const
    { return dynamicEddyViscosity() / NavierStokesParentType::density(); }

    Scalar turbulentPrandtlNumber() const
    { return 0.7; }

    Scalar turbulentSchmidtNumber() const
    { return 0.6; }

    Scalar karmanConstant() const
    { return 0.41; }

    Scalar wallDistance() const
    { return wallDistance_; }

protected:

    Scalar stressTensorScalarProduct_;
    Scalar vorticityTensorScalarProduct_;
    VelocityGradientMatrix ccVelocityGradients_;

    Scalar dynamicEddyViscosity_ = 0.0; // With no turbulence model, this remains zero.
    Scalar eddyDiffusivity_ = 0.0;
    Scalar eddyThermalConductivity_ = 0.0;

    Scalar wallDistance_;
};
} // end namespace Dumux

#endif
