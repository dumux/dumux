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
 * \ingroup ZeroEqModel
 * \copydoc Dumux::ZeroEqVolumeVariables
 */
#ifndef DUMUX_ZEROEQ_VOLUME_VARIABLES_HH
#define DUMUX_ZEROEQ_VOLUME_VARIABLES_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup ZeroEqModel
 * \brief Volume variables for the single-phase 0-Eq. model.
 */
template <class Traits, class NSVolumeVariables>
class ZeroEqVolumeVariables
: public RANSVolumeVariables< Traits, NSVolumeVariables>
{
    using RANSParentType = RANSVolumeVariables<Traits, NSVolumeVariables>;

    using Scalar = typename Traits::PrimaryVariables::value_type;

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
        additionalRoughnessLength_ = problem.additionalRoughnessLength(RANSParentType::elementIdx());
        yPlusRough_ = wallDistanceRough() * RANSParentType::uStar() / RANSParentType::kinematicViscosity();
        RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity(elemSol, problem, element, scv, problem.eddyViscosityModel()));
        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Calculate and set the dynamic eddy viscosity.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param modelName The name of the used model
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    Scalar calculateEddyViscosity(const ElementSolution &elemSol,
                                  const Problem &problem,
                                  const Element &element,
                                  const SubControlVolume& scv,
                                  const std::string modelName)
    {
        using std::abs;
        using std::exp;
        using std::sqrt;
        Scalar kinematicEddyViscosity = 0.0;
        unsigned int flowDirectionAxis = problem.flowDirectionAxis(RANSParentType::elementIdx());
        unsigned int wallNormalAxis = problem.wallNormalAxis(RANSParentType::elementIdx());
        Scalar velGrad = abs(RANSParentType::velocityGradients()[flowDirectionAxis][wallNormalAxis]);

        if (modelName.compare("none") == 0)
        {
            // kinematicEddyViscosity = 0.0
        }
        else if (modelName.compare("prandtl") == 0)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough();
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (modelName.compare("vanDriest") == 0)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough()
                                  * (1.0 - exp(-yPlusRough() / 26.0))
                                  / sqrt(1.0 - exp(-0.26 * yPlusRough()));
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (modelName.compare("baldwinLomax") == 0)
        {
            kinematicEddyViscosity = problem.kinematicEddyViscosity(RANSParentType::elementIdx());
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,
                       "The eddy viscosity model \"" << modelName << "\" is not implemented.");
        }

        return kinematicEddyViscosity * RANSParentType::density();
    }

    /*!
     * \brief Return the wall distance \f$\mathrm{[m]}\f$ including an additional roughness length
     */
    Scalar wallDistanceRough() const
    { return RANSParentType::wallDistance() + additionalRoughnessLength_; }

    /*!
     * \brief Return the dimensionless wall distance \f$\mathrm{[-]}\f$  including an additional roughness length
     */
    Scalar yPlusRough() const
    { return yPlusRough_; }

protected:
    Scalar additionalRoughnessLength_ = 0.0;
    Scalar yPlusRough_ = 0.0;
};

} // end namespace Dumux

#endif
