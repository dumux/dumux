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
 * \ingroup PorousmediumflowModels
 * \copydoc Dumux::PorousMediumFlowVolumeVariables
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_VOLUME_VARIABLES_HH
#define DUMUX_POROUSMEDIUMFLOW_VOLUME_VARIABLES_HH

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

namespace Dumux {
namespace Concepts {

template<typename T>
concept PMFlowModelTraits = requires {
    typename T::Indices;

    { T::numFluidPhases() };
    { T::numFluidComponents() };
};

template<typename T>
concept PMFlowVolVarsTraits = requires {
    typename T::PrimaryVariables;
    Indexable<typename T::PrimaryVariables;

    typename T::ModelTraits;
    PMFlowModelTraits<typename T::ModelTraits>;
};

} // namespace Concepts





/*!
 * \ingroup PorousmediumflowModels
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 * \tparam Traits The volume variables traits
 */
template<Concepts::PMFlowVolVarsTraits Traits>
class PorousMediumFlowVolumeVariables
{
public:
    using ModelTraits = typename Traits::ModelTraits;
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Scalar = IndexedType<PrimaryVariables>;

    //! Return number of phases considered by the model
    static constexpr int numFluidPhases()
    { return ModelTraits::numFluidPhases(); }

    //! Return number of components considered by the model
    static constexpr int numFluidComponents()
    { return ModelTraits::numFluidComponents(); }

    /*!
     * \brief Updates all quantities at a given integration point.
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<Concepts::Indexable ElemSol,
             typename Parameters,
             typename LocalGeometry,
             typename IntegrationPoint>
    void update(const ElemSol& elemSol,
                const Parameters& params,
                const LocalGeometry& localGeometry,
                const IntegrationPoint& ip)
    {
        // TODO: place priVars extraction somewhere to distinguish FEM/FV
        priVars_ = interpolate(elemSol, localGeom, ip);
        extrusionFactor_ = params.extrusionFactor(element, scv, elemSol);
    }

    /*!
     * \brief Returns the vector of primary variables.
     */
    const PrimaryVariables &priVars() const
    { return priVars_; }

    /*!
     * \brief Returns a component of primary variable vector.
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Returns how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
