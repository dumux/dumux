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
 * \ingroup PorousmediumFlow
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_VOLUME_VARIABLES_HH
#define DUMUX_POROUSMEDIUMFLOW_VOLUME_VARIABLES_HH

namespace Dumux {

//! helper struct detecting if the volVars allow for a special boundary treatment
template<class ElemSol, class Problem, class Element, class Scv, class Scvf>
struct supportsBoundaryUpdates
{
    template<class VolumeVariables>
    auto operator()(VolumeVariables&& v)
    -> decltype(v.update(std::declval<ElemSol>(),
                         std::declval<Problem>(),
                         std::declval<Element>(),
                         std::declval<Scv>(),
                         std::declval<Scvf>()))
    {}
};

/*!
 * \ingroup PorousmediumFlow
 * \brief The isothermal base class
 *
 * \tparam Traits The volume variables traits
 * \tparam Impl The implementation of the volume variables
 */
template<class Traits>
class PorousMediumFlowVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the type encapsulating primary variable indices
    using Indices = typename Traits::ModelTraits::Indices;

    //! return number of phases considered by the model
    static constexpr int numPhases() { return Traits::ModelTraits::numPhases(); }
    //! return number of components considered by the model
    static constexpr int numComponents() { return Traits::ModelTraits::numComponents(); }

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        priVars_ = elemSol[scv.localDofIndex()];
        extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &priVars() const
    { return priVars_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    template<class VolumeVariables, class ElemSol, class Problem, class Element, class Scv, class Scvf>
    static auto updateBoundaryVolVars_(VolumeVariables& volVars,
                                       const ElemSol& elemSol,
                                       const Problem& problem,
                                       const Element& element,
                                       const Scv& scv,
                                       const Scvf& scvf)
     -> typename std::enable_if_t<decltype(isValid(supportsBoundaryUpdates<ElemSol, Problem, Element, Scv, Scvf>())
                                           (volVars))::value, void>
     {
         volVars.update(elemSol,
                        problem,
                        element,
                        scv,
                        scvf);
     }

     template<class VolumeVariables, class ElemSol, class Problem, class Element, class Scv, class Scvf>
     static auto updateBoundaryVolVars_(VolumeVariables& volVars,
                                       const ElemSol& elemSol,
                                       const Problem& problem,
                                       const Element& element,
                                       const Scv& scv,
                                       const Scvf& scvf)
     -> typename std::enable_if_t<!decltype(isValid(supportsBoundaryUpdates<ElemSol, Problem, Element, Scv, Scvf>())
                                            (volVars))::value, void>
     {
         volVars.update(elemSol,
                        problem,
                        element,
                        scv);
     }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
