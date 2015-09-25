// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan                                        *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Quantities required by the one-phase box model defined on a vertex.
 */
#ifndef DUMUX_1P_VOLUME_VARIABLES_HH
#define DUMUX_1P_VOLUME_VARIABLES_HH

#include "1pproperties.hh"
#include <dumux/boxmodels/common/boxvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class OnePVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The local primary variable vector
     * \param problem The problem object
     * \param element The current element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars, problem, element, elemGeom, scvIdx, isOldSol);

        this->asImp_().updateTemperature_(priVars,
                                          element,
                                          elemGeom,
                                          scvIdx,
                                          problem);

        pressure_ = priVars[Indices::pressureIdx];
        density_ = Fluid::density(temperature(), pressure());
        viscosity_ = Fluid::viscosity(temperature(), pressure());

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
    };

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return pressure_; }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume.
     */
    Scalar viscosity() const
    { return viscosity_; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

protected:
    void updateTemperature_(const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem)
    { temperature_ = problem.temperature(element, elemGeom, scvIdx); }

    Scalar temperature_;
    Scalar pressure_;
    Scalar density_;
    Scalar viscosity_;
    Scalar porosity_;
};

}

#endif
