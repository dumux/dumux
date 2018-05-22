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
 * \brief Quantities required by the linear elasticity box
 *        model defined on a vertex.
 */
#ifndef DUMUX_ELASTIC_VOLUME_VARIABLES_HH
#define DUMUX_ELASTIC_VOLUME_VARIABLES_HH

#include <dumux/implicit/volumevariables.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ElasticBoxModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the linear elasticity model.
 */
template <class TypeTag>
class ElasticVolumeVariablesBase : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum{
        dim = GridView::dimension,
    };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar,dim> DimVector;

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {   std::cout << "KAMEN HIER VORBEI volvars" << std::endl;
        ParentType::update(priVars,
                           problem,
                           element,
                           scv);

        primaryVars_ = priVars;

        for (int i = 0; i < dim; ++i)
            displacement_[i] = priVars[Indices::u(i)];

        // retrieve Lame parameters and rock density from spatialParams
        const Dune::FieldVector<Scalar, 2> &lameParams = problem.spatialParams().lameParams(element, scv);
        lambda_ = lameParams[0];
        mu_ = lameParams[1];
        rockDensity_ = problem.spatialParams().rockDensity(element, scv);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    /*!
     * \brief Sets the evaluation point used in the by the local jacobian.
     */
    void setEvalPoint(const Implementation *ep)
    { }

    /*!
      * \brief Return the Lame parameter lambda \f$\mathrm{[Pa]}\f$ within the control volume.
      */
    Scalar lambda() const
    { return lambda_; }

    /*!
      * \brief Return the Lame parameter mu \f$\mathrm{[Pa]}\f$ within the control volume.
      */
    Scalar mu() const
    { return mu_; }

    /*!
     * \brief Returns the rock density \f$\mathrm{[kg / m^3]}\f$ within the control volume .
     */
    Scalar rockDensity() const
    { return rockDensity_; }

    /*!
     * \brief Returns the solid displacement \f$\mathrm{[m]}\f$ in space
     * directions dimIdx within the control volume.
     */
    Scalar displacement(int dimIdx) const
    { return displacement_[dimIdx]; }

    /*!
     * \brief Returns the solid displacement vector \f$\mathrm{[m]}\f$
     *  within the control volume.
     */
    const DimVector &displacement() const
    { return displacement_; }

protected:
    PrimaryVariables primaryVars_;
    DimVector displacement_;
    Scalar lambda_;
    Scalar mu_;
    Scalar rockDensity_;
};

}

#endif
