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
 *
 * \brief Quantities required by the two-phase linear-elastic model which
 *           are defined on a vertex.
 */
#ifndef DUMUX_DECOUPLED_ELASTIC_VOLUME_VARIABLES_HH
#define DUMUX_DECOUPLED_ELASTIC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/2p/implicit/volumevariables.hh>

#include "dumux/geomechanics/el2p/properties.hh"

namespace Dumux {
/*!
 * \ingroup ElTwoPModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase linear-elastic model.
 *
 *        This class inherits from the vertexdata of the two-phase
 *        model and from the vertexdata of the simple
 *        linear-elastic model
 */
template<class TypeTag>
class DecoupledVolumeVariables: public TwoPVolumeVariables<TypeTag> {

    typedef TwoPVolumeVariables<TypeTag> TwoPBase;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {  dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx,
                bool isOldSol)
    {
        primaryVars_ = priVars;

        TwoPBase::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
        primaryVars_ = priVars;

        for (int coordDir = 0; coordDir < dim; ++coordDir)
            displacement_[coordDir] = priVars[Indices::u(coordDir)];

        effFluidDensity_ = this->density(wPhaseIdx) * this->saturation(wPhaseIdx)
                        + this->density(nPhaseIdx) * this->saturation(nPhaseIdx);

        const Dune::FieldVector<Scalar, 3> &lameParams =
                problem.spatialParams().lameParams(element, fvGeometry, scvIdx);

        lambda_ = lameParams[0];
        mu_ = lameParams[1];

        rockDensity_ = problem.spatialParams().rockDensity(element, scvIdx);

        // porosity
        initialPorosity_ = problem.spatialParams().porosity(element,
                                                     fvGeometry,
                                                     scvIdx);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    /*!
     * \brief Return the vector of primary variables
     */
    const Scalar &priVar(int idx) const
    { return primaryVars_[idx]; }

    /*!
     * \brief Sets the evaluation point used in the by the local jacobian.
     */
    void setEvalPoint(const Implementation *ep)
    { }

    /*!
     * \brief Returns the effective effective fluid density within
     *        the control volume.
     */
    Scalar effFluidDensity() const
    { return effFluidDensity_; }

    /*!
       * \brief Returns the Lame parameter lambda within the control volume.
       */
     Scalar lambda() const
     { return lambda_; }

     /*!
       * \brief Returns the Lame parameter mu within the control volume.
       */
     Scalar mu() const
     { return mu_; }

     /*!
      * \brief Returns the rock density within the control volume.
      */
     Scalar rockDensity() const
     { return rockDensity_; }

     /*!
      * \brief Returns the solid displacement in all space
      * directions within the control volume.
      */
     Scalar displacement(int dimIdx) const
     { return displacement_[dimIdx]; }

     /*!
      * \brief Returns the solid displacement vector
      * within the control volume.
      */
     DimVector displacement() const
     { return displacement_; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar initialPorosity() const
    { return initialPorosity_; }

    mutable Scalar divU;
    mutable Scalar effPorosity;
    mutable Scalar volumetricStrain;

protected:
    Scalar effFluidDensity_;
    PrimaryVariables primaryVars_, prevPrimaryVars_;
    DimVector displacement_, prevDisplacement_;
    Scalar lambda_;
    Scalar mu_;
    Scalar rockDensity_;
    Scalar initialPorosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
