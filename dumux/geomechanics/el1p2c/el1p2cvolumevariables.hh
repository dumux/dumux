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
 * \brief Quantities required by the single-phase, two-component
 *           linear elasticity model which are defined on a vertex.
 */
#ifndef DUMUX_ELASTIC1P2C_VOLUME_VARIABLES_HH
#define DUMUX_ELASTIC1P2C_VOLUME_VARIABLES_HH


#include <dumux/implicit/1p2c/1p2cvolumevariables.hh>
#include <dumux/implicit/volumevariables.hh>

#include "el1p2cproperties.hh"

namespace Dumux {
/*!
 * \ingroup ElOnePTwoCBoxModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the single-phase, two-component, linear elasticity model.
 *
 *        This class inherits from the volumevariables of the one-phase
 *        two-component model
 */
template<class TypeTag>
class ElOnePTwoCVolumeVariables : public OnePTwoCVolumeVariables<TypeTag>{

    typedef OnePTwoCVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {  dim = GridView::dimension };

    enum {  phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar,dim> DimVector;

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

        ParentType::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
        int vIdxGlobal = problem.vertexMapper().subIndex(element, scvIdx, dim);

        primaryVars_ = priVars;
        prevPrimaryVars_ = problem.model().prevSol()[vIdxGlobal];

        ParentType prev1p2cVolVars;
        prev1p2cVolVars.update(problem.model().prevSol()[vIdxGlobal],
                               problem,
                               element,
                               fvGeometry,
                               scvIdx,
                               true);

        dPressure_ = this->pressure() - prev1p2cVolVars.pressure();

        for (int i = 0; i < dim; ++i){
            displacement_[i] = primaryVars_[Indices::u(i)];
            prevDisplacement_[i] = prevPrimaryVars_[Indices::u(i)];
        }

        dU_ = displacement_ - prevDisplacement_;

        const Dune::FieldVector<Scalar, 2> &lameParams =
                problem.spatialParams().lameParams(element, fvGeometry, scvIdx);

        lambda_ = lameParams[0];
        mu_ = lameParams[1];

        rockDensity_ = problem.spatialParams().rockDensity(element, scvIdx);
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
      * \brief Returns the change in solid displacement \f$\mathrm{[m]}\f$ between
      * the last and the current timestep for the space direction dimIdx within the control volume.
      */
     Scalar dU(int dimIdx) const
     { return dU_[dimIdx]; }

     /*!
      * \brief Returns the change in pressure \f$\mathrm{[Pa]}\f$ between the last and the
      * current timestep within the control volume.
      */
     Scalar dPressure() const
     { return dPressure_; }

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


     /*!
      * \brief the effective porosity and volumetric strain divU is defined as mutable variable here since it
      * is updated within the elementVolumeVariables.
      */
     mutable Scalar effPorosity;
     mutable Scalar divU;

     /*!
     * \brief Returns the mass fraction of a given component in the
     *        given fluid phase within the control volume.
     *
     * \param compIdx The component index
     */
     Scalar massFraction(const int compIdx) const
     { return this->fluidState_.massFraction(phaseIdx, compIdx); }

     /*!
     * \brief Returns the mole fraction of a given component in the
     *        given fluid phase within the control volume.
     *
     * \param compIdx The component index
     */
     Scalar moleFraction(const int compIdx) const
     { return this->fluidState_.moleFraction(phaseIdx, compIdx); }

protected:
    PrimaryVariables primaryVars_, prevPrimaryVars_;
    DimVector displacement_, prevDisplacement_;
    DimVector dU_;
    Scalar dPressure_;
    Scalar lambda_;
    Scalar mu_;
    Scalar rockDensity_;
};

}

#endif
