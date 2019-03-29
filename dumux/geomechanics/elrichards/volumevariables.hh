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
 * \brief Quantities required by the richards linear-elastic model which
 *           are defined on a vertex.
 */
#ifndef DUMUX_ELRICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_ELRICHARDS_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/richards/implicit/volumevariables.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include "properties.hh"

namespace Dumux {
/*!
 * \ingroup ElRichardsModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the richards linear-elastic model.
 *
 *        This class inherits from the vertexdata of the richards
 *        model and from the vertexdata of the simple
 *        linear-elastic model
 */
template<class TypeTag>
class ElRichardsVolumeVariables: public RichardsVolumeVariables<TypeTag> {//this name is inherited from the equivalent file in richards folder
//class RichardsVolumeVariables : public ImplicitVolumeVariables<TypeTag>


    typedef RichardsVolumeVariables<TypeTag> RichardsBase;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
 //from richards
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;


     enum{hIdx = Indices::hIdx,
         pwIdx = Indices::pwIdx,
         wPhaseIdx = Indices::wPhaseIdx,
         nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {  dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;//from richards
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
         RichardsBase::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
         primaryVars_ = priVars;

        for (int coordDir = 0; coordDir < dim; ++coordDir)
            displacement_[coordDir] = priVars[Indices::u(coordDir)];
//At the end I should once calculate the fluid density with this method and compare the results with the richards methods which is brought in Density and viscosity part
        //effFluidDensity_ = this->density(wPhaseIdx) * this->saturation(wPhaseIdx)
                        + this->density(nPhaseIdx) * this->saturation(nPhaseIdx);

        const Dune::FieldVector<Scalar, 2> &lameParams =
                problem.spatialParams().lameParams(element, fvGeometry, scvIdx);

        lambda_ = lameParams[0];
        mu_ = lameParams[1];

        rockDensity_ = problem.spatialParams().rockDensity(element, scvIdx);
        Valgrind::CheckDefined(primaryVars_);
        Valgrind::CheckDefined(prevPrimaryVars_);
        Valgrind::CheckDefined(displacement_);
        Valgrind::CheckDefined(prevDisplacement_);
        Valgrind::CheckDefined(lambda_);
        Valgrind::CheckDefined(mu_);
        Valgrind::CheckDefined(rockDensity_);

        fluidState_.setPressure(wPhaseIdx, primaryVars_[pwIdx]);
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
//      */
    //Scalar effFluidDensity() const
    //{ return effFluidDensity_; }


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
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * For the non-wetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the non-wetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

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

    mutable Scalar divU;
    mutable Scalar effPorosity;

protected:
    //Scalar effFluidDensity_;
    PrimaryVariables primaryVars_, prevPrimaryVars_;
    DimVector displacement_, prevDisplacement_;
    Scalar lambda_;
    Scalar mu_;
    Scalar rockDensity_;
    FluidState fluidState_;


private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};


}

#endif
