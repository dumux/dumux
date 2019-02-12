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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the Stokes box model.
 */
#ifndef DUMUX_STOKES_SECONDARY_VARIABLES_HH
#define DUMUX_STOKES_SECONDARY_VARIABLES_HH

#include "properties.hh"

//from elastic
#include <dumux/discretization/fem/secondaryvariablesbase.hh>

//#include <dumux/implicit/volumevariables.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the Stokes box model.
 */
template <class TypeTag>
class StokesSecondaryVariables : public FemSecondaryVariablesBase<TypeTag>
{
using ParentType = FemSecondaryVariablesBase<TypeTag>;

using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

using Element = typename GridView::template Codim<0>::Entity;
using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

//copied from secondaryvariablesbase
using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        momentumXIdx = Indices::momentumXIdx,
        lastMomentumIdx = Indices::lastMomentumIdx,
        pressureIdx = Indices::pressureIdx
    };

    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    /*!
     * \copydoc SecondaryVariablesBase::update
     */
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const IpData& ipData)
    {
    ParentType::update(elemSol, problem, element, ipData);

    completeFluidState(ParentType::priVars(), problem, element, fluidState_);
    for (int dimIdx=0; dimIdx<dim; ++dimIdx)
        velocity_[dimIdx] = ParentType::priVars()[Indices::momentum(dimIdx)];
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState()
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   FluidState& fluidState
                                   )
    {
        Scalar temperature = problem.temperature();

        fluidState.setTemperature(temperature);
        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        // create NullParameterCache and do dummy update
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        fluidState.setDensity(phaseIdx,
                              FluidSystem::density(fluidState,
                                                   paramCache,
                                                   phaseIdx));

        fluidState.setViscosity(phaseIdx,
                                FluidSystem::viscosity(fluidState,
                                                       paramCache,
                                                       phaseIdx));;
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }
    FluidState &fluidState()
    { return fluidState_; }


    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the fluid within the
     *        sub-control volume.
     */
    Scalar density() const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ of the fluid within the
     *        sub-control volume.
     */
    Scalar molarDensity() const
    { return this->fluidState_.density(phaseIdx) / this->fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the fluid pressure \f$\mathrm{[Pa]}\f$ within
     *        the sub-control volume.
     */
    Scalar pressure() const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature\f$\mathrm{[T]}\f$ inside the sub-control volume.
     */
    Scalar temperature() const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$ \mathrm{[Pa s]} \f$ of the fluid in
     *        the sub-control volume.
     */
    Scalar dynamicViscosity() const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the kinematic viscosity \f$ \frac{m^2}{s} \f$ of the fluid in
     *        the sub-control volume.
     */
    Scalar kinematicViscosity() const
    { return fluidState_.viscosity(phaseIdx) / fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the velocity vector in the sub-control volume.
     */
    const DimVector &velocity() const
    { return velocity_; }



    Scalar artificialCompr() const
    {
        return 10;
    }



    Scalar epsGeneral(DimVector velocity, Scalar pressure) const {
    Scalar eps(0.0),k(0.0), normU(0.0);
    //set specific k, usually 1e-9 < k < 1e-3
    k = 1e-3;

    //calculate norm of velocity U
    for(int i=0; i<dim; i++){
        normU += velocity[i]*velocity[i];
    }
    normU = sqrt(normU);


    if(k*normU/pressure < GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PenEps)){
        return GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PenEps);
    }
    else{
        //absolute value of pressure
        if(pressure < 0){
        pressure *= -1;
        }

        return k*normU/pressure;
     }

    }


    Scalar meshWidthX() const
    {
        Scalar hx(0.0);
        DimVector upperRight(0.0), numbCells(0.0);

        upperRight = GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.UpperRight);
        numbCells = GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.Cells);

        hx = upperRight[0]/numbCells[0];

        return hx;
    }


    Scalar meshWidthY() const
    {
        Scalar hy(0.0);
        DimVector upperRight(0.0), numbCells(0.0);

        upperRight = GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.UpperRight);
        numbCells = GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.Cells);

        hy = upperRight[1]/numbCells[1];

        return hy;
    }



    //viscosity = dynamicViscosity? -> viscosity = kinematicViscosity (Wiki: ReynoldsNumber)
    Scalar stabAlpha(Scalar velocity, Scalar viscosity, Scalar meshWidth) const
    {
        Scalar alpha(0.0);
        Scalar reynolds(0.0);
        Scalar reynoldsHalf(0.0);
        Scalar funcCosH(0.0);

        Scalar H(0.0), alphaDash(0.0);
        reynolds = (velocity*meshWidth) / viscosity;

        H = reynolds/2;

        if(H>=-3 && H<=3){
            alphaDash = H/3;
        }else{
            if(H>0){
                alphaDash = 1;
            }else{
                alphaDash = -1;
            }
        }

        alpha = 0.5*alphaDash;


        return alpha;
    }

    Scalar stabAlpha2(Scalar velocity, Scalar viscosity, Scalar meshWidth) const
    {
        Scalar alpha(0.0);
        Scalar reynolds(0.0);
        Scalar reynoldsHalf(0.0);
        Scalar funcCosH(0.0);

        Scalar H(0.0), alphaDash(0.0);
        reynolds = (velocity*meshWidth) / viscosity;

        H = reynolds/2;

        if(H>=-3 && H<=3){
            alphaDash = H/3;
        }else{
            if(H>0){
                alphaDash = 1;
            }else{
                alphaDash = -1;
            }
        }


        alpha = 0.5*alphaDash;


        return alpha;
    }




    Scalar penaltyEps() const
    {
      return GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PenEps);
    }


    Scalar penaltyEpsTimesP(DimVector velocity) const
    { Scalar k = 1000;
      Scalar normU(0.0);

      for(int i=0; i<dim; i++){
          normU += velocity[i]*velocity[i];
      }
      normU = sqrt(normU);

      Scalar epsP = (k*normU);
      return epsP;
    }


protected:

    DimVector velocity_;
    FluidState fluidState_;

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
