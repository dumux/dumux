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
#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/grid.hh>

#include <dumux/implicit/model.hh>
#include "properties.hh"

#include <dune/common/fmatrix.hh>

namespace Dumux
{
/*!
 * \ingroup FemModel, instead of BoxStokesModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local Jacobian matrix for problems
 *        using the Stokes box model.
 *
 * This class is also used for the non-isothermal and the two-component Stokes
 * model (static polymorphism).
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);

using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    //added to resemble hookes law
    static constexpr int dimWorld = GridView::dimensionworld;

    //added to resemble elastic
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    enum {
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        lastMomentumIdx = Indices::lastMomentumIdx //!< Index of the last component of the momentum balance
    };
    enum { pressureIdx = Indices::pressureIdx }; //!< Index of the pressure in a solution vector

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;
    using ReferenceElement = typename Dune::ReferenceElement<Scalar, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);


    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);


    static const bool enableUnsymmetrizedVelocityGradient = GET_PROP_VALUE(TypeTag, EnableUnsymmetrizedVelocityGradient);
    static const bool calculateNavierStokes = GET_PROP_VALUE(TypeTag, EnableNavierStokes);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    //copied from hookesLaw
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

 public:
    //added to resemble elastic
    using typename ParentType::FluxTermType;


    /*!
     *  method signature taken from elastic (4 parameters)
     *
     * \brief Evaluate the amount of all conservation quantities
     *        within a finite volume.
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     */
    PrimaryVariables computeStorage(const Element& element,
    const IpData& ipData,
    const SecondaryVariables& secVars,
    const ElementSolution& elemSol) const
    {
    PrimaryVariables storage(0.0);

    for (int dir = 0; dir < dim; ++dir)
        storage[Indices::momentum(dir)] += secVars.velocity()[dir];

        return storage; //changed from void to PrimaryVariables as its required in geomechanics
    }

    /*!
     *  method signature taken from elastic (4 parameters)
     *
     * \brief Evaluate the stresses.
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     */
    FluxTermType computeFlux(const Element& element,
    const IpData& ipData,
    const SecondaryVariables& secVars,
    const ElementSolution& elemSol) const
    {
    FluxTermType flux(0.0);

    StressTensor gradV(0.0);
    for (int dir = 0; dir < dim; ++dir){
        for (unsigned int i = 0; i < elemSol.size(); ++i){
            gradV[dir].axpy(elemSol[i][Indices::momentum(dir)], ipData.shapeGradients(i));
        }
    }

    StressTensor gradSym(0.0);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dimWorld; j++){
            gradSym[i][j] = gradV[i][j] + gradV[j][i];
        }
    }

    StressTensor sigma(0.0);
    sigma = gradSym;
    sigma *= secVars.kinematicViscosity();

//    sigma *= -1;

    Scalar divV(0.0);
    for (unsigned int i = 0; i < elemSol.size(); ++i){
        for (int dir = 0; dir < dim; ++dir){
            divV += ipData.shapeGradients(i)[dir]*elemSol[i][Indices::momentum(dir)];
        }
    }

//    Scalar eps = secVars.penaltyEps();
    Scalar eps = secVars.epsGeneral(secVars.velocity(), secVars.pressure());
//    eps /= secVars.density();

    divV /= eps;
    divV /= secVars.density();

//    divV = -secVars.pressure();
//    divV /= secVars.density();
    /*******************
     ***************
     */
//mit dichte multiplizieren?
   // divV *= secVars.density();



    for(int i=0; i<dim; i++){
        sigma[i][i] += divV;
    }


    StressTensor advTerm(0.0);

    advTerm[0][0] = secVars.velocity()[0]*secVars.velocity()[0];
    advTerm[0][1] = secVars.velocity()[0]*secVars.velocity()[1];
    advTerm[1][0] = secVars.velocity()[1]*secVars.velocity()[0];
    advTerm[1][1] = secVars.velocity()[1]*secVars.velocity()[1];

    for(int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; ++momentumIdx)
    {
        for(int col=0; col<dim; col++)
        {
            flux[momentumIdx][col] -= sigma[momentumIdx][col];
            flux[momentumIdx][col] += advTerm[momentumIdx][col];
        }
    }

   auto velocityU = secVars.velocity();

   // velocityU *= secVars.density();

    //add massBalance to flux
    for(int col=0; col<dim; col++)
    {
         flux[dim][col]+= velocityU[col];
    }

    return flux;
    }




    /*!
     *  method signature taken from elastic (4 parameters)
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     *
     */
    PrimaryVariables computeSource(const Element& element,
    const IpData& ipData,
    const SecondaryVariables& secVars,
    const ElementSolution& elemSol)
    {
      PrimaryVariables source(0.0);

      Scalar epsTimesP = secVars.epsGeneral(secVars.velocity(), secVars.pressure());
      epsTimesP*= secVars.pressure();
   //   epsTimesP*= secVars.kinematicViscosity();
      source[dim] -= epsTimesP;

      return source;
    }


    Scalar massStab(const Element& element,
            const IpData& ipData,
            const SecondaryVariables& secVars,
            const ElementSolution& elemSol)
    {
        //for mass balance
        Scalar massStab(0.0);
        for (unsigned int i = 0; i < elemSol.size(); ++i){
            massStab += elemSol[i][0]*ipData.shapeGradients(i)[0] + elemSol[i][1]*ipData.shapeGradients(i)[1];
        }

        return massStab;
    }




    DimVector computeStabilizationTerms(const Element& element,
            const IpData& ipData,
            const SecondaryVariables& secVars,
            const ElementSolution& elemSol)
    {
        DimVector stabTerms(0.0);

        StressTensor advTerm(0.0);
        DimVector divAdvTerm(0.0);

        for (unsigned int i = 0; i < elemSol.size(); ++i){
            divAdvTerm[0] += 2*elemSol[i][0]*elemSol[i][0]*ipData.shapeValues(i)*ipData.shapeGradients(i)[0] + 2*elemSol[i][1]*elemSol[i][0]*ipData.shapeValues(i)*ipData.shapeGradients(i)[1];
            divAdvTerm[1] += 2*elemSol[i][0]*elemSol[i][1]*ipData.shapeValues(i)*ipData.shapeGradients(i)[0] + 2*elemSol[i][1]*elemSol[i][1]*ipData.shapeValues(i)*ipData.shapeGradients(i)[1];
        }

        //for diffusive part, doesnt fall off entirely
        Scalar meshQuot = secVars.meshWidthX()*secVars.meshWidthY();
        meshQuot = 1/meshQuot;

        DimVector divGradVel(0.0);
        for (unsigned int i = 0; i < elemSol.size(); ++i){

            divGradVel[0] += elemSol[i][1]*meshQuot;
            divGradVel[1] += elemSol[i][0]*meshQuot;
        }


        Scalar tau(0.0), S(0.0);

        tau = 1*secVars.stabAlpha(secVars.velocity()[0], secVars.kinematicViscosity(), 1)*secVars.velocity()[0] +
              1*secVars.stabAlpha(secVars.velocity()[1], secVars.kinematicViscosity(), 1)*secVars.velocity()[1];
        tau /= 2;

        S = secVars.velocity()[0]*secVars.velocity()[0] + secVars.velocity()[1]*secVars.velocity()[1];

        if(S != 0){
        tau /= S;

        divAdvTerm *= tau;

        stabTerms[0] = divAdvTerm[0];
        stabTerms[1] = divAdvTerm[1];

        }
        return stabTerms;
    }


};

}

#endif
