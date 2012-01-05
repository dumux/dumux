// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2007-2009 by Jochen Fritz                                 *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef DUMUX_FVPRESSURE2P2C_MULTIPHYSICS_HH
#define DUMUX_FVPRESSURE2P2C_MULTIPHYSICS_HH

// dumux environment
#include <dumux/decoupled/2p2c/fvpressure2p2c.hh>
#include <dumux/decoupled/2p2c/pseudo1p2cfluidstate.hh>

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{
//! The finite volume model for the solution of the compositional pressure equation
/*! \ingroup multiphysics
 *  Provides a Finite Volume implementation for the pressure equation of a gas-liquid
 *  system with two components. An IMPES-like method is used for the sequential
 *  solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 *  Isothermal conditions and local thermodynamic
 *  equilibrium are assumed.  Gravity is included.
 *  \f[
         c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right)
          = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 *  \f]
 *  where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 *  \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 *
 *  The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * The model domain is automatically divided
 * in a single-phase and a two-phase domain. The full 2p2c model is only evaluated within the
 * two-phase subdomain, whereas a single-phase transport model is computed in the rest of the
 * domain.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVPressure2P2CMultiPhysics : public FVPressure2P2C<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<Scalar, 2> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;


//! Access functions to the current problem object
    Problem& problem()
    {    return this->problem_;   }
    const Problem& problem() const
    {    return this->problem_;   }
public:
    //function which assembles the system of equations to be solved
    void assemble(bool first);

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    //! \brief Write data files
     /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem().variables().addOutputVtkFields(writer);

        int size_ = problem().variables().gridSize();
        // add multiphysics stuff
        Dune::BlockVector<Dune::FieldVector<int,1> > *subdomainPtr = writer.template allocateManagedBuffer<int, 1> (size_);
        *subdomainPtr = problem().variables().subdomain();
        writer.attachCellData(*subdomainPtr, "subdomain");

#if DUNE_MINIMAL_DEBUG_LEVEL <= 2
        // add debug stuff
        typename SolutionTypes::ScalarSolution *numdensityW = writer.allocateManagedBuffer (size_);
        typename SolutionTypes::ScalarSolution *numdensityNW = writer.allocateManagedBuffer (size_);
        *numdensityW = problem().variables().numericalDensity(wPhaseIdx);
        *numdensityNW = problem().variables().numericalDensity(nPhaseIdx);
        writer.attachCellData(*numdensityW, "numerical density (mass/volume) w_phase");
        writer.attachCellData(*numdensityNW, "numerical density (mass/volume) nw_phase");

        typename SolutionTypes::ScalarSolution *errorCorrPtr = writer.allocateManagedBuffer (size_);
        *errorCorrPtr = problem().variables().errorCorrection();
//        int size = subdomainPtr.size();
        writer.attachCellData(*errorCorrPtr, "Error Correction");

        typename SolutionTypes::ScalarSolution *dv_dpPtr = writer.allocateManagedBuffer (size_);
        *dv_dpPtr = problem().variables().dv_dp();
        writer.attachCellData(*dv_dpPtr, "dv_dp");

        typename SolutionTypes::ScalarSolution *dV_dC1Ptr = writer.allocateManagedBuffer (size_);
        typename SolutionTypes::ScalarSolution *dV_dC2Ptr = writer.allocateManagedBuffer (size_);
        *dV_dC1Ptr = problem().variables().dv()[0];
        *dV_dC2Ptr = problem().variables().dv()[1];
        writer.attachCellData(*dV_dC1Ptr, "dV_dC1");
        writer.attachCellData(*dV_dC2Ptr, "dV_dC2");

        Dune::BlockVector<Dune::FieldVector<double,1> > *updEstimate1 = writer.allocateManagedBuffer (size_);
        Dune::BlockVector<Dune::FieldVector<double,1> > *updEstimate2 = writer.allocateManagedBuffer (size_);
        *updEstimate1 = problem().variables().updateEstimate()[0];
        *updEstimate2 = problem().variables().updateEstimate()[1];
        writer.attachCellData(*updEstimate1, "updEstimate comp 1");
        writer.attachCellData(*updEstimate2, "updEstimate comp 2");
#endif

        return;
    }

    //! Constructs a FVPressure2P2CPC object
    /**
     * \param problem a problem class object
     */
    FVPressure2P2CMultiPhysics(Problem& problem) : FVPressure2P2C<TypeTag>(problem),
            gravity(problem.gravity()), timer_(false)
    {
        // resize multiphysics part
        problem.variables().subdomain().resize(problem.variables().gridSize());
        nextSubdomain.resize(problem.variables().gridSize());

        // assign whole domain to most complex subdomain, => 2p
        problem.variables().subdomain() = 2;
    }

private:
    // subdomain map
    Dune::BlockVector<Dune::FieldVector<int,1> > nextSubdomain;  //! vector holding next subdomain

protected:
    const GlobalPosition& gravity; //!< vector including the gravity constant
    static constexpr int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    Dune::Timer timer_;
};

//! function which assembles the system of equations to be solved
/** for first == true, this function assembles the matrix and right hand side for
 * the solution of the pressure field in the same way as in the class FVPressure2P.
 * for first == false, the approach is changed to \f[-\frac{\partial V}{\partial p}
 * \frac{\partial p}{\partial t}+\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}\nabla\cdot
 * \left(\sum_{\alpha}C_{\alpha}^{\kappa}\mathbf{v}_{\alpha}\right)
 * =\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}q^{\kappa} \f]. See Paper SPE 99619.
 * This is done to account for the volume effects which appear when gas and liquid are dissolved in each other.
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    // initialization: set the fluid volume derivatives to zero
    Scalar dv_dC1(0.), dv_dC2(0.);
    for (int i = 0; i < problem().variables().gridSize(); i++)
    {
        problem().variables().dv(i, wPhaseIdx) = 0;     // dv_dC1 = dV/dm1
        problem().variables().dv(i, nPhaseIdx) = 0;     // dv / dC2
        problem().variables().dv_dp(i) = 0;      // dv / dp
    }

    // determine maximum error to scale error-term
    Scalar timestep_ = problem().timeManager().timeStepSize();
    Scalar maxErr = fabs(problem().variables().volErr().infinity_norm()) / timestep_;

    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // acess geometry
        const typename ElementIterator::Entity::Geometry &geo = eIt->geometry();

        // get global coordinate of cell center
        const GlobalPosition& globalPos = geo.center();

        // cell index
        int globalIdxI = problem().variables().index(*eIt);

        // cell volume & perimeter, assume linear map here
        Scalar volume = eIt->geometry().volume();
        Scalar perimeter = problem().variables().perimeter(globalIdxI);

        Scalar densityWI = problem().variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem().variables().densityNonwetting(globalIdxI);

        /****************implement sources and sinks************************/
        PrimaryVariables source(NAN);
        problem().source(source, *eIt);
        if(first || problem().variables().subdomain(globalIdxI) == 1)
        {
                source[Indices::contiWEqIdx] /= densityWI;
                source[Indices::contiNEqIdx] /= densityNWI;
        }
        else
        {
                // derivatives of the fluid volume with respect to concentration of components, or pressure
                if (problem().variables().dv(globalIdxI, wPhaseIdx) == 0)
                    this->volumeDerivatives(globalPos, *eIt,
                            problem().variables().dv(globalIdxI, wPhaseIdx),
                            problem().variables().dv(globalIdxI, nPhaseIdx),
                            problem().variables().dv_dp(globalIdxI));

                source[Indices::contiWEqIdx] *= problem().variables().dv(globalIdxI, wPhaseIdx);        // note: dV_[i][1] = dv_dC1 = dV/dm1
                source[Indices::contiNEqIdx] *= problem().variables().dv(globalIdxI, nPhaseIdx);
        }
        this->f_[globalIdxI] = volume * (source[Indices::contiWEqIdx] + source[Indices::contiNEqIdx]);
        /********************************************************************/

        // get absolute permeability
        FieldMatrix permeabilityI(problem().spatialParameters().intrinsicPermeability(globalPos, *eIt));

        // get mobilities and fractional flow factors
        Scalar lambdaWI = problem().variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem().variables().mobilityNonwetting(globalIdxI);
        Scalar fractionalWI=0, fractionalNWI=0;
        if (first)
        {
            fractionalWI = lambdaWI / (lambdaWI+ lambdaNWI);
            fractionalNWI = lambdaNWI / (lambdaWI+ lambdaNWI);
        }

        Scalar pcI = problem().variables().capillaryPressure(globalIdxI);

        // prepare meatrix entries
        Scalar entry=0.;
        Scalar rightEntry=0.;

        // iterate over all faces of the cell
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // get normal vector
            const GlobalPosition& unitOuterNormal = isIt->centerUnitOuterNormal();

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                // access neighbor
                const ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem().variables().index(*neighborPointer);

                // gemotry info of neighbor
                const typename ElementIterator::Entity::Geometry &geoNeighbor = neighborPointer->geometry();
                const GlobalPosition& globalPosNeighbor = geoNeighbor.center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                FieldMatrix permeabilityJ
                    = problem().spatialParameters().intrinsicPermeability(globalPosNeighbor,
                                                                            *neighborPointer);

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);
                Dumux::harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

                Dune::FieldVector<Scalar, dim> permeability(0);
                meanPermeability.mv(unitDistVec, permeability);

                // get mobilities in neighbor
                Scalar lambdaWJ = problem().variables().mobilityWetting(globalIdxJ);
                Scalar lambdaNWJ = problem().variables().mobilityNonwetting(globalIdxJ);

                // phase densities in cell in neighbor
                Scalar densityWJ = problem().variables().densityWetting(globalIdxJ);
                Scalar densityNWJ = problem().variables().densityNonwetting(globalIdxJ);

                Scalar pcJ = problem().variables().capillaryPressure(globalIdxJ);

                Scalar rhoMeanW = 0.5 * (densityWI + densityWJ);
                Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWJ);

                // reset potential gradients, density
                Scalar potentialW = 0;
                Scalar potentialNW = 0;
                Scalar densityW = 0;
                Scalar densityNW = 0;

                if (first)     // if we are at the very first iteration we can't calculate phase potentials
                {
                    // get fractional flow factors in neigbor
                    Scalar fractionalWJ = lambdaWJ / (lambdaWJ+ lambdaNWJ);
                    Scalar fractionalNWJ = lambdaNWJ / (lambdaWJ+ lambdaNWJ);

                    // perform central weighting
                    Scalar lambda = (lambdaWI + lambdaWJ) * 0.5 + (lambdaNWI + lambdaNWJ) * 0.5;

                    entry = fabs(lambda*faceArea*(permeability*unitOuterNormal)/(dist));

                    Scalar factor = (fractionalWI + fractionalWJ) * (rhoMeanW) * 0.5 + (fractionalNWI + fractionalNWJ) * (rhoMeanNW) * 0.5;
                    rightEntry = factor * lambda * faceArea * (permeability * gravity) * (unitOuterNormal * unitDistVec);
                }
                else if (problem().variables().subdomain(globalIdxI) == 1) //the current cell in the 1p domain
                {
                    // 1p => no pC => only 1 pressure, potential
                    Scalar potential = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI]
                                 - problem().variables().pressure()[globalIdxJ]) / (dist*dist);
                    potentialW = potentialNW = potential;

                    potentialW += densityW * (unitDistVec * gravity);
                    potentialNW += densityNW * (unitDistVec * gravity);

                    double lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                        {
                        lambdaW = lambdaWI;
                        densityW = densityWI;
                        }
                    else
                        {
                        lambdaW = lambdaWJ;
                        densityW = densityWJ;
                        }

                    if (potentialNW >= 0.)
                        {
                        lambdaNW = lambdaNWI;
                        densityNW = densityNWI;
                        }
                    else
                        {
                        lambdaNW = lambdaNWJ;
                        densityNW = densityNWJ;
                        }

                    entry  = (lambdaW + lambdaNW) * faceArea * (permeability * unitOuterNormal) / (dist);
                    rightEntry = (densityW * lambdaW + densityNW * lambdaNW) * faceArea * (permeability * gravity) * (unitOuterNormal * unitDistVec);
                }
                else    //current cell is in the 2p domain
                {
                    Scalar graddv_dC1(0.), graddv_dC2(0.);

                    //check neighboring cell
                    if (problem().variables().subdomain(globalIdxJ)== 2)
                    {
                        // determine volume derivatives
                        if (problem().variables().dv(globalIdxJ, wPhaseIdx) == 0)
                        this->volumeDerivatives(globalPosNeighbor, *neighborPointer,
                                problem().variables().dv(globalIdxJ, wPhaseIdx),
                                problem().variables().dv(globalIdxJ, nPhaseIdx),
                                problem().variables().dv_dp(globalIdxJ));
                    dv_dC1 = (problem().variables().dv(globalIdxJ, wPhaseIdx)
                                + problem().variables().dv(globalIdxI, wPhaseIdx)) / 2; // dV/dm1= dV/dC^1
                    dv_dC2 = (problem().variables().dv(globalIdxJ, nPhaseIdx)
                                + problem().variables().dv(globalIdxI, nPhaseIdx)) / 2;

                    graddv_dC1 = (problem().variables().dv(globalIdxJ, wPhaseIdx)
                                            + problem().variables().dv(globalIdxI, wPhaseIdx)) / dist;
                    graddv_dC2 = (problem().variables().dv(globalIdxJ, nPhaseIdx)
                                            + problem().variables().dv(globalIdxI, nPhaseIdx)) / dist;
                    }
                    else
                    {
                        dv_dC1 = problem().variables().dv(globalIdxI, wPhaseIdx); //only regard volume changes in 2p area => no averageing
                        dv_dC2 = problem().variables().dv(globalIdxI, nPhaseIdx);

                    }

//                    potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                    potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
//
//                    densityW = (potentialW > 0.) ? densityWI : densityWJ;
//                    densityNW = (potentialNW > 0.) ? densityNWI : densityNWJ;
//
//                    densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                    densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                    //jochen: central weighting for gravity term
                    densityW = rhoMeanW; densityNW = rhoMeanNW;

                    switch (pressureType)    //Markus: hab (unitOuterNormal * distVec)/dist hinzugefuegt
                    {
                    case pw:
                    {
                        potentialW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI]
                                - problem().variables().pressure()[globalIdxJ]) / (dist*dist);
                        potentialNW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI]
                                - problem().variables().pressure()[globalIdxJ] + pcI - pcJ) / (dist*dist);
                        break;
                    }
                    case pn:
                    {
                        potentialW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI]
                                - problem().variables().pressure()[globalIdxJ] - pcI + pcJ) / (dist*dist);
                        potentialNW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI]
                                - problem().variables().pressure()[globalIdxJ]) / (dist*dist);
                        break;
                    }
                    }
                    potentialW += densityW * (unitDistVec * gravity);
                    potentialNW += densityNW * (unitDistVec * gravity);

                    // initialize convenience shortcuts
                    Scalar lambdaW, lambdaN;
                    Scalar dV_w(0.), dV_n(0.);        // dV_a = \sum_k \rho_a * dv/dC^k * X^k_a
                    Scalar gV_w(0.), gV_n(0.);        // multipaper eq(3.3) zeile 3 analogon dV_w


                    //do the upwinding of the mobility & density depending on the phase potentials
                    if (potentialW >= 0.)
                    {
                        dV_w = (dv_dC1 * problem().variables().wet_X1(globalIdxI) + dv_dC2 * (1. - problem().variables().wet_X1(globalIdxI)));
                        dV_w *= densityWI;
                        lambdaW = problem().variables().mobilityWetting(globalIdxI);

                        if (problem().variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_w = (graddv_dC1 * problem().variables().wet_X1(globalIdxI) + graddv_dC2 * (1. - problem().variables().wet_X1(globalIdxI)));
                            gV_w *= densityWI;
                        }
                    }
                    else
                    {
                        dV_w = (dv_dC1 * problem().variables().wet_X1(globalIdxJ) + dv_dC2 * (1. - problem().variables().wet_X1(globalIdxJ)));
                        dV_w *= densityWJ;
                        lambdaW = problem().variables().mobilityWetting(globalIdxJ);
                        if (problem().variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_w = (graddv_dC1 * problem().variables().wet_X1(globalIdxJ) + graddv_dC2 * (1. - problem().variables().wet_X1(globalIdxJ)));
                            gV_w *= densityWJ;
                        }
                    }
                    if (potentialNW >= 0.)
                    {
                        dV_n = (dv_dC1 * problem().variables().nonwet_X1(globalIdxI) + dv_dC2 * (1. - problem().variables().nonwet_X1(globalIdxI)));
                        dV_n *= densityNWI;
                        lambdaN = problem().variables().mobilityNonwetting(globalIdxI);
                        if (problem().variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_n = (graddv_dC1 * problem().variables().nonwet_X1(globalIdxI) + graddv_dC2 * (1. - problem().variables().nonwet_X1(globalIdxI)));
                            gV_n *= densityNWI;
                        }
                    }
                    else
                    {
                        dV_n = (dv_dC1 * problem().variables().nonwet_X1(globalIdxJ) + dv_dC2 * (1. - problem().variables().nonwet_X1(globalIdxJ)));
                        dV_n *= densityNWJ;
                        lambdaN = problem().variables().mobilityNonwetting(globalIdxJ);
                        if (problem().variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_n = (graddv_dC1 * problem().variables().nonwet_X1(globalIdxJ) + graddv_dC2 * (1. - problem().variables().nonwet_X1(globalIdxJ)));
                            gV_n *= densityNWJ;
                        }
                    }
                    //calculate current matrix entry
                    entry = faceArea * (lambdaW * dV_w + lambdaN * dV_n);

                    //calculate right hand side
                    rightEntry = faceArea  * (unitOuterNormal * unitDistVec) * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);

                    if (problem().variables().subdomain(globalIdxJ)!= 1) // complex 2p subdomain
                    {
                        //subtract area integral
                        entry -= volume * faceArea / perimeter * (lambdaW * gV_w + lambdaN * gV_n);
                        rightEntry -= volume * faceArea / perimeter * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n);
                    }
                    entry *= fabs((permeability*unitOuterNormal)/(dist));
                    rightEntry *= (permeability * gravity);

                    // include capillary pressure fluxes
                    switch (pressureType)
                    {
                    case pw:
                        {
                            // calculate capillary pressure gradient
                            Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                            pCGradient *= (pcI - pcJ) / dist;

                            //add capillary pressure term to right hand side
                            rightEntry += lambdaN * dV_n * (permeability * pCGradient) * faceArea
                                         - lambdaN * gV_n * (permeability * pCGradient) * volume * faceArea / perimeter;
                            break;
                        }
                    case pn:
                        {
                            // calculate capillary pressure gradient
                            Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                            pCGradient *= (pcI - pcJ) / dist;

                            //add capillary pressure term to right hand side
                            rightEntry -= lambdaW * dV_w * (permeability * pCGradient) * faceArea
                                            - lambdaW * gV_w * (permeability * pCGradient) * volume * faceArea / perimeter;
                            break;
                        }
                    }
                }   // end !first

                //set right hand side
                this->f_[globalIdxI] -= rightEntry;
                // set diagonal entry
                this->A_[globalIdxI][globalIdxI] += entry;
                // set off-diagonal entry
                this->A_[globalIdxI][globalIdxJ] = -entry;
            }   // end neighbor
            /****************************************************/


            /************* boundary face ************************/
            else
            {
                // get volume derivatives inside the cell
                dv_dC1 = problem().variables().dv(globalIdxI, wCompIdx);
                dv_dC2 = problem().variables().dv(globalIdxI, nCompIdx);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                // geometrical information
                GlobalPosition distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary condition for boundary face center
                typename GET_PROP_TYPE(TypeTag, BoundaryTypes) bcType;
                problem().boundaryTypes(bcType, *isIt);

                // prepare pressure boundary condition
                PhaseVector pressBC(0.);
                Scalar pcBound (0.);

                /**********         Dirichlet Boundary        *************/
                if (bcType.isDirichlet(Indices::pressureEqIdx))
                {
                    //permeability vector at boundary
                    Dune::FieldVector<Scalar, dim> permeability(0);
                    permeabilityI.mv(unitDistVec, permeability);

                    // create a fluid state for the boundary
                    FluidState BCfluidState;

                    //read boundary values
                    PrimaryVariables primaryVariablesOnBoundary(NAN);
                    problem().dirichlet(primaryVariablesOnBoundary, *isIt);

                    if (first)
                    {
                        Scalar lambda = lambdaWI+lambdaNWI;
                        this->A_[globalIdxI][globalIdxI] += lambda * faceArea * (permeability * unitOuterNormal) / (dist);
                        Scalar pressBC = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                        this->f_[globalIdxI] += lambda * faceArea * pressBC * (permeability * unitOuterNormal) / (dist);
                        rightEntry = (fractionalWI * densityWI
                                             + fractionalNWI * densityNWI)
                                             * lambda * faceArea * (unitOuterNormal * unitDistVec) *(permeability*gravity);
                        this->f_[globalIdxI] -= rightEntry;
                    }
                    else    //not first
                    {
                        // read boundary values
                        problem().transportModel().evalBoundary(globalPosFace,
                                                                    isIt,
                                                                    eIt,
                                                                    BCfluidState,
                                                                    pressBC);
                        pcBound = pressBC[nPhaseIdx] - pressBC[wPhaseIdx];

                        // determine fluid properties at the boundary
                        Scalar lambdaWBound = 0.;
                        Scalar lambdaNWBound = 0.;

                        Scalar densityWBound =
                            FluidSystem::density(BCfluidState, wPhaseIdx);
                        Scalar densityNWBound =
                            FluidSystem::density(BCfluidState, nPhaseIdx);
                        Scalar viscosityWBound =
                            FluidSystem::viscosity(BCfluidState, wPhaseIdx);
                        Scalar viscosityNWBound =
                            FluidSystem::viscosity(BCfluidState, nPhaseIdx);

                        // mobility at the boundary
                        switch (GET_PROP_VALUE(TypeTag, BoundaryMobility))
                        {
                        case Indices::satDependent:
                            {
                            lambdaWBound = BCfluidState.saturation(wPhaseIdx)
                                    / viscosityWBound;
                            lambdaNWBound = BCfluidState.saturation(nPhaseIdx)
                                    / viscosityNWBound;
                            break;
                            }
                        case Indices::permDependent:
                            {
                            lambdaWBound = MaterialLaw::krw(
                                    problem().spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityWBound;
                            lambdaNWBound = MaterialLaw::krn(
                                    problem().spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityNWBound;
                            break;
                            }
                        }

                        Scalar rhoMeanW = 0.5 * (densityWI + densityWBound);
                        Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWBound);

                        Scalar potentialW = 0, potentialNW = 0;

//                            potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                            potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
//
//                            // do potential upwinding according to last potGradient vs Jochen: central weighting
//                            densityW = (potentialW > 0.) ? densityWI : densityWBound;
//                            densityNW = (potentialNW > 0.) ? densityNWI : densityNWBound;
//
//                            densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                            densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                        Scalar densityW=rhoMeanW;
                        Scalar densityNW=rhoMeanNW;

                        //calculate potential gradient
                        switch (pressureType)
                        {
                        case pw:
                            {
                                potentialW = (unitOuterNormal * distVec) *
                                        (problem().variables().pressure()[globalIdxI] - pressBC[wPhaseIdx]) / (dist * dist);
                                potentialNW = (unitOuterNormal * distVec) *
                                        (problem().variables().pressure()[globalIdxI] + pcI
                                                - pressBC[wPhaseIdx] - pcBound) / (dist * dist);
                                break;
                            }
                        case pn:
                            {
                                potentialW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI] - pcI -
                                        pressBC[nPhaseIdx] + pcBound)
                                        / (dist * dist);
                                potentialNW = (unitOuterNormal * distVec) * (problem().variables().pressure()[globalIdxI] -
                                        pressBC[nPhaseIdx]) / (dist * dist);
                                break;
                            }
                        }

                        potentialW += densityW * (unitDistVec * gravity);
                        potentialNW += densityNW * (unitDistVec * gravity);

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW, lambdaNW;
                        Scalar dV_w, dV_n;     // no gV_a necessary, because dV/dc is assumed spatially constant at the boarder cell

                        if(problem().variables().subdomain(globalIdxI)==1)    // easy 1p subdomain
                        {
                            if (potentialW >= 0.)
                            {
                                densityW = (potentialW == 0) ? rhoMeanW : densityWI;
                                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaWI;
                            }
                            else
                            {
                                densityW = densityWBound;
                                lambdaW = lambdaWBound;
                            }
                            if (potentialNW >= 0.)
                            {
                                densityNW = (potentialNW == 0) ? rhoMeanNW : densityNWI;
                                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNWI;
                            }
                            else
                            {
                                densityNW = densityNWBound;
                                lambdaNW = lambdaNWBound;
                            }

                            //calculate current matrix entry
                            entry = (lambdaW + lambdaNW) * ((permeability * unitDistVec) / dist) * faceArea
                                    * (unitOuterNormal * unitDistVec);

                            //calculate right hand side
                            rightEntry = (lambdaW * densityW + lambdaNW * densityNW) * (permeability * gravity)
                                    * faceArea ;
                        }
                        else    // 2p subdomain
                        {
                            if (potentialW >= 0.)
                            {
                                densityW = (potentialW == 0) ? rhoMeanW : densityWI;
                                dV_w = (dv_dC1 * problem().variables().wet_X1(globalIdxI)
                                         + dv_dC2 * (1. - problem().variables().wet_X1(globalIdxI)));
                                dV_w *= densityW;
                                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaWI;
                            }
                            else
                            {
                                densityW = densityWBound;
                                dV_w = (dv_dC1 * BCfluidState.massFraction(wPhaseIdx, wCompIdx) + dv_dC2 * BCfluidState.massFraction(wPhaseIdx, nCompIdx));
                                dV_w *= densityW;
                                lambdaW = lambdaWBound;
                            }
                            if (potentialNW >= 0.)
                            {
                                densityNW = (potentialNW == 0) ? rhoMeanNW : densityNWI;
                                dV_n = (dv_dC1 * problem().variables().nonwet_X1(globalIdxI)
                                        + dv_dC2 * (1. - problem().variables().nonwet_X1(globalIdxI)));
                                dV_n *= densityNW;
                                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNWI;
                            }
                            else
                            {
                                densityNW = densityNWBound;
                                dV_n = (dv_dC1 * BCfluidState.massFraction(nPhaseIdx, wCompIdx) + dv_dC2 * BCfluidState.massFraction(nPhaseIdx, nCompIdx));
                                dV_n *= densityNW;
                                lambdaNW = lambdaNWBound;
                            }

                            //calculate current matrix entry
                            entry = (lambdaW * dV_w + lambdaNW * dV_n) * ((permeability * unitDistVec) / dist) * faceArea
                                    * (unitOuterNormal * unitDistVec);

                            //calculate right hand side
                            rightEntry = (lambdaW * densityW * dV_w + lambdaNW * densityNW * dV_n) * (permeability * gravity)
                                    * faceArea ;

                            // include capillary pressure fluxes
                            switch (pressureType)
                            {
                            case pw:
                                {
                                    // calculate capillary pressure gradient
                                    Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                                    pCGradient *= (pcI - pcBound) / dist;

                                    //add capillary pressure term to right hand side
                                    rightEntry += lambdaNW * dV_n * (permeability * pCGradient) * faceArea;
                                    break;
                                }
                            case pn:
                                {
                                    // calculate capillary pressure gradient
                                    Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                                    pCGradient *= (pcI - pcBound) / dist;

                                    //add capillary pressure term to right hand side
                                    rightEntry -= lambdaW * dV_w * (permeability * pCGradient) * faceArea;
                                    break;
                                }
                            }
                        } //end 2p subdomain


                        // set diagonal entry and right hand side entry
                        this->A_[globalIdxI][globalIdxI] += entry;
                        this->f_[globalIdxI] += entry * primaryVariablesOnBoundary[Indices::pressureEqIdx];
                        this->f_[globalIdxI] -= rightEntry * (unitOuterNormal * unitDistVec);
                    }    //end of if(first) ... else{...
                }   // end dirichlet

                /**********************************
                 * set neumann boundary condition
                 **********************************/
                else
                {
                    PrimaryVariables J(NAN);
                    problem().neumann(J, *isIt);
                    if (first || problem().variables().subdomain(globalIdxI)==1)
                    {
                        J[contiWEqIdx] /= densityWI;
                        J[contiNEqIdx] /= densityNWI;
                    }
                    else
                    {
                        J[contiWEqIdx] *= dv_dC1;
                        J[contiNEqIdx] *= dv_dC2;
                    }

                    this->f_[globalIdxI] -= (J[contiWEqIdx] + J[contiNEqIdx]) * faceArea;
                }
                /*************************************************/
            }
        } // end all intersections

        // compressibility term
        if (!first && timestep_ != 0.)
        {
            if (problem().variables().dv_dp(globalIdxI) == 0.)
            {
                // if incompressible fluids are used, the following has to be changed
                assert(problem().variables().subdomain(globalIdxI)==1);

                Scalar incp = 1e-2;
                // numerical derivative of fluid volume with respect to pressure
                Scalar p_ = problem().variables().pressure()[globalIdxI] + incp;
                Scalar Z1 = problem().variables().totalConcentration(globalIdxI, wCompIdx)
                        / (problem().variables().totalConcentration(globalIdxI, wCompIdx)
                                + problem().variables().totalConcentration(globalIdxI, nCompIdx));
                PseudoOnePTwoCFluidState<TypeTag> pseudoFluidState;
                pseudoFluidState.update(Z1, p_, problem().variables().saturation(globalIdxI), problem().temperatureAtPos(globalPos));
                if(problem().variables().saturation(globalIdxI)==1)  //only w-phase
                {
                    Scalar v_w_ = 1. / FluidSystem::density(pseudoFluidState, wPhaseIdx);
                    problem().variables().dv_dp(globalIdxI) = (problem().variables().totalConcentration(globalIdxI, wCompIdx)
                            + problem().variables().totalConcentration(globalIdxI, nCompIdx))
                            * ( v_w_  - 1./densityWI) / incp;

                    if (problem().variables().dv_dp(globalIdxI) > 0) // this is not physically possible!!
                    {
                        p_ -= 2*incp;
                        Scalar v_w_ = 1. / FluidSystem::density(pseudoFluidState, wPhaseIdx);
                        problem().variables().dv_dp(globalIdxI) = (problem().variables().totalConcentration(globalIdxI, wCompIdx)
                                                    + problem().variables().totalConcentration(globalIdxI, nCompIdx))
                                                    * ( v_w_  - 1./densityWI) / -incp;
                    }
                }
                if(problem().variables().saturation(globalIdxI)==0.)  //only nw-phase
                {
                    Scalar v_n_ = 1. / FluidSystem::density(pseudoFluidState, nPhaseIdx);
                    problem().variables().dv_dp(globalIdxI) = (problem().variables().totalConcentration(globalIdxI, wCompIdx)
                            + problem().variables().totalConcentration(globalIdxI, nCompIdx))
                            * ( v_n_  - 1./densityNWI) / incp;

                    if (problem().variables().dv_dp(globalIdxI) > 0) // this is not physically possible!!
                    {
                        p_ -= 2*incp;
                        Scalar v_n_ = 1. / FluidSystem::density(pseudoFluidState, nPhaseIdx);
                        problem().variables().dv_dp(globalIdxI) = (problem().variables().totalConcentration(globalIdxI, wCompIdx)
                                + problem().variables().totalConcentration(globalIdxI, nCompIdx))
                                * ( v_n_  - 1./densityNWI) / -incp;
                    }
                }
            }    //end calculation of 1p compress_term
            Scalar compress_term = problem().variables().dv_dp(globalIdxI) / timestep_;

            this->A_[globalIdxI][globalIdxI] -= compress_term*volume;
            this->f_[globalIdxI] -= problem().variables().pressure()[globalIdxI] * compress_term * volume;

            if (isnan(compress_term) || isinf(compress_term))
                DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << globalIdxI);

            if(!GET_PROP_VALUE(TypeTag, EnableCompressibility))
                DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
        }

        // error reduction routine: volumetric error is damped and inserted to right hand side
        // if damping is not done, the solution method gets unstable!
        problem().variables().volErr()[globalIdxI] /= timestep_;
        Scalar erri = fabs(problem().variables().volErr()[globalIdxI]);
        Scalar x_lo = this->ErrorTermLowerBound_;
        Scalar x_mi = this->ErrorTermUpperBound_;
        Scalar fac  = this->ErrorTermFactor_;
        if (pressureType == pw)
            fac = 0.1*this->ErrorTermFactor_;
        Scalar lofac = 0.;
//            lofac/=fac;
        Scalar hifac = 0.;
//        hifac /= fac;

        if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxErr) && (!problem().timeManager().willBeFinished()))
        {
            if (erri <= x_mi * maxErr)
                this->f_[globalIdxI] +=
                        problem().variables().errorCorrection(globalIdxI) =
                                fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/maxErr)
                                    * problem().variables().volErr()[globalIdxI] * volume;
            else
                this->f_[globalIdxI] +=
                        problem().variables().errorCorrection(globalIdxI) =
                                fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr)
                                    * problem().variables().volErr()[globalIdxI] * volume;
        }
    } // end grid traversal
    return;
}

//! constitutive functions are updated once if new concentrations are calculated and stored in the variables container
/*!
 * In contrast to the standard sequential 2p2c model, this method also holds routines
 * to adapt the subdomain. The subdomain indicates weather we are in 1p domain (value = 1)
 * or in the two phase subdomain (value = 2).
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::updateMaterialLaws()
{
    problem().variables().communicateTransportedQuantity();
    problem().variables().communicatePressure();

    // instantiate standard 2p2c and pseudo1p fluid state objects
    FluidState fluidState;
    PseudoOnePTwoCFluidState<TypeTag> pseudoFluidState;

    //get timestep for error term
    Scalar dt = problem().timeManager().timeStepSize();

    // next subdomain map
    if (problem().timeManager().time() == 0.)
        nextSubdomain = 2;  // start with complicated sub in initialization
    else
        nextSubdomain = 1;  // reduce complexity after first TS

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        const typename ElementIterator::Entity::Geometry &geo = eIt->geometry();

        // get global coordinate of cell center
        GlobalPosition globalPos = geo.center();

        int globalIdx = problem().variables().index(*eIt);

        Scalar temperature_ = problem().temperatureAtPos(globalPos);
        // reset volErr
        problem().variables().volErr()[globalIdx] = 0;

        // get the overall mass of component 1:  Z1 = C^k / (C^1+C^2) [-]
        Scalar Z1 = problem().variables().totalConcentration(globalIdx, wCompIdx)
                / (problem().variables().totalConcentration(globalIdx, wCompIdx)
                        + problem().variables().totalConcentration(globalIdx, nCompIdx));
        // make shure only physical quantities enter flash calculation
        #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
        if(Z1<0. || Z1 > 1.)
        {
            std::cout << "Feed mass fraction unphysical: Z1 = " << Z1
                   << " at global Idx " << globalIdx
                   << " , because totalConcentration(globalIdx, wCompIdx) = "
                   << problem().variables().totalConcentration(globalIdx, wCompIdx)
                   << " and totalConcentration(globalIdx, nCompIdx) = "
                   << problem().variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            if(Z1<0.)
            {
            Z1 = 0.;
            // add this error to volume error term for correction in next TS
            problem().variables().volErr()[globalIdx] +=
                    problem().variables().totalConcentration(globalIdx, wCompIdx)
                    / problem().variables().densityWetting(globalIdx);
            //regul!
            problem().variables().totalConcentration(globalIdx, wCompIdx) = 0.;
            Dune::dgrave << "Regularize totalConcentration(globalIdx, wCompIdx) = "
                << problem().variables().totalConcentration(globalIdx, wCompIdx)<< std::endl;
            }
        else
            {
            Z1 = 1.;
            // add this error to volume error term for correction in next TS
            problem().variables().volErr()[globalIdx] +=
                    problem().variables().totalConcentration(globalIdx, nCompIdx)
                    / problem().variables().densityNonwetting(globalIdx);
            //regul!
            problem().variables().totalConcentration(globalIdx, nCompIdx) = 0.;
            Dune::dgrave << "Regularize totalConcentration(globalIdx, nCompIdx) = "
                << problem().variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            }
        }
        #endif

        if (problem().variables().subdomain(globalIdx)==2)   //=> 2p domain
        {
            //determine phase pressures from primary pressure variable and pc of last TS
            PhaseVector pressure(0.);
            switch (pressureType)
            {
            case pw:
            {
                pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx];
                pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx]
                          + problem().variables().capillaryPressure(globalIdx);
                break;
            }
            case pn:
            {
                pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx]
                         - problem().variables().capillaryPressure(globalIdx);
                pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx];
                break;
            }
            }

            //complete fluid state
            fluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, *eIt), temperature_);

            // iterations part in case of enabled capillary pressure
            Scalar pc(0.), oldPc(0.);
            if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
            {
                pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, *eIt),
                        fluidState.saturation(wPhaseIdx));
                int maxiter = 5; int iterout = -1;
                //start iteration loop
                for(int iter=0; iter < maxiter; iter++)
                {
                    //prepare pressures to enter flash calculation
                    switch (pressureType)
                    {
                    case pw:
                    {
                        // pressure[w] does not change, since it is primary variable
                        pressure[nPhaseIdx] = pressure[wPhaseIdx] + pc;
                        break;
                    }
                    case pn:
                    {
                        // pressure[n] does not change, since it is primary variable
                        pressure[wPhaseIdx] = pressure[nPhaseIdx] - pc;
                        break;
                    }
                    }

                    //store old pc
                    oldPc = pc;
                    //update with better pressures
                    fluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, *eIt),
                                        problem().temperatureAtPos(globalPos));
                    pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, *eIt),
                                        fluidState.saturation(wPhaseIdx));
                    // TODO: get right criterion, do output for evaluation
                    //converge criterion
                    if (abs(oldPc-pc)<10)
                        maxiter = 1;
                    iterout = iter;
                }
                if(iterout !=0)
                Dune::dinfo << iterout << "times iteration of pc was applied at Idx " << globalIdx << ", pc delta still " << abs(oldPc-pc) << std::endl;
            }

            /******** update variables in variableclass **********/
            // initialize saturation
            problem().variables().saturation(globalIdx) = fluidState.saturation(wPhaseIdx);

            // the following makes shure pC is neglected
            problem().variables().capillaryPressure(globalIdx) = pc; // pc=0 if EnableCapillarity=false

            // initialize viscosities
            problem().variables().viscosityWetting(globalIdx)
                    = FluidSystem::viscosity(fluidState, wPhaseIdx);
            problem().variables().viscosityNonwetting(globalIdx)
                    = FluidSystem::viscosity(fluidState, nPhaseIdx);

            // initialize mobilities
            problem().variables().mobilityWetting(globalIdx) =
                    MaterialLaw::krw(problem().spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                        / problem().variables().viscosityWetting(globalIdx);
            problem().variables().mobilityNonwetting(globalIdx) =
                    MaterialLaw::krn(problem().spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                        / problem().variables().viscosityNonwetting(globalIdx);

            // initialize mass fractions
            problem().variables().wet_X1(globalIdx) = fluidState.massFraction(wPhaseIdx, wCompIdx);
            problem().variables().nonwet_X1(globalIdx) = fluidState.massFraction(nPhaseIdx, wCompIdx);

            // initialize densities
            problem().variables().densityWetting(globalIdx)
                    = FluidSystem::density(fluidState, wPhaseIdx);
            problem().variables().densityNonwetting(globalIdx)
                    = FluidSystem::density(fluidState, nPhaseIdx);

            // determine volume mismatch between actual fluid volume and pore volume
            Scalar sumConc = (problem().variables().totalConcentration(globalIdx, wCompIdx)
                    + problem().variables().totalConcentration(globalIdx, nCompIdx));
            Scalar massw = problem().variables().numericalDensity(globalIdx, wPhaseIdx) = sumConc * fluidState.phaseMassFraction(wPhaseIdx);
            Scalar massn = problem().variables().numericalDensity(globalIdx, nPhaseIdx) =sumConc * fluidState.phaseMassFraction(nPhaseIdx);
            Scalar vol = massw / problem().variables().densityWetting(globalIdx)
                       + massn / problem().variables().densityNonwetting(globalIdx);
            if (dt != 0)
            {
                problem().variables().volErr()[globalIdx] += (vol - problem().spatialParameters().porosity(globalPos, *eIt));

                Scalar volErrI = problem().variables().volErr(globalIdx);
                if (std::isnan(volErrI))
                {
                    DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate:\n"
                            << "volErr[" << globalIdx << "] isnan: vol = " << vol
                            << ", massw = " << massw << ", rho_l = " << problem().variables().densityWetting(globalIdx)
                            << ", massn = " << massn << ", rho_g = " << problem().variables().densityNonwetting(globalIdx)
                            << ", poro = " << problem().spatialParameters().porosity(globalPos, *eIt) << ", dt = " << dt);
                }
            }
            else
                problem().variables().volErr()[globalIdx] = 0;

            // check subdomain consistency
            timer_.start();
            if (problem().variables().saturation(globalIdx) != (1. || 0.)) // cell still 2p
            {
                // mark this element
                nextSubdomain[globalIdx] = 2;

                // mark neighbors
                IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt!=isItEnd; ++isIt)
                {
                    if (isIt->neighbor())
                    {
                        int globalIdxJ = problem().variables().index(*(isIt->outside()));
                        // mark neighbor Element
                        nextSubdomain[globalIdxJ] = 2;
                    }
                }
                }
            timer_.stop();
            // end subdomain check
        }
        else    // simple
        {
            // Todo: only variables of present phase need to be updated
            Scalar press1p = problem().variables().pressure()[globalIdx];
            pseudoFluidState.update(Z1, press1p, problem().variables().saturation(globalIdx), problem().temperatureAtPos(globalPos));

            // initialize viscosities
            problem().variables().viscosityWetting(globalIdx)
                    = FluidSystem::viscosity(pseudoFluidState, wPhaseIdx);
            problem().variables().viscosityNonwetting(globalIdx)
                    = FluidSystem::viscosity(pseudoFluidState, nPhaseIdx);

            // initialize mobilities
            problem().variables().mobilityWetting(globalIdx) =
                    MaterialLaw::krw(problem().spatialParameters().materialLawParams(globalPos, *eIt), pseudoFluidState.saturation(wPhaseIdx))
                        / problem().variables().viscosityWetting(globalIdx);
            problem().variables().mobilityNonwetting(globalIdx) =
                    MaterialLaw::krn(problem().spatialParameters().materialLawParams(globalPos, *eIt), pseudoFluidState.saturation(wPhaseIdx))
                        / problem().variables().viscosityNonwetting(globalIdx);

            // initialize mass fractions
            problem().variables().wet_X1(globalIdx) = pseudoFluidState.massFraction(wPhaseIdx, wCompIdx);
            problem().variables().nonwet_X1(globalIdx) = pseudoFluidState.massFraction(nPhaseIdx, wCompIdx);

            // initialize densities
            problem().variables().densityWetting(globalIdx)
                    = FluidSystem::density(pseudoFluidState, wPhaseIdx);
            problem().variables().densityNonwetting(globalIdx)
                    = FluidSystem::density(pseudoFluidState, nPhaseIdx);

            // error term handling
            Scalar sumConc = (problem().variables().totalConcentration(globalIdx, wCompIdx)
                    + problem().variables().totalConcentration(globalIdx, nCompIdx));
            Scalar vol(0.);
            if(problem().variables().saturation(globalIdx) == 1.)    //only w-phase
                vol = sumConc / problem().variables().densityWetting(globalIdx);
            else    //only nw-phase
                vol = sumConc / problem().variables().densityNonwetting(globalIdx);

            if (dt != 0)
                problem().variables().volErr()[globalIdx] = (vol - problem().spatialParameters().porosity(globalPos, *eIt));


        }
    }// end grid traversal

    problem().variables().subdomain() = nextSubdomain;
    Dune::dinfo << "Subdomain routines took " << timer_.elapsed() << " seconds" << std::endl;

    return;
}

}//end namespace Dumux
#endif
