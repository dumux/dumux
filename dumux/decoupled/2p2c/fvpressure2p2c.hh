// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
#ifndef DUMUX_FVPRESSURE2P2C_HH
#define DUMUX_FVPRESSURE2P2C_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include <dumux/decoupled/2p2c/fvpressurecompositional.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Benjamin Faigle, Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{
//! The finite volume model for the solution of the compositional pressure equation
/*! \ingroup multiphase
 *  Provides a Finite Volume implementation for the pressure equation of a gas-liquid
 *  system with two components. An IMPES-like method is used for the sequential
 *  solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 *  Isothermal conditions and local thermodynamic
 *  equilibrium are assumed.  Gravity is included.
 *  \f[
         c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
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
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2P2C
: public FVPressureCompositional<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
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
    typedef typename GridView::Intersection Intersection;

    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<Scalar, 2> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    // the typenames used for the stiffness matrix and solution vector
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) RHSVector;

protected:
    Problem& problem()
    {
        return problem_;
    }
    const Problem& problem() const
    {
        return problem_;
    }

public:
    void getSource(Dune::FieldVector<Scalar, 2>&, const Element&, const CellData&, const bool);

    void getStorage(Dune::FieldVector<Scalar, 2>&, const Element&, const CellData&, const bool);

    void getFlux(Dune::FieldVector<Scalar, 2>&, const Intersection&, const CellData&, const bool);

    void getFluxOnBoundary(Dune::FieldVector<Scalar, 2>&,
                            const Intersection&, const CellData&, const bool);

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();
    //updates secondary variables for one cell and stores in the variables object
    void updateMaterialLawsInElement(const Element&);

#if 0
    //the variables object is initialized, non-compositional before and compositional after first pressure calculation
    void initialMaterialLaws(bool compositional);

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    //numerical volume derivatives wrt changes in mass, pressure
    void volumeDerivatives(GlobalPosition globalPos, ElementPointer ep, Scalar& dv_dC1, Scalar& dv_dC2, Scalar& dv_dp);
#endif
    //! Constructs a FVPressure2P2C object
    /**
     * \param problem a problem class object
     */
    FVPressure2P2C(Problem& problem) : FVPressureCompositional<TypeTag>(problem),
        problem_(problem)
    {
        ErrorTermFactor_ = GET_PARAM(TypeTag, Scalar, ErrorTermFactor);
        ErrorTermLowerBound_ = GET_PARAM(TypeTag, Scalar, ErrorTermLowerBound);
        ErrorTermUpperBound_ = GET_PARAM(TypeTag, Scalar, ErrorTermUpperBound);

        maxError_=0.;
        if (pressureType != pw && pressureType != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
    }

protected:
    Problem& problem_;
    Scalar maxError_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening
    static constexpr int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
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
void FVPressure2P2C<TypeTag>::getSource(Dune::FieldVector<Scalar, 2>& sourceEntry, const Element& elementI, const CellData& cellDataI, const bool first)
{
    sourceEntry=0.;

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementI.geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementI.geometry().volume();

    Scalar densityWI = cellDataI.density(wPhaseIdx);
    Scalar densityNWI = cellDataI.density(nPhaseIdx);

    /****************implement sources************************/
    PrimaryVariables source(NAN);
    problem().source(source, elementI);
    if(first)
    {
            source[contiWEqIdx] /= densityWI;
            source[contiNEqIdx] /= densityNWI;
    }
    else
    {
            // derivatives of the fluid volume with respect to concentration of components, or pressure
            if (!cellDataI.hasVolumeDerivatives())
                this->volumeDerivatives(globalPos, elementI);

            source[contiWEqIdx] *= cellDataI.dv(wCompIdx);        // note: dV_[i][1] = dv_dC1 = dV/dm1
            source[contiNEqIdx] *= cellDataI.dv(nCompIdx);
    }
    sourceEntry[1] = volume * (source[contiWEqIdx] + source[contiNEqIdx]);

    return;
}

template<class TypeTag>
void FVPressure2P2C<TypeTag>::getStorage(Dune::FieldVector<Scalar, 2>& storageEntry, const Element& elementI, const CellData& cellDataI,
                                                                                const bool first)
{
    storageEntry = 0.;
    // cell index
    int globalIdxI = problem().variables().index(elementI);
    Scalar volume = elementI.geometry().volume();

    storageEntry = 0.;
    // determine maximum error to scale error-term
    Scalar timestep_ = problem().timeManager().timeStepSize();
//    maxErr /= timestep_;

    // compressibility term
    if (!first && timestep_ != 0)
    {
        Scalar compress_term = cellDataI.dv_dp() / timestep_;

        storageEntry[0] -= compress_term*volume;
        storageEntry[1] -= this->pressure(globalIdxI) * compress_term * volume;

        if (isnan(compress_term) || isinf(compress_term))
            DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << globalIdxI);

        if(!GET_PROP_VALUE(TypeTag, EnableCompressibility))
            DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
    }


    // error reduction routine: volumetric error is damped and inserted to right hand side
    // if damping is not done, the solution method gets unstable!
    problem().variables().cellData(globalIdxI).volumeError() /= timestep_;
    Scalar maxError = this->maxError_;
    Scalar erri = fabs(cellDataI.volumeError());
    Scalar x_lo = ErrorTermLowerBound_;
    Scalar x_mi = ErrorTermUpperBound_;
    Scalar fac  = ErrorTermFactor_;
    if (pressureType == pw)
        fac = 0.1*ErrorTermFactor_;
    Scalar lofac = 0.;
    Scalar hifac = 0.;

    if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxError) && (!problem().timeManager().willBeFinished()))
    {
        if (erri <= x_mi * maxError)
            storageEntry[1] +=
                    problem().variables().cellData(globalIdxI).errorCorrection() =
                            fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/maxError)
                                * cellDataI.volumeError() * volume;
        else
            storageEntry[1] +=
                    problem().variables().cellData(globalIdxI).errorCorrection() =
                            fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxError)
                                * cellDataI.volumeError() * volume;
    }

    return;
}



template<class TypeTag>
void FVPressure2P2C<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI, const bool first)
{
    entries = 0.;
    ElementPointer elementPointerI = intersection.inside();

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // cell index
    int globalIdxI = problem().variables().index(*elementPointerI);

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementPointerI->geometry().volume();
    Scalar perimeter = cellDataI.perimeter();

    const GlobalPosition& gravity_ = problem().gravity();

    Scalar densityWI = cellDataI.density(wPhaseIdx);
    Scalar densityNWI = cellDataI.density(nPhaseIdx);

    // get absolute permeability
    FieldMatrix permeabilityI(problem().spatialParameters().intrinsicPermeability(globalPos, *elementPointerI));

    // get mobilities and fractional flow factors
    Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
    Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);
    Scalar fractionalWI=0, fractionalNWI=0;
    if (first)
    {
        fractionalWI = lambdaWI / (lambdaWI+ lambdaNWI);
        fractionalNWI = lambdaNWI / (lambdaWI+ lambdaNWI);
    }

    Scalar pcI = cellDataI.capillaryPressure();

    // get normal vector
    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();

    // get face volume
    Scalar faceArea = intersection.geometry().volume();

        // access neighbor
        ElementPointer neighborPointer = intersection.outside();
        int globalIdxJ = problem().variables().index(*neighborPointer);
        CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

        // gemotry info of neighbor
        const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

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

        // get densities & mobilities in neighbor
        Scalar densityWJ = cellDataJ.density(wPhaseIdx);
        Scalar densityNWJ = cellDataJ.density(nPhaseIdx);
        Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
        Scalar lambdaNWJ = cellDataJ.mobility(nPhaseIdx);

        Scalar pcJ = cellDataJ.capillaryPressure();

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

            entries[0] = fabs(lambda*faceArea*(permeability*unitOuterNormal)/(dist));

            Scalar factor = (fractionalWI + fractionalWJ) * (rhoMeanW) * 0.5 + (fractionalNWI + fractionalNWJ) * (rhoMeanNW) * 0.5;
            entries[1] = factor * lambda * faceArea * (permeability * gravity_) * (unitOuterNormal * unitDistVec);
        }
        else
        {
            // determine volume derivatives
            if (!cellDataJ.hasVolumeDerivatives())
                this->volumeDerivatives(globalPosNeighbor, *neighborPointer);

            Scalar dv_dC1 = (cellDataJ.dv(wPhaseIdx)
                        + cellDataI.dv(wPhaseIdx)) / 2; // dV/dm1= dv/dC^1
            Scalar dv_dC2 = (cellDataJ.dv(nPhaseIdx)
                        + cellDataI.dv(nPhaseIdx)) / 2;

            Scalar graddv_dC1 = (cellDataJ.dv(wPhaseIdx)
                                    - cellDataI.dv(wPhaseIdx)) / dist;
            Scalar graddv_dC2 = (cellDataJ.dv(nPhaseIdx)
                                    - cellDataI.dv(nPhaseIdx)) / dist;


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

            switch (pressureType)
            {
            case pw:
            {
                potentialW = (unitOuterNormal * unitDistVec) * (this->pressure(globalIdxI)
                        - this->pressure(globalIdxJ)) / dist;
                potentialNW = (unitOuterNormal * unitDistVec) * (this->pressure(globalIdxI)
                        - this->pressure(globalIdxJ) + pcI - pcJ) / dist;
                break;
            }
            case pn:
            {
                potentialW = (unitOuterNormal * unitDistVec) * (this->pressure(globalIdxI)
                        - this->pressure(globalIdxJ) - pcI + pcJ) / dist;
                potentialNW = (unitOuterNormal * unitDistVec) * (this->pressure(globalIdxI)
                        - this->pressure(globalIdxJ)) / dist;
                break;
            }
            }

            potentialW += densityW * (unitDistVec * gravity_);
            potentialNW += densityNW * (unitDistVec * gravity_);

            // initialize convenience shortcuts
            Scalar lambdaW, lambdaN;
            Scalar dV_w(0.), dV_n(0.);        // dV_a = \sum_k \rho_a * dv/dC^k * X^k_a
            Scalar gV_w(0.), gV_n(0.);        // multipaper eq(3.3) line 3 analogon dV_w


            //do the upwinding of the mobility depending on the phase potentials
            if (potentialW >= 0.)
            {
                dV_w = (dv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
                lambdaW = cellDataI.mobility(wPhaseIdx);
                gV_w = (graddv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                        + graddv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
                dV_w *= densityWI; gV_w *= densityWI;
            }
            else
            {
                dV_w = (dv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
                lambdaW = cellDataJ.mobility(wPhaseIdx);
                gV_w = (graddv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                        + graddv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
                dV_w *= densityWJ; gV_w *= densityWJ;
            }
            if (potentialNW >= 0.)
            {
                dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
                lambdaN = cellDataI.mobility(nPhaseIdx);
                gV_n = (graddv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                        + graddv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
                dV_n *= densityNWI; gV_n *= densityNWI;
            }
            else
            {
                dV_n = (dv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
                lambdaN = cellDataJ.mobility(nPhaseIdx);
                gV_n = (graddv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                        + graddv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
                dV_n *= densityNWJ; gV_n *= densityNWJ;
            }

            //calculate current matrix entry
            entries[0] = faceArea * (lambdaW * dV_w + lambdaN * dV_n);
            entries[0] -= volume * faceArea / perimeter * (lambdaW * gV_w + lambdaN * gV_n);     // = boundary integral - area integral
            entries[0] *= fabs((permeability*unitOuterNormal)/(dist));

            //calculate right hand side
            entries[1] = faceArea  * (unitOuterNormal * unitDistVec) * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);
            entries[1] -= volume * faceArea / perimeter * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n);
            entries[1] *= (permeability * gravity_);         // = multipaper eq(3.3) line 2+3

            // include capillary pressure fluxes
            switch (pressureType)
            {
            case pw:
                {
                    // calculate capillary pressure gradient
                    Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                    pCGradient *= (pcI - pcJ) / dist;

                    //add capillary pressure term to right hand side
                    entries[1] += lambdaN * dV_n * (permeability * pCGradient) * faceArea
                                 - lambdaN * gV_n * (permeability * pCGradient) * volume * faceArea / perimeter;
                    break;
                }
            case pn:
                {
                    // calculate capillary pressure gradient
                    Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                    pCGradient *= (pcI - pcJ) / dist;

                    //add capillary pressure term to right hand side
                    entries[1] -= lambdaW * dV_w * (permeability * pCGradient) * faceArea
                                    - lambdaW * gV_w * (permeability * pCGradient) * volume * faceArea / perimeter;
                    break;
                }
            }
        }   // end !first
}

template<class TypeTag>
void FVPressure2P2C<TypeTag>::getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI, const bool first)
{
    entries = 0.;
    // get global coordinate of cell center
    ElementPointer elementPointerI = intersection.inside();
    const GlobalPosition& globalPos = elementPointerI->geometry().center();
    int globalIdxI = problem().variables().index(*elementPointerI);

    // get normal vector
    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();
    // get face volume
    Scalar faceArea = intersection.geometry().volume();
                // get volume derivatives inside the cell
                Scalar dv_dC1 = cellDataI.dv(wCompIdx);
                Scalar dv_dC2 = cellDataI.dv(nCompIdx);
                Scalar densityWI = cellDataI.density(wPhaseIdx);
                Scalar densityNWI = cellDataI.density(nPhaseIdx);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = intersection.geometry().center();

                // geometrical information
                GlobalPosition distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary condition for boundary face center
                BoundaryTypes bcType;
                problem().boundaryTypes(bcType, intersection);

                // prepare pressure boundary condition
                PhaseVector pressBC(0.);
                Scalar pcBound (0.);

                /**********         Dirichlet Boundary        *************/
                if (bcType.isDirichlet(Indices::pressureEqIdx))
                {
                    // get absolute permeability
                    FieldMatrix permeabilityI(problem().spatialParameters().intrinsicPermeability(globalPos, *elementPointerI));
                    const GlobalPosition& gravity_ = problem().gravity();
                    // get mobilities and fractional flow factors
                    Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
                    Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);

                    //permeability vector at boundary
                    Dune::FieldVector<Scalar, dim> permeability(0);
                    permeabilityI.mv(unitDistVec, permeability);

                    // create a fluid state for the boundary
                    FluidState BCfluidState;

                    //read boundary values
                    PrimaryVariables primaryVariablesOnBoundary(NAN);
                    problem().dirichlet(primaryVariablesOnBoundary, intersection);

                    if (first)
                    {

                        Scalar fractionalWI=0, fractionalNWI=0;
                        fractionalWI = lambdaWI / (lambdaWI+ lambdaNWI);
                        fractionalNWI = lambdaNWI / (lambdaWI+ lambdaNWI);

                        Scalar lambda = lambdaWI+lambdaNWI;
                        entries[0] += lambda * faceArea * (permeability * unitOuterNormal) / (dist);
                        Scalar pressBC = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                        entries[1] += lambda * faceArea * pressBC * (permeability * unitOuterNormal) / (dist);
                        Scalar rightentry = (fractionalWI * densityWI
                                             + fractionalNWI * densityNWI)
                                             * lambda * faceArea * (unitOuterNormal * unitDistVec) *(permeability*gravity_);
                        entries[1] -= rightentry;
                    }
                    else    //not first
                    {
                        // read boundary values
                        problem().transportModel().evalBoundary(globalPosFace,
                                                                    intersection,
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
                                    problem().spatialParameters().materialLawParams(globalPos, *elementPointerI), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityWBound;
                            lambdaNWBound = MaterialLaw::krn(
                                    problem().spatialParameters().materialLawParams(globalPos, *elementPointerI), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityNWBound;
                            break;
                            }
                        }

                        Scalar rhoMeanW = 0.5 * (densityWI + densityWBound);
                        Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWBound);

                        Scalar potentialW = 0;
                        Scalar potentialNW = 0;

                        Scalar densityW = 0;
                        Scalar densityNW = 0;

                        if (!first)
                        {
//                            potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                            potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
//
//                            // do potential upwinding according to last potGradient vs Jochen: central weighting
//                            densityW = (potentialW > 0.) ? densityWI : densityWBound;
//                            densityNW = (potentialNW > 0.) ? densityNWI : densityNWBound;
//
//                            densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                            densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                            densityW=rhoMeanW; densityNW=rhoMeanNW;

                            //calculate potential gradient
                            if (pressureType == pw)
                            {
                                potentialW = ((unitOuterNormal * distVec)  / (dist * dist)) *
                                        (this->pressure(globalIdxI) - pressBC[wPhaseIdx]);
                                potentialNW = ((unitOuterNormal * distVec)  / (dist * dist)) *
                                        (this->pressure(globalIdxI) + cellDataI.capillaryPressure()
                                                - pressBC[wPhaseIdx] - pcBound);
                            }
                            else if (pressureType == pn)
                            {
                                potentialW = ((unitOuterNormal * distVec) / (dist * dist)) *
                                        (this->pressure(globalIdxI) - cellDataI.capillaryPressure() - pressBC[nPhaseIdx] + pcBound);
                                potentialNW = ((unitOuterNormal * distVec) / (dist * dist)) *
                                        (this->pressure(globalIdxI) - pressBC[nPhaseIdx]);
                            }
                            potentialW += densityW * (unitDistVec * gravity_);
                            potentialNW += densityNW * (unitDistVec * gravity_);
                       }   //end !first

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW, lambdaNW;
                        Scalar dV_w, dV_n;     // gV_a weglassen, da dV/dc am Rand ortsunabhängig angenommen -> am rand nicht bestimmbar -> nur Randintegral ohne Gebietsintegral

                        if (potentialW >= 0.)
                        {
                            densityW = (potentialW == 0) ? rhoMeanW : densityWI;
                            dV_w = (dv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                                               + dv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
                            dV_w *= densityW;
                            lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaWI;
                        }
                        else
                        {
                            densityW = densityWBound;
                            dV_w = (dv_dC1 * BCfluidState.massFraction(wPhaseIdx, wCompIdx)
                                     + dv_dC2 * BCfluidState.massFraction(wPhaseIdx, nCompIdx));
                            dV_w *= densityW;
                            lambdaW = lambdaWBound;
                        }
                        if (potentialNW >= 0.)
                        {
                            densityNW = (potentialNW == 0) ? rhoMeanNW : densityNWI;
                            dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                                    + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
                            dV_n *= densityNW;
                            lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNWI;
                        }
                        else
                        {
                            densityNW = densityNWBound;
                            dV_n = (dv_dC1 * BCfluidState.massFraction(nPhaseIdx, wCompIdx)
                                    + dv_dC2 * BCfluidState.massFraction(nPhaseIdx, nCompIdx));
                            dV_n *= densityNW;
                            lambdaNW = lambdaNWBound;
                        }

                        //calculate current matrix entry
                        Scalar entry = (lambdaW * dV_w + lambdaNW * dV_n) * ((permeability * unitDistVec) / dist) * faceArea
                                * (unitOuterNormal * unitDistVec);

                        //calculate right hand side
                        Scalar rightEntry = (lambdaW * densityW * dV_w + lambdaNW * densityNW * dV_n) * (permeability * gravity_)
                                * faceArea ;


                        // include capillary pressure fluxes
                        switch (pressureType)
                        {
                        case pw:
                            {
                                // calculate capillary pressure gradient
                                Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                                pCGradient *= (cellDataI.capillaryPressure() - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                rightEntry += lambdaNW * dV_n * (permeability * pCGradient) * faceArea;
                                break;
                            }
                        case pn:
                            {
                                // calculate capillary pressure gradient
                                Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
                                pCGradient *= (cellDataI.capillaryPressure() - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                rightEntry -= lambdaW * dV_w * (permeability * pCGradient) * faceArea;
                                break;
                            }
                        }


                        // set diagonal entry and right hand side entry
                        entries[0] += entry;
                        entries[1] += entry * primaryVariablesOnBoundary[Indices::pressureEqIdx];
                        entries[1] -= rightEntry * (unitOuterNormal * unitDistVec);
                    }    //end of if(first) ... else{...
                }   // end dirichlet

                /**********************************
                 * set neumann boundary condition
                 **********************************/
                else if(bcType.isNeumann(Indices::pressureEqIdx))
                {
                    PrimaryVariables J(NAN);
                    problem().neumann(J, intersection);
                    if (first)
                    {
                        J[contiWEqIdx] /= densityWI;
                        J[contiNEqIdx] /= densityNWI;
                    }
                    else
                    {
                        J[contiWEqIdx] *= dv_dC1;
                        J[contiNEqIdx] *= dv_dC2;
                    }

                    entries[1] -= (J[contiWEqIdx] + J[contiNEqIdx]) * faceArea;
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Boundary Condition neither Dirichlet nor Neumann!");


    return;
}

//! constitutive functions are updated once if new concentrations are calculated and stored in the variables container
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLaws()
{

    // instantiate a brandnew fluid state object
    FluidState fluidState;

    Scalar maxError = 0.;
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem().variables().index(*eIt);

        CellData& cellData = problem().variables().cellData(globalIdx);

        this->updateMaterialLawsInElement(*eIt);

        maxError = std::max(maxError, fabs(cellData.volumeError()));
    }
    maxError_ = maxError/problem().timeManager().timeStepSize();
    return;
}

template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLawsInElement(const Element& elementI)
{
    // instantiate a brandnew fluid state object
    FluidState fluidState;

    // get global coordinate of cell center
    GlobalPosition globalPos = elementI.geometry().center();

    int globalIdx = problem().variables().index(elementI);
    CellData& cellData = problem().variables().cellData(globalIdx);

    Scalar temperature_ = problem().temperatureAtPos(globalPos);
    // reset
    cellData.reset();

    // make shure total concentrations from solution vector are exact in fluidstate
    fluidState.setMassConcentration(wCompIdx,
            problem().transportModel().transportedQuantity()[wCompIdx][globalIdx]);
    fluidState.setMassConcentration(nCompIdx,
            problem().transportModel().transportedQuantity()[nCompIdx][globalIdx]);
    // get the overall mass of component 1 Z1 = C^k / (C^1+C^2) [-]
    Scalar Z1 = fluidState.massConcentration(wCompIdx)
            / (fluidState.massConcentration(wCompIdx)
                    + fluidState.massConcentration(nCompIdx));

    // make shure only physical quantities enter flash calculation
    #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
    if(Z1<0. || Z1 > 1.)
    {
        Dune::dgrave << "Feed mass fraction unphysical: Z1 = " << Z1
               << " at global Idx " << globalIdx
               << " , because totalConcentration(wCompIdx) = "
               << cellData.totalConcentration(wCompIdx)
               << " and totalConcentration(nCompIdx) = "
               << cellData.totalConcentration(nCompIdx)<< std::endl;
        if(Z1<0.)
            {
            Z1 = 0.;
            cellData.setTotalConcentration(wCompIdx, 0.);
            Dune::dgrave << "Regularize totalConcentration(wCompIdx) = "
                << cellData.totalConcentration(wCompIdx)<< std::endl;
            }
        else
            {
            Z1 = 1.;
            cellData.setTotalConcentration(nCompIdx, 0.);
            Dune::dgrave << "Regularize totalConcentration(globalIdx, nCompIdx) = "
                << cellData.totalConcentration(nCompIdx)<< std::endl;
            }
    }
    #endif

    //determine phase pressures from primary pressure variable and pc of last TS
    PhaseVector pressure(0.);
    switch (pressureType)
    {
    case pw:
    {
        pressure[wPhaseIdx] = this->pressure(globalIdx);
        pressure[nPhaseIdx] = this->pressure(globalIdx)
                  + cellData.capillaryPressure();
        break;
    }
    case pn:
    {
        pressure[wPhaseIdx] = this->pressure(globalIdx)
                 - cellData.capillaryPressure();
        pressure[nPhaseIdx] = this->pressure(globalIdx);
        break;
    }
    }

    //complete fluid state
    fluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, elementI), temperature_);

    // iterations part in case of enabled capillary pressure
    Scalar pc(0.), oldPc(0.);
    if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
    {
        pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, elementI),
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
            fluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, elementI),
                                problem().temperatureAtPos(globalPos));
            pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, elementI),
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

    // secondary variables
    // initialize saturation, capillary pressure
    cellData.setFluidState(fluidState);

    // initialize phase properties not stored in fluidstate
    cellData.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState, wPhaseIdx));
    cellData.setViscosity(nPhaseIdx, FluidSystem::viscosity(fluidState, nPhaseIdx));

    // initialize mobilities
    cellData.setMobility(wPhaseIdx, MaterialLaw::krw(problem().spatialParameters().materialLawParams(globalPos, elementI), fluidState.saturation(wPhaseIdx))
                / cellData.viscosity(wPhaseIdx));
    cellData.setMobility(nPhaseIdx, MaterialLaw::krn(problem().spatialParameters().materialLawParams(globalPos, elementI), fluidState.saturation(wPhaseIdx))
                / cellData.viscosity(nPhaseIdx));

    // determine volume mismatch between actual fluid volume and pore volume
    Scalar sumConc = (cellData.totalConcentration(wCompIdx)
            + cellData.totalConcentration(nCompIdx));
    Scalar massw = cellData.numericalDensity(wPhaseIdx) = sumConc * fluidState.phaseMassFraction(wPhaseIdx);
    Scalar massn = cellData.numericalDensity(nPhaseIdx) = sumConc * fluidState.phaseMassFraction(nPhaseIdx);

    if ((cellData.density(wPhaseIdx)*cellData.density(nPhaseIdx)) == 0)
        DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate: try to divide by 0 density");
    Scalar vol = massw / cellData.density(wPhaseIdx) + massn / cellData.density(nPhaseIdx);
    if (problem().timeManager().timeStepSize() != 0)
    {
        cellData.volumeError()=(vol - problem().spatialParameters().porosity(globalPos, elementI));

        if (std::isnan(cellData.volumeError()))
        {
            DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate:\n"
                    << "volErr[" << globalIdx << "] isnan: vol = " << vol
                    << ", massw = " << massw << ", rho_l = " << cellData.density(wPhaseIdx)
                    << ", massn = " << massn << ", rho_g = " << cellData.density(nPhaseIdx)
                    << ", poro = " << problem().spatialParameters().porosity(globalPos, elementI)
                    << ", dt = " << problem().timeManager().timeStepSize());
        }
    }
    else
        cellData.volumeError()=0.;

    cellData.reset();

    return;
}


#if 0
/*!
 *  \brief initializes the fluid distribution and hereby the variables container
 *
 *  This function equals the method initialguess() and transportInitial() in the old model.
 *  It differs from updateMaterialLaws because there are two possible initial conditions:
 *  saturations and concentration.
 *  \param compositional flag that determines if compositional effects are regarded, i.e.
 *      a reasonable pressure field is known.
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::initialMaterialLaws(bool compositional)
{
//    problem().variables().communicateTransportedQuantity();
//    problem().variables().communicatePressure();

    // initialize the fluid system
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem().gridView().template end<0>();
    ElementIterator eIt = problem().gridView().template begin<0>();
    for (; eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        // assign an Index for convenience
        int globalIdx = problem().variables().index(*eIt);

        // get the temperature
        Scalar temperature_ = problem().temperatureAtPos(globalPos);

        // initial conditions
        PhaseVector pressure(0.);
        Scalar sat_0=0.;

        typename Indices::BoundaryFormulation icFormulation;
        problem().initialFormulation(icFormulation, *eIt);            // get type of initial condition

        if(!compositional) //means that we do the first approximate guess without compositions
        {
            // phase pressures are unknown, so start with an exemplary
            Scalar exemplaryPressure = problem().referencePressure(*eIt);
            pressure[wPhaseIdx] = pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx] = exemplaryPressure;
            problem().variables().capillaryPressure(globalIdx) = 0.;
            if (icFormulation == Indices::saturation)  // saturation initial condition
            {
                sat_0 = problem().initSat(*eIt);
                fluidState.satFlash(sat_0, pressure, problem().spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
            else if (icFormulation == Indices::concentration) // concentration initial condition
            {
                Scalar Z1_0 = problem().initConcentration(*eIt);
                fluidState.update(Z1_0, pressure, problem().spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
        }
        else if(compositional)    //means we regard compositional effects since we know an estimate pressure field
        {
            if (icFormulation == Indices::saturation)  // saturation initial condition
            {
                //get saturation, determine pc
                sat_0 = problem().initSat(*eIt);
                if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
                {
                    problem().variables().capillaryPressure(globalIdx)
                            = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, *eIt),
                                    sat_0);
                }
                else
                    problem().variables().capillaryPressure(globalIdx) = 0.;

                //determine phase pressures from primary pressure variable
                switch (pressureType)
                {
                    case pw:
                    {
                        pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx];
                        pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx] + problem().variables().capillaryPressure(globalIdx);
                        break;
                    }
                    case pn:
                    {
                        pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx] - problem().variables().capillaryPressure(globalIdx);
                        pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx];
                        break;
                    }
                }

                fluidState.satFlash(sat_0, pressure, problem().spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
            else if (icFormulation == Indices::concentration) // concentration initial condition
            {
                Scalar Z1_0 = problem().initConcentration(*eIt);
                // If total concentrations are given at the boundary, saturation is unknown.
                // This may affect pc and hence p_alpha and hence again saturation -> iteration.

                // iterations in case of enabled capillary pressure
                if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
                {
                    //start with pc from last TS
                    Scalar pc(problem().variables().capillaryPressure(globalIdx));

                    int maxiter = 3;
                    //start iteration loop
                    for(int iter=0; iter < maxiter; iter++)
                    {
                        //determine phase pressures from primary pressure variable
                        switch (pressureType)
                        {
                            case pw:
                            {
                                pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx];
                                pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx] + pc;
                                break;
                            }
                            case pn:
                            {
                                pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx] - pc;
                                pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx];
                                break;
                            }
                        }

                        //store old pc
                        Scalar oldPc = pc;
                        //update with better pressures
                        fluidState.update(Z1_0, pressure, problem().spatialParameters().porosity(globalPos, *eIt),
                                            problem().temperatureAtPos(globalPos));
                        pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, *eIt),
                                            fluidState.saturation(wPhaseIdx));
                        // TODO: get right criterion, do output for evaluation
                        //converge criterion
                        if (abs(oldPc-pc)<10)
                            iter = maxiter;

                        pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPos, *eIt),
                                fluidState.saturation(wPhaseIdx));
                    }
                    problem().variables().capillaryPressure(globalIdx) = pc; //complete iteration procedure
                }
                else  // capillary pressure neglected
                {
                    problem().variables().capillaryPressure(globalIdx) = 0.;
                    pressure[wPhaseIdx] = pressure[nPhaseIdx]
                        = problem().variables().pressure()[globalIdx];
                    fluidState.update(Z1_0, pressure, problem().spatialParameters().porosity(globalPos, *eIt), temperature_);
                }
            } //end conc initial condition
        } //end compositional

        // initialize densities
        problem().variables().densityWetting(globalIdx) = FluidSystem::density(fluidState, wPhaseIdx);
        problem().variables().densityNonwetting(globalIdx) = FluidSystem::density(fluidState, nPhaseIdx);

        // initialize mass fractions
        problem().variables().wet_X1(globalIdx) = fluidState.massFraction(wPhaseIdx, wCompIdx);
        problem().variables().nonwet_X1(globalIdx) = fluidState.massFraction(nPhaseIdx, wCompIdx);


        // initialize viscosities
        problem().variables().viscosityWetting(globalIdx) = FluidSystem::viscosity(fluidState, wPhaseIdx);
        problem().variables().viscosityNonwetting(globalIdx) = FluidSystem::viscosity(fluidState, nPhaseIdx);

        // initialize mobilities
        problem().variables().mobilityWetting(globalIdx) = MaterialLaw::krw(problem().spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem().variables().viscosityWetting(globalIdx);
        problem().variables().mobilityNonwetting(globalIdx) = MaterialLaw::krn(problem().spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem().variables().viscosityNonwetting(globalIdx);

        // initialize cell concentration
        problem().variables().totalConcentration(globalIdx, wCompIdx) = fluidState.massConcentration(wCompIdx);
        problem().variables().totalConcentration(globalIdx, nCompIdx) = fluidState.massConcentration(nCompIdx);
        problem().variables().saturation()[globalIdx] = fluidState.saturation(wPhaseIdx);

        // to prevent output errors, set derivatives 0
        problem().variables().dv(globalIdx, wCompIdx) = 0.;     // dv_dC1 = dV/dm1
        problem().variables().dv(globalIdx, nCompIdx) = 0.;     // dv / dC2
        problem().variables().dv_dp(globalIdx) = 0.;      // dv / dp
    }
    return;
}

//! constitutive functions are updated once if new concentrations are calculated and stored in the variables container
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLaws()
{
//    problem().variables().communicateTransportedQuantity();
//    problem().variables().communicatePressure();

    // instantiate a brandnew fluid state object
    FluidState fluidState;

    //get timestep for error term
    Scalar dt = problem().timeManager().timeStepSize();

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        int globalIdx = problem().variables().index(*eIt);

        Scalar temperature_ = problem().temperatureAtPos(globalPos);
        // reset volume error
        problem().variables().volErr()[globalIdx] = 0;


        // get the overall mass of component 1 Z1 = C^k / (C^1+C^2) [-]
        Scalar Z1 = problem().variables().totalConcentration(globalIdx, wCompIdx)
                / (problem().variables().totalConcentration(globalIdx, wCompIdx)
                        + problem().variables().totalConcentration(globalIdx, nCompIdx));

        // make shure only physical quantities enter flash calculation
        #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
        if(Z1<0. || Z1 > 1.)
        {
            Dune::dgrave << "Feed mass fraction unphysical: Z1 = " << Z1
                   << " at global Idx " << globalIdx
                   << " , because totalConcentration(globalIdx, wCompIdx) = "
                   << problem().variables().totalConcentration(globalIdx, wCompIdx)
                   << " and totalConcentration(globalIdx, nCompIdx) = "
                   << problem().variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            if(Z1<0.)
                {
                Z1 = 0.;
                problem().variables().totalConcentration(globalIdx, wCompIdx) = 0.;
                Dune::dgrave << "Regularize totalConcentration(globalIdx, wCompIdx) = "
                    << problem().variables().totalConcentration(globalIdx, wCompIdx)<< std::endl;
                }
            else
                {
                Z1 = 1.;
                problem().variables().totalConcentration(globalIdx, nCompIdx) = 0.;
                Dune::dgrave << "Regularize totalConcentration(globalIdx, nCompIdx) = "
                    << problem().variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
                }
        }
        #endif

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
        // initialize saturation, capillary pressure
        problem().variables().saturation(globalIdx) = fluidState.saturation(wPhaseIdx);
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
        Scalar massn = problem().variables().numericalDensity(globalIdx, nPhaseIdx) = sumConc * fluidState.phaseMassFraction(nPhaseIdx);

        if ((problem().variables().densityWetting(globalIdx)*problem().variables().densityNonwetting(globalIdx)) == 0)
            DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate: try to divide by 0 density");
        Scalar vol = massw / problem().variables().densityWetting(globalIdx) + massn / problem().variables().densityNonwetting(globalIdx);
        if (dt != 0)
        {
            problem().variables().volErr()[globalIdx] = (vol - problem().spatialParameters().porosity(globalPos, *eIt));

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
    }
    return;
}

//! partial derivatives of the volumes w.r.t. changes in total concentration and pressure
/*!
 * This method calculates the volume derivatives via a secant method, where the
 * secants are gained in a pre-computational step via the transport equation and
 * the last TS size.
 * The partial derivatives w.r.t. mass are defined as
 * \f$ \frac{\partial v}{\partial C^{\kappa}} = \frac{\partial V}{\partial m^{\kappa}}\f$
 *
 * \param globalPos The global position of the current element
 * \param ep A pointer to the current element
 * \param[out] dv_dC1 partial derivative of fluid volume w.r.t. mass of component 1 [m^3/kg]
 * \param[out] dv_dC2 partial derivative of fluid volume w.r.t. mass of component 2 [m^3/kg]
 * \param[out] dv_dp partial derivative of fluid volume w.r.t. pressure [1/Pa]
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::volumeDerivatives(GlobalPosition globalPos, ElementPointer ep, Scalar& dv_dC1, Scalar& dv_dC2, Scalar& dv_dp)
{
    // cell index
    int globalIdx = problem().variables().index(*ep);

    // get cell temperature
    Scalar temperature_ = problem().temperatureAtPos(globalPos);

    // initialize an Fluid state for the update
    FluidState updFluidState;

    /**********************************
     * a) get necessary variables
     **********************************/
    //determine phase pressures from primary pressure variable
    PhaseVector pressure(0.);
    switch (pressureType)
    {
        case pw:
        {
            pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx];
            pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx] + problem().variables().capillaryPressure(globalIdx);
            break;
        }
        case pn:
        {
            pressure[wPhaseIdx] = problem().variables().pressure()[globalIdx] - problem().variables().capillaryPressure(globalIdx);
            pressure[nPhaseIdx] = problem().variables().pressure()[globalIdx];
            break;
        }
    }

    Scalar v_w = 1. / problem().variables().densityWetting(globalIdx);
    Scalar v_g = 1. / problem().variables().densityNonwetting(globalIdx);
    Scalar sati = problem().variables().saturation()[globalIdx];
    // mass of components inside the cell
    Scalar m1 = problem().variables().totalConcentration(globalIdx, wCompIdx);
    Scalar m2 = problem().variables().totalConcentration(globalIdx, nCompIdx);
    // mass fraction of wetting phase
    Scalar nuw1 =  sati/ v_w / (sati/v_w + (1-sati)/v_g);
    // actual fluid volume
    Scalar volalt = (m1+m2) * (nuw1 * v_w + (1-nuw1) * v_g);

    /**********************************
     * b) define increments
     **********************************/
    // increments for numerical derivatives
    Scalar inc1 = (fabs(problem().variables().updateEstimate(globalIdx, wCompIdx)) > 1e-8 / v_w) ?  problem().variables().updateEstimate(globalIdx,wCompIdx) : 1e-8/v_w;
    Scalar inc2 =(fabs(problem().variables().updateEstimate(globalIdx, nCompIdx)) > 1e-8 / v_g) ?  problem().variables().updateEstimate(globalIdx,nCompIdx) : 1e-8 / v_g;
    Scalar incp = 1e-2;


    /**********************************
     * c) Secant method for derivatives
     **********************************/

    // numerical derivative of fluid volume with respect to pressure
    PhaseVector p_(incp);
    p_ += pressure;
    Scalar Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1,
            p_, problem().spatialParameters().porosity(globalPos, *ep), temperature_);
    Scalar v_w_ = 1. / FluidSystem::density(updFluidState, wPhaseIdx);
    Scalar v_g_ = 1. / FluidSystem::density(updFluidState, nPhaseIdx);
    dv_dp = ((m1+m2) * (nuw1 * v_w_ + (1-nuw1) * v_g_) - volalt) /incp;
    if (dv_dp>0)
    {
        // dV_dp > 0 is unphysical: Try inverse increment for secant
        Dune::dinfo << "dv_dp larger 0 at Idx " << globalIdx << " , try and invert secant"<< std::endl;

        p_ -= 2*incp;
        updFluidState.update(Z1,
                    p_, problem().spatialParameters().porosity(globalPos, *ep), temperature_);
        Scalar v_w_ = 1. / FluidSystem::density(updFluidState, wPhaseIdx);
        Scalar v_g_ = 1. / FluidSystem::density(updFluidState, nPhaseIdx);
        dv_dp = ((m1+m2) * (nuw1 * v_w_ + (1-nuw1) * v_g_) - volalt) /-incp;
        // dV_dp > 0 is unphysical: Try inverse increment for secant
        if (dv_dp>0)
        {
            Dune::dinfo << "dv_dp still larger after inverting secant"<< std::endl;
        }
    }

    // numerical derivative of fluid volume with respect to mass of component 1
    m1 +=  inc1;
    Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, *ep), temperature_);
    Scalar satt = updFluidState.saturation(wPhaseIdx);
    Scalar nuw = satt / v_w / (satt/v_w + (1-satt)/v_g);
    dv_dC1 = ((m1+m2) * (nuw * v_w + (1-nuw) * v_g) - volalt) /inc1;
    m1 -= inc1;

    // numerical derivative of fluid volume with respect to mass of component 2
    m2 += inc2;
    Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1, pressure, problem().spatialParameters().porosity(globalPos, *ep), temperature_);
    satt = updFluidState.saturation(wPhaseIdx);
    nuw = satt / v_w / (satt/v_w + (1-satt)/v_g);
    dv_dC2 = ((m1+m2) * (nuw * v_w + (1-nuw) * v_g) - volalt)/ inc2;
    m2 -= inc2;

    //check routines if derivatives are meaningful
    if (isnan(dv_dC1) || isinf(dv_dC1) || isnan(dv_dC2) || isinf(dv_dC2)|| isnan(dv_dp) || isinf(dv_dp))
    {
        DUNE_THROW(Dune::MathError, "NAN/inf of dV_dm. If that happens in first timestep, try smaller firstDt!");
    }
}
#endif

}//end namespace Dumux
#endif
