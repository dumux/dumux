// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * @brief  Finite Volume 2p2c Diffusion Model
 * @author Benjamin Faigle, Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{
//! The finite volume model for the solution of the compositional pressure equation
/*! \ingroup multiphase
 *  Provides a Finite Volume implementation for the pressure equation of a compressible
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
 * The pressure base class FVPressure assembles the matrix and right-hand-side vector and solves for the pressure vector,
 * whereas this class provides the actual entries for the matrix and RHS vector.
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2P2C
: public FVPressureCompositional<TypeTag>
{
    typedef FVPressure<TypeTag> BaseType;

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
    //! Indices of matrix and rhs entries
    /** During the assembling of the global system of equations get-functions are called
     * (getSource(), getFlux(), etc.), which return global matrix or right hand side entries
     * in a vector. These can be accessed using following indices:
    */
    enum
    {
        rhs = 1,//!<index for the right hand side entry
        matrix = 0//!<index for the global matrix entry

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

        enableVolumeIntegral = GET_PARAM(TypeTag,bool, EnableVolumeIntegral);

        maxError_=0.;
        if (pressureType != pw && pressureType != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
    }

protected:
    Problem& problem_;
    Scalar maxError_;
    bool enableVolumeIntegral;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening
    static constexpr int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
};


//! function which assembles the system of equations to be solved
/** for first == true, a source is implemented as in FVPressure2P.
 * for first == false, the source is translated into a volumentric source term:
 * \f[ V_i \sum_{\kappa} \frac{\partial v_{t}}{\partial C^{\kappa}} q^{\kappa}_i  \f].
 * \param sourceEntry The Matrix and RHS entries
 * \param elementI The element I
 * \param cellDataI Data of cell I
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::getSource(Dune::FieldVector<Scalar, 2>& sourceEntry,
											const Element& elementI,
											const CellData& cellDataI,
											const bool first)
{
    sourceEntry=0.;
    // cell volume & perimeter, assume linear map here
    Scalar volume = elementI.geometry().volume();

    // get sources from problem
    PrimaryVariables source(NAN);
    problem().source(source, elementI);

    if(first)
    {
        source[contiWEqIdx] /= cellDataI.density(wPhaseIdx);
        source[contiNEqIdx] /= cellDataI.density(nPhaseIdx);
    }
    else
    {
        // get global coordinate of cell center
        const GlobalPosition& globalPos = elementI.geometry().center();
        // derivatives of the fluid volume with respect to concentration of components, or pressure
        if (!cellDataI.hasVolumeDerivatives())
            this->volumeDerivatives(globalPos, elementI);

        source[contiWEqIdx] *= cellDataI.dv(wCompIdx);        // dV_[i][1] = dv_dC1 = dV/dm1
        source[contiNEqIdx] *= cellDataI.dv(nCompIdx);
    }
    sourceEntry[rhs] = volume * (source[contiWEqIdx] + source[contiNEqIdx]);

    return;
}
//! function which assembles the system of equations to be solved
/** for first == true, there is no storage contribution.
 * for first == false, the storage term comprises the compressibility (due to a change in
 * pressure from last timestep):
 *  \f[ V_i c_{t,i} \frac{p^t_i - p^{t-\Delta t}_i}{\Delta t} \f]
 * and the damped error introduced by the incorrect transport of the last timestep:
 *  \f[ V_i \alpha_r \frac{v_{t} - \phi}{\Delta t} \f].
 * The latter is damped according to Fritz 2011.
 * \param storageEntry The Matrix and RHS entries
 * \param elementI The element I
 * \param cellDataI Data of cell I
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::getStorage(Dune::FieldVector<Scalar, 2>& storageEntry,
											const Element& elementI,
											const CellData& cellDataI,
                                            const bool first)
{
    storageEntry = 0.;
    // cell index
    int globalIdxI = problem().variables().index(elementI);
    Scalar volume = elementI.geometry().volume();

    // determine maximum error to scale error-term
    Scalar timestep_ = problem().timeManager().timeStepSize();

    // compressibility term
    if (!first && timestep_ != 0)
    {
        Scalar compress_term = cellDataI.dv_dp() / timestep_;

        storageEntry[matrix] -= compress_term*volume;
        // cellData has data from last TS, and pressurType points to
        // the pressure Index used as a Primary Variable
        storageEntry[rhs] -= cellDataI.pressure(pressureType) * compress_term * volume;

        if (isnan(compress_term) || isinf(compress_term))
            DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << globalIdxI);

        if(!GET_PROP_VALUE(TypeTag, EnableCompressibility))
            DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
    }

    // Abort error damping if there will be a possibly tiny timestep compared with last one
    // This might be the case if the episode or simulation comes to an end.
    if( problem().timeManager().episodeWillBeOver()
            || problem().timeManager().willBeFinished())
    {
        problem().variables().cellData(globalIdxI).errorCorrection() = 0.;
        return;
    }

    // error reduction routine: volumetric error is damped and inserted to right hand side
    // if damping is not done, the solution method gets unstable!
    problem().variables().cellData(globalIdxI).volumeError() /= timestep_;
    Scalar maxError = maxError_;
    Scalar erri = fabs(cellDataI.volumeError());
    Scalar x_lo = ErrorTermLowerBound_;
    Scalar x_mi = ErrorTermUpperBound_;
    Scalar fac  = ErrorTermFactor_;
    if (pressureType == pw)
        fac = 0.1*ErrorTermFactor_;
    Scalar lofac = 0.;
    Scalar hifac = 0.;

    if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxError))
    {
        if (erri <= x_mi * maxError)
            storageEntry[rhs] +=
                    problem().variables().cellData(globalIdxI).errorCorrection() =
                            fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/maxError)
                                * cellDataI.volumeError() * volume;
        else
            storageEntry[rhs] +=
                    problem().variables().cellData(globalIdxI).errorCorrection() =
                            fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxError)
                                * cellDataI.volumeError() * volume;
    }

    return;
}


//! Get flux at an interface between two cells
/** for first == true, the flux is calculated in traditional fractional-flow forn as in FVPressure2P.
 * for first == false, the flux thorugh \f$ \gamma{ij} \f$  is calculated via a volume balance formulation
 *  \f[ - A_{\gamma_{ij}} \cdot \mathbf{K} \cdot \mathbf{u} \cdot (\mathbf{n}_{\gamma_{ij}} \cdot \mathbf{u})
      \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa} \frac{\partial v_{t}}{\partial C^{\kappa}} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}\right)
    + V_i \frac{A_{\gamma_{ij}}}{U_i} \mathbf{K}
    \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa}
\frac{\frac{\partial v_{t,j}}{\partial C^{\kappa}_j}-\frac{\partial v_{t,i}}{\partial C^{\kappa}_i}}{\Delta x} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}\right) \f]
 * This includes a boundary integral and a volume integral, because
 *  \f$ \frac{\partial v_{t,i}}{\partial C^{\kappa}_i} \f$ is not constant.
 * Here, \f$ \mathbf{u} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma_{ij}} \f$
 * represents the normal of the face \f$ \gamma{ij} \f$.
 * \param entries The Matrix and RHS entries
 * \param intersection Intersection between cell I and J
 * \param cellDataI Data of cell I
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI, const bool first)
{
    entries = 0.;
    ElementPointer elementPointerI = intersection.inside();

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementPointerI->geometry().volume();
    Scalar perimeter = cellDataI.perimeter();

    const GlobalPosition& gravity_ = problem().gravity();

    // get absolute permeability
    FieldMatrix permeabilityI(problem().spatialParameters().intrinsicPermeability(*elementPointerI));

    // get mobilities and fractional flow factors
    Scalar fractionalWI=0, fractionalNWI=0;
    if (first)
    {
        fractionalWI = cellDataI.mobility(wPhaseIdx)
                / (cellDataI.mobility(wPhaseIdx)+ cellDataI.mobility(nPhaseIdx));
        fractionalNWI = cellDataI.mobility(nPhaseIdx)
                / (cellDataI.mobility(wPhaseIdx)+ cellDataI.mobility(nPhaseIdx));
    }

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
        = problem().spatialParameters().intrinsicPermeability(*neighborPointer);

    // compute vectorized permeabilities
    FieldMatrix meanPermeability(0);
    Dumux::harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

    Dune::FieldVector<Scalar, dim> permeability(0);
    meanPermeability.mv(unitDistVec, permeability);

    // get average density for gravity flux
    Scalar rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
    Scalar rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));

    // reset potential gradients
    Scalar potentialW = 0;
    Scalar potentialNW = 0;

    if (first)     // if we are at the very first iteration we can't calculate phase potentials
    {
        // get fractional flow factors in neigbor
        Scalar fractionalWJ = cellDataJ.mobility(wPhaseIdx)
                / (cellDataJ.mobility(wPhaseIdx)+ cellDataJ.mobility(nPhaseIdx));
        Scalar fractionalNWJ = cellDataJ.mobility(nPhaseIdx)
                / (cellDataJ.mobility(wPhaseIdx)+ cellDataJ.mobility(nPhaseIdx));

        // perform central weighting
        Scalar lambda = (cellDataI.mobility(wPhaseIdx) + cellDataJ.mobility(wPhaseIdx)) * 0.5
                + (cellDataI.mobility(nPhaseIdx) + cellDataJ.mobility(nPhaseIdx)) * 0.5;

        entries[0] = fabs(lambda*faceArea*(permeability*unitOuterNormal)/(dist));

        Scalar factor = (fractionalWI + fractionalWJ) * (rhoMeanW) * 0.5
                    + (fractionalNWI + fractionalNWJ) * (rhoMeanNW) * 0.5;
        entries[1] = factor * lambda * faceArea * (permeability * gravity_)
                * (unitOuterNormal * unitDistVec);
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
        Scalar densityW = rhoMeanW;
        Scalar densityNW = rhoMeanNW;

        potentialW = (unitOuterNormal * unitDistVec) *(cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx))/dist;
        potentialNW = (unitOuterNormal * unitDistVec) *(cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx))/dist;

        potentialW += densityW * (unitDistVec * gravity_);
        potentialNW += densityNW * (unitDistVec * gravity_);

        // initialize convenience shortcuts
        Scalar lambdaW(0.), lambdaN(0.);
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
            dV_w *= cellDataI.density(wPhaseIdx);
            gV_w *= cellDataI.density(wPhaseIdx);
        }
        else
        {
            dV_w = (dv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
            lambdaW = cellDataJ.mobility(wPhaseIdx);
            gV_w = (graddv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
            dV_w *= cellDataJ.density(wPhaseIdx);
            gV_w *= cellDataJ.density(wPhaseIdx);
        }
        if (potentialNW >= 0.)
        {
            dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
            lambdaN = cellDataI.mobility(nPhaseIdx);
            gV_n = (graddv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
            dV_n *= cellDataI.density(nPhaseIdx);
            gV_n *= cellDataI.density(nPhaseIdx);
        }
        else
        {
            dV_n = (dv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
            lambdaN = cellDataJ.mobility(nPhaseIdx);
            gV_n = (graddv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
            dV_n *= cellDataJ.density(nPhaseIdx);
            gV_n *= cellDataJ.density(nPhaseIdx);
        }

        //calculate current matrix entry
        entries[0] = faceArea * (lambdaW * dV_w + lambdaN * dV_n);
        if(enableVolumeIntegral)
            entries[0] -= volume * faceArea / perimeter * (lambdaW * gV_w + lambdaN * gV_n);     // = boundary integral - area integral
        entries[0] *= fabs((permeability*unitOuterNormal)/(dist));

        //calculate right hand side
        entries[rhs] = faceArea  * (unitOuterNormal * unitDistVec) * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);
        if(enableVolumeIntegral)
            entries[rhs] -= volume * faceArea / perimeter * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n);
        entries[rhs] *= (permeability * gravity_);         // = multipaper eq(3.3) line 2+3

        // calculate capillary pressure gradient
        Dune::FieldVector<Scalar, dim> pCGradient = unitDistVec;
        pCGradient *= (cellDataI.capillaryPressure() - cellDataJ.capillaryPressure()) / dist;

        // include capillary pressure fluxes
        switch (pressureType)
        {
        case pw:
        {
            //add capillary pressure term to right hand side
            entries[rhs] += lambdaN * dV_n * (permeability * pCGradient) * faceArea;
            if(enableVolumeIntegral)
                entries[rhs]-= lambdaN * gV_n * (permeability * pCGradient) * volume * faceArea / perimeter;
            break;
        }
        case pn:
        {
            //add capillary pressure term to right hand side
            entries[rhs] -= lambdaW * dV_w * (permeability * pCGradient) * faceArea;
            if(enableVolumeIntegral)
                entries[rhs]+= lambdaW * gV_w * (permeability * pCGradient) * volume * faceArea / perimeter;
            break;
        }
        }
    }   // end !first
}

//! Get flux on Boundary
/** for first == true, the flux is calculated in traditional fractional-flow forn as in FVPressure2P.
 * for first == false, the flux thorugh \f$ \gamma{ij} \f$  is calculated via a volume balance formulation
 *  \f[ - A_{\gamma_{ij}} \cdot \mathbf{K} \cdot \mathbf{u} \cdot (\mathbf{n}_{\gamma_{ij}} \cdot \mathbf{u})
      \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa} \frac{\partial v_{t}}{\partial C^{\kappa}} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}\right) \;, \f]
 * where we skip the volume integral assuming  \f$ \frac{\partial v_{t,i}}{\partial C^{\kappa}_i} \f$
 * to be constant at the boundary.
 * Here, \f$ \mathbf{u} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma_{ij}} \f$
 * represents the normal of the face \f$ \gamma{ij} \f$.
 * \param entries The Matrix and RHS entries
 * \param intersection Intersection between cell I and J
 * \param cellDataI Data of cell I
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI, const bool first)
{
    entries = 0.;
    // get global coordinate of cell center
    ElementPointer elementPointerI = intersection.inside();
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // get normal vector
    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();
    // get face volume
    Scalar faceArea = intersection.geometry().volume();

    // get volume derivatives inside the cell
    Scalar dv_dC1 = cellDataI.dv(wCompIdx);
    Scalar dv_dC2 = cellDataI.dv(nCompIdx);

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
        FieldMatrix permeabilityI(problem().spatialParameters().intrinsicPermeability(*elementPointerI));
        const GlobalPosition& gravity_ = problem().gravity();

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
            fractionalWI = cellDataI.mobility(wPhaseIdx)
                    / (cellDataI.mobility(wPhaseIdx)+ cellDataI.mobility(nPhaseIdx));
            fractionalNWI = cellDataI.mobility(nPhaseIdx)
                    / (cellDataI.mobility(wPhaseIdx)+ cellDataI.mobility(nPhaseIdx));

            Scalar lambda = cellDataI.mobility(wPhaseIdx)+cellDataI.mobility(nPhaseIdx);
            entries[matrix] += lambda * faceArea * (permeability * unitOuterNormal) / (dist);
            Scalar pressBC = primaryVariablesOnBoundary[Indices::pressureEqIdx];
            entries[rhs] += lambda * faceArea * pressBC * (permeability * unitOuterNormal) / (dist);
            Scalar rightentry = (fractionalWI * cellDataI.density(wPhaseIdx)
                                 + fractionalNWI * cellDataI.density(nPhaseIdx))
                                 * lambda * faceArea * (unitOuterNormal * unitDistVec) *(permeability*gravity_);
            entries[rhs] -= rightentry;
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
            if(GET_PROP_VALUE(TypeTag, BoundaryMobility) == Indices::satDependent)
            {
                lambdaWBound = BCfluidState.saturation(wPhaseIdx)
                        / viscosityWBound;
                lambdaNWBound = BCfluidState.saturation(nPhaseIdx)
                        / viscosityNWBound;
            }
            else if(GET_PROP_VALUE(TypeTag, BoundaryMobility) == Indices::permDependent)
            {
                lambdaWBound
                    = MaterialLaw::krw(problem().spatialParameters().materialLawParams(*elementPointerI),
                            BCfluidState.saturation(wPhaseIdx)) / viscosityWBound;
                lambdaNWBound
                    = MaterialLaw::krn(problem().spatialParameters().materialLawParams(*elementPointerI),
                            BCfluidState.saturation(wPhaseIdx)) / viscosityNWBound;
            }
            // get average density
            Scalar rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + densityWBound);
            Scalar rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + densityNWBound);

            Scalar potentialW = 0;
            Scalar potentialNW = 0;
            if (!first)
            {
//                            potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                            potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
//
//                            // do potential upwinding according to last potGradient vs Jochen: central weighting
//                            densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : densityWBound;
//                            densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : densityNWBound;
//
//                            densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                            densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                Scalar densityW=rhoMeanW;
                Scalar densityNW=rhoMeanNW;

                //calculate potential gradient
                potentialW = (cellDataI.pressure(wPhaseIdx) - pressBC[wPhaseIdx])/dist;
                potentialNW = (cellDataI.pressure(nPhaseIdx) - pressBC[nPhaseIdx])/dist;

                potentialW += densityW * (unitDistVec * gravity_);
                potentialNW += densityNW * (unitDistVec * gravity_);
           }   //end !first

            //do the upwinding of the mobility depending on the phase potentials
            Scalar lambdaW, lambdaNW;
            Scalar densityW(0.), densityNW(0.);
            Scalar dV_w, dV_n;     // no gV, because dv/dC assumed to be constant at the boundary
                                   // => no area integral, only boundary integral

            if (potentialW >= 0.)
            {
                densityW = (potentialW == 0) ? rhoMeanW : cellDataI.density(wPhaseIdx);
                dV_w = (dv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                                   + dv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
                dV_w *= densityW;
                lambdaW = (potentialW == 0) ? 0.5 * (cellDataI.mobility(wPhaseIdx) + lambdaWBound)
                                            : cellDataI.mobility(wPhaseIdx);
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
                densityNW = (potentialNW == 0) ? rhoMeanNW : cellDataI.density(nPhaseIdx);
                dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
                dV_n *= densityNW;
                lambdaNW = (potentialNW == 0) ? 0.5 * (cellDataI.mobility(nPhaseIdx) + lambdaNWBound)
                                              : cellDataI.mobility(nPhaseIdx);
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
                    * faceArea ; // * (unitOuterNormal * unitDistVec) is done after pc gradient inclusion

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
            entries[matrix] += entry;
            entries[rhs] += entry * primaryVariablesOnBoundary[Indices::pressureEqIdx];
            entries[rhs] -= rightEntry * (unitOuterNormal * unitDistVec);
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
            J[contiWEqIdx] /= cellDataI.density(wPhaseIdx);
            J[contiNEqIdx] /= cellDataI.density(nPhaseIdx);
        }
        else
        {
            J[contiWEqIdx] *= dv_dC1;
            J[contiNEqIdx] *= dv_dC2;
        }

        entries[rhs] -= (J[contiWEqIdx] + J[contiNEqIdx]) * faceArea;
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Boundary Condition neither Dirichlet nor Neumann!");

    return;
}

//! updates secondary variables
/*
 * A loop through all elements updates the secondary variables stored in the variableclass
 * by using the updated primary variables.
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLaws()
{
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
//! updates secondary variables of one cell
/* For each element, the secondary variables are updated according to the
 * primary variables.
 * \param elementI The element
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLawsInElement(const Element& elementI)
{
    // get global coordinate of cell center
    GlobalPosition globalPos = elementI.geometry().center();

    // cell Index and cell data
    int globalIdx = problem().variables().index(elementI);
    CellData& cellData = problem().variables().cellData(globalIdx);

    // acess the fluid state and prepare for manipulation
    FluidState& fluidState = cellData.manipulateFluidState();

    Scalar temperature_ = problem().temperatureAtPos(globalPos);
    // reset
    cellData.reset();

//    // make shure total concentrations from solution vector are exact in fluidstate
//    fluidState.setMassConcentration(wCompIdx,
//            problem().transportModel().totalConcentration(wCompIdx,globalIdx));
//    fluidState.setMassConcentration(nCompIdx,
//            problem().transportModel().totalConcentration(nCompIdx,globalIdx));
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
            problem().transportModel().totalConcentration(wCompIdx, globalIdx) = 0.;
            Dune::dgrave << "Regularize totalConcentration(wCompIdx) = "
                << cellData.totalConcentration(wCompIdx)<< std::endl;
            }
        else
            {
            Z1 = 1.;
            cellData.setTotalConcentration(nCompIdx, 0.);
            problem().transportModel().totalConcentration(nCompIdx,globalIdx) = 0.;
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
    fluidState.update(Z1, pressure, problem().spatialParameters().porosity(elementI), temperature_);

    // iterations part in case of enabled capillary pressure
    Scalar pc(0.), oldPc(0.);
    if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
    {
        pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(elementI),
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
            fluidState.update(Z1, pressure, problem().spatialParameters().porosity(elementI),
                                problem().temperatureAtPos(globalPos));
            pc = MaterialLaw::pC(problem().spatialParameters().materialLawParams(elementI),
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
    // initialize phase properties not stored in fluidstate
    cellData.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState, wPhaseIdx));
    cellData.setViscosity(nPhaseIdx, FluidSystem::viscosity(fluidState, nPhaseIdx));

    // initialize mobilities
    cellData.setMobility(wPhaseIdx, MaterialLaw::krw(problem().spatialParameters().materialLawParams(elementI), fluidState.saturation(wPhaseIdx))
                / cellData.viscosity(wPhaseIdx));
    cellData.setMobility(nPhaseIdx, MaterialLaw::krn(problem().spatialParameters().materialLawParams(elementI), fluidState.saturation(wPhaseIdx))
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
        cellData.volumeError()=(vol - problem().spatialParameters().porosity(elementI));

        if (std::isnan(cellData.volumeError()))
        {
            DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate:\n"
                    << "volErr[" << globalIdx << "] isnan: vol = " << vol
                    << ", massw = " << massw << ", rho_l = " << cellData.density(wPhaseIdx)
                    << ", massn = " << massn << ", rho_g = " << cellData.density(nPhaseIdx)
                    << ", poro = " << problem().spatialParameters().porosity(elementI)
                    << ", dt = " << problem().timeManager().timeStepSize());
        }
    }
    else
        cellData.volumeError()=0.;

    cellData.reset();

    return;
}

}//end namespace Dumux
#endif
