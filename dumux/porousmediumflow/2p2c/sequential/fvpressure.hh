// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Finite volume 2p2c pressure model
 */
#ifndef DUMUX_FVPRESSURE2P2C_HH
#define DUMUX_FVPRESSURE2P2C_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/porousmediumflow/2p2c/sequential/fvpressurecompositional.hh>
#include <dumux/material/constraintsolvers/compositionalflash.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief The finite volume model for the solution of the compositional pressure equation.
 *
 * Provides a Finite Volume implementation for the pressure equation of a compressible
 * system with two components. An IMPES-like method is used for the sequential
 * solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 * Isothermal conditions and local thermodynamic
 * equilibrium are assumed.  Gravity is included.
 * \f[
         c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}}
         \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
          = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 *  \f]
 *  where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields
 *  \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 *  \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability,
 *  \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and
 *  \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 *
 * The pressure base class FVPressure assembles the matrix and right-hand-side vector and solves for the pressure vector,
 * whereas this class provides the actual entries for the matrix and RHS vector.
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 */
template<class TypeTag> class FVPressure2P2C
: public FVPressureCompositional<TypeTag>
{
    //the model implementation
    using Implementation = GetPropType<TypeTag, Properties::PressureModel>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureN,
        pGlobal = Indices::pressureGlobal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };

    /*!
     * \brief Indices of matrix and rhs entries
     * During the assembling of the global system of equations get-functions are called
     * (getSource(), getFlux(), etc.), which return global matrix or right hand side entries
     * in a vector. These can be accessed using following indices:
     */
    enum
    {
        rhs = 1,//!<index for the right hand side entry
        matrix = 0//!<index for the global matrix entry

    };

    // using declarations to abbreviate several dune classes...
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    // convenience shortcuts for Vectors/Matrices
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using PhaseVector = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumPhases>()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    // the typenames used for the stiffness matrix and solution vector
    using Matrix = GetPropType<TypeTag, Properties::PressureCoefficientMatrix>;

protected:
    //! @copydoc FVPressure::EntryType
    using EntryType = Dune::FieldVector<Scalar, 2>;

    Problem& problem()
    {
        return problem_;
    }
    const Problem& problem() const
    {
        return problem_;
    }

public:

    void getSource(EntryType& sourceEntry, const Element& elementI, const CellData& cellDataI, const bool first);

    void getStorage(EntryType& storageEntry, const Element& elementI, const CellData& cellDataI, const bool first);

    void getFlux(EntryType& entries, const Intersection& intersection, const CellData& cellDataI, const bool first);

    void getFluxOnBoundary(EntryType& entries, const Intersection& intersection, const CellData& cellDataI, const bool first);

    //updates secondary variables for one cell and stores in the variables object
    void updateMaterialLawsInElement(const Element& elementI, bool postTimeStep);

    /*!
     * \brief Constructs a FVPressure2P2C object
     * \param problem a problem class object
     */
    FVPressure2P2C(Problem& problem) : FVPressureCompositional<TypeTag>(problem),
        problem_(problem)
    {
        ErrorTermFactor_ = getParam<Scalar>("Impet.ErrorTermFactor");
        ErrorTermLowerBound_ = getParam<Scalar>("Impet.ErrorTermLowerBound", 0.2);
        ErrorTermUpperBound_ = getParam<Scalar>("Impet.ErrorTermUpperBound");

        enableVolumeIntegral = getParam<bool>("Impet.EnableVolumeIntegral");
        regulateBoundaryPermeability = getPropValue<TypeTag, Properties::RegulateBoundaryPermeability>();
        if(regulateBoundaryPermeability)
        {
            minimalBoundaryPermeability = getParam<Scalar>("SpatialParams.MinBoundaryPermeability");
            Dune::dinfo << " Warning: Regulating Boundary Permeability requires correct subface indices on reference elements!"
                        << std::endl;
        }
    }

protected:
    Problem& problem_;
    bool enableVolumeIntegral; //!< Enables the volume integral of the pressure equation
    bool regulateBoundaryPermeability; //!< Enables regulation of permeability in the direction of a Dirichlet Boundary Condition
    Scalar minimalBoundaryPermeability; //!< Minimal limit for the boundary permeability
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening
    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static constexpr int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();
private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}
};

/*!
 * \brief Assembles the source term
 *
 * for first == true, a source is implemented as in FVPressure2P.
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
            asImp_().volumeDerivatives(globalPos, elementI);

        source[contiWEqIdx] *= cellDataI.dv(wCompIdx);        // dV_[i][1] = dv_dC1 = dV/dm1
        source[contiNEqIdx] *= cellDataI.dv(nCompIdx);
    }
    sourceEntry[rhs] = volume * (source[contiWEqIdx] + source[contiNEqIdx]);

    return;
}

/*!
 * \brief Assembles the storage term
 *
 * for first == true, there is no storage contribution.
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
    int eIdxGlobalI = problem().variables().index(elementI);
    Scalar volume = elementI.geometry().volume();

    // determine maximum error to scale error-term
    Scalar timestep_ = problem().timeManager().timeStepSize();

    // compressibility term
    if (!first && Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(timestep_, 0.0, 1.0e-30))
    {
        Scalar compress_term = cellDataI.dv_dp() / timestep_;

        storageEntry[matrix] -= compress_term*volume;
        // cellData has data from last TS, and pressurType points to
        // the pressure Index used as a Primary Variable
        storageEntry[rhs] -= cellDataI.pressure(pressureType) * compress_term * volume;

        using std::isnan;
        using std::isinf;
        if (isnan(compress_term) || isinf(compress_term))
            DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << eIdxGlobalI);

        if(!getPropValue<TypeTag, Properties::EnableCompressibility>())
            DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
    }

    // Abort error damping if there will be a possibly tiny timestep compared with last one
    // This might be the case if the episode or simulation comes to an end.
    if( problem().timeManager().episodeWillBeFinished()
            || problem().timeManager().willBeFinished())
    {
        problem().variables().cellData(eIdxGlobalI).errorCorrection() = 0.;
        return;
    }

    // error reduction routine: volumetric error is damped and inserted to right hand side
    // if damping is not done, the solution method gets unstable!
    problem().variables().cellData(eIdxGlobalI).volumeError() /= timestep_;
    Scalar erri = fabs(cellDataI.volumeError());
    Scalar x_lo = ErrorTermLowerBound_;
    Scalar x_mi = ErrorTermUpperBound_;
    Scalar fac  = ErrorTermFactor_;
    Scalar lofac = 0.;
    Scalar hifac = 1.-x_mi;

    if ((erri*timestep_ > 5e-5) && (erri > x_lo * this->maxError_))
    {
        if (erri <= x_mi * this->maxError_)
            storageEntry[rhs] +=
                    problem().variables().cellData(eIdxGlobalI).errorCorrection() =
                            fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/this->maxError_)
                                * cellDataI.volumeError() * volume;
        else
            storageEntry[rhs] +=
                    problem().variables().cellData(eIdxGlobalI).errorCorrection() =
                            fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/this->maxError_)
                                * cellDataI.volumeError() * volume;
    }
    else
        problem().variables().cellData(eIdxGlobalI).errorCorrection() = 0.;

    return;
}

/*!
 * \brief Get flux at an interface between two cells
 *
 * for first == true, the flux is calculated in traditional fractional-flow forn as in FVPressure2P.
 * for first == false, the flux thorugh \f$ \gamma \f$  is calculated via a volume balance formulation
 *  \f[ - A_{\gamma} \mathbf{n}^T_{\gamma} \mathbf{K}  \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
     \mathbf{d}_{ij}  \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}^T \mathbf{d}_{ij} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha} \frac{\partial v_{t}}{\partial C^{\kappa}}
    + V_i \frac{A_{\gamma}}{U_i} \mathbf{d}^T \mathbf{K} \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
     \mathbf{d}_{ij}  \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}^T \mathbf{d}_{ij} \right)
          \sum_{\kappa} X^{\kappa}_{\alpha}
          \frac{\frac{\partial v_{t,j}}{\partial C^{\kappa}_j}-\frac{\partial v_{t,i}}{\partial C^{\kappa}_i}}{\Delta x} \f]
 * This includes a boundary integral and a volume integral, because
 *  \f$ \frac{\partial v_{t,i}}{\partial C^{\kappa}_i} \f$ is not constant.
 * Here, \f$ \mathbf{d}_{ij} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma} \f$
 * represents the normal of the face \f$ \gamma \f$.
 * \param entries The Matrix and RHS entries
 * \param intersection Intersection between cell I and J
 * \param cellDataI Data of cell I
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entries,
                                      const Intersection& intersection,
                                      const CellData& cellDataI,
                                      const bool first)
{
    entries = 0.;
    auto elementI = intersection.inside();
    int eIdxGlobalI = problem().variables().index(elementI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementI.geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementI.geometry().volume();
    Scalar perimeter = cellDataI.perimeter();
//#warning perimeter hack 2D!
//    perimeter = intersection.geometry().volume()*2;

    const GlobalPosition& gravity_ = problem().gravity();

    // get absolute permeability
    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(elementI));

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
    auto neighbor = intersection.outside();
    int eIdxGlobalJ = problem().variables().index(neighbor);
    CellData& cellDataJ = problem().variables().cellData(eIdxGlobalJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    DimMatrix permeabilityJ
        = problem().spatialParams().intrinsicPermeability(neighbor);

    // compute vectorized permeabilities
    DimMatrix meanPermeability(0);
    harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

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

        entries[0] = fabs(lambda*faceArea*fabs(permeability*unitOuterNormal)/(dist));

        Scalar factor = (fractionalWI + fractionalWJ) * (rhoMeanW) * 0.5
                    + (fractionalNWI + fractionalNWJ) * (rhoMeanNW) * 0.5;
        entries[1] = factor * lambda * faceArea * fabs(unitOuterNormal*permeability)
                * (gravity_ * unitDistVec);
    }
    else
    {
        // determine volume derivatives
        if (!cellDataJ.hasVolumeDerivatives())
            asImp_().volumeDerivatives(globalPosNeighbor, neighbor);

        Scalar dv_dC1 = (cellDataJ.dv(wPhaseIdx)
                    + cellDataI.dv(wPhaseIdx)) / 2; // dV/dm1= dv/dC^1
        Scalar dv_dC2 = (cellDataJ.dv(nPhaseIdx)
                    + cellDataI.dv(nPhaseIdx)) / 2;

        Scalar graddv_dC1 = (cellDataJ.dv(wPhaseIdx)
                                - cellDataI.dv(wPhaseIdx)) / dist;
        Scalar graddv_dC2 = (cellDataJ.dv(nPhaseIdx)
                                - cellDataI.dv(nPhaseIdx)) / dist;

//                    potentialW = problem().variables().potentialWetting(eIdxGlobalI, isIndex);
//                    potentialNW = problem().variables().potentialNonwetting(eIdxGlobalI, isIndex);
//
//                    densityW = (potentialW > 0.) ? densityWI : densityWJ;
//                    densityNW = (potentialNW > 0.) ? densityNWI : densityNWJ;
//
//                    densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                    densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
        //jochen: central weighting for gravity term
        Scalar densityW = rhoMeanW;
        Scalar densityNW = rhoMeanNW;

        potentialW = (cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx))/dist;
        potentialNW = (cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx))/dist;

        potentialW += densityW * (unitDistVec * gravity_);
        potentialNW += densityNW * (unitDistVec * gravity_);

        // initialize convenience shortcuts
        Scalar lambdaW(0.), lambdaN(0.);
        Scalar dV_w(0.), dV_n(0.);        // dV_a = \sum_k \rho_a * dv/dC^k * X^k_a
        Scalar gV_w(0.), gV_n(0.);        // multipaper eq(3.3) line 3 analogon dV_w


        //do the upwinding of the mobility depending on the phase potentials
        const CellData* upwindWCellData(0);
        const CellData* upwindNCellData(0);
        if (potentialW > 0.)
            upwindWCellData = &cellDataI;
        else if (potentialW < 0.)
            upwindWCellData = &cellDataJ;
        else
        {
            if(cellDataI.isUpwindCell(intersection.indexInInside(), contiWEqIdx))
                upwindWCellData = &cellDataI;
            else if(cellDataJ.isUpwindCell(intersection.indexInOutside(), contiWEqIdx))
                upwindWCellData = &cellDataJ;
            //else
            //  upwinding is not done!
        }

        if (potentialNW > 0.)
            upwindNCellData = &cellDataI;
        else if (potentialNW < 0.)
            upwindNCellData = &cellDataJ;
        else
        {
            if(cellDataI.isUpwindCell(intersection.indexInInside(), contiNEqIdx))
                upwindNCellData = &cellDataI;
            else if(cellDataJ.isUpwindCell(intersection.indexInOutside(), contiNEqIdx))
                upwindNCellData = &cellDataJ;
            //else
            //  upwinding is not done!
        }

        //perform upwinding if desired
        if(!upwindWCellData or (cellDataI.wasRefined() && cellDataJ.wasRefined() && elementI.father() == neighbor.father()))
        {
            if (cellDataI.wasRefined() && cellDataJ.wasRefined())
            {
                problem().variables().cellData(eIdxGlobalI).setUpwindCell(intersection.indexInInside(), contiWEqIdx, false);
                cellDataJ.setUpwindCell(intersection.indexInOutside(), contiWEqIdx, false);
            }

            Scalar averagedMassFraction[2];
            averagedMassFraction[wCompIdx]
               = harmonicMean(cellDataI.massFraction(wPhaseIdx, wCompIdx), cellDataJ.massFraction(wPhaseIdx, wCompIdx));
            averagedMassFraction[nCompIdx]
               = harmonicMean(cellDataI.massFraction(wPhaseIdx, nCompIdx), cellDataJ.massFraction(wPhaseIdx, nCompIdx));
            Scalar averageDensity = harmonicMean(cellDataI.density(wPhaseIdx), cellDataJ.density(wPhaseIdx));

            //compute means
            dV_w = dv_dC1 * averagedMassFraction[wCompIdx] + dv_dC2 * averagedMassFraction[nCompIdx];
            dV_w *= averageDensity;
            gV_w = graddv_dC1 * averagedMassFraction[wCompIdx] + graddv_dC2 * averagedMassFraction[nCompIdx];
            gV_w *= averageDensity;
            lambdaW = harmonicMean(cellDataI.mobility(wPhaseIdx), cellDataJ.mobility(wPhaseIdx));
        }
        else //perform upwinding
        {
            dV_w = (dv_dC1 * upwindWCellData->massFraction(wPhaseIdx, wCompIdx)
                    + dv_dC2 * upwindWCellData->massFraction(wPhaseIdx, nCompIdx));
            lambdaW = upwindWCellData->mobility(wPhaseIdx);
            gV_w = (graddv_dC1 * upwindWCellData->massFraction(wPhaseIdx, wCompIdx)
                    + graddv_dC2 * upwindWCellData->massFraction(wPhaseIdx, nCompIdx));
            dV_w *= upwindWCellData->density(wPhaseIdx);
            gV_w *= upwindWCellData->density(wPhaseIdx);
        }

        if(!upwindNCellData or (cellDataI.wasRefined() && cellDataJ.wasRefined()))
        {
            if (cellDataI.wasRefined() && cellDataJ.wasRefined())
            {
                problem().variables().cellData(eIdxGlobalI).setUpwindCell(intersection.indexInInside(), contiNEqIdx, false);
                cellDataJ.setUpwindCell(intersection.indexInOutside(), contiNEqIdx, false);
            }
            Scalar averagedMassFraction[2];
            averagedMassFraction[wCompIdx]
               = harmonicMean(cellDataI.massFraction(nPhaseIdx, wCompIdx), cellDataJ.massFraction(nPhaseIdx, wCompIdx));
            averagedMassFraction[nCompIdx]
               = harmonicMean(cellDataI.massFraction(nPhaseIdx, nCompIdx), cellDataJ.massFraction(nPhaseIdx, nCompIdx));
            Scalar averageDensity = harmonicMean(cellDataI.density(nPhaseIdx), cellDataJ.density(nPhaseIdx));

            //compute means
            dV_n = dv_dC1 * averagedMassFraction[wCompIdx] + dv_dC2 * averagedMassFraction[nCompIdx];
            dV_n *= averageDensity;
            gV_n = graddv_dC1 * averagedMassFraction[wCompIdx] + graddv_dC2 * averagedMassFraction[nCompIdx];
            gV_n *= averageDensity;
            lambdaN = harmonicMean(cellDataI.mobility(nPhaseIdx), cellDataJ.mobility(nPhaseIdx));
        }
        else
       {
            dV_n = (dv_dC1 * upwindNCellData->massFraction(nPhaseIdx, wCompIdx)
                    + dv_dC2 * upwindNCellData->massFraction(nPhaseIdx, nCompIdx));
            lambdaN = upwindNCellData->mobility(nPhaseIdx);
            gV_n = (graddv_dC1 * upwindNCellData->massFraction(nPhaseIdx, wCompIdx)
                    + graddv_dC2 * upwindNCellData->massFraction(nPhaseIdx, nCompIdx));
            dV_n *= upwindNCellData->density(nPhaseIdx);
            gV_n *= upwindNCellData->density(nPhaseIdx);
        }

        //calculate current matrix entry
        entries[matrix] = faceArea * (lambdaW * dV_w + lambdaN * dV_n)
                * fabs((unitOuterNormal*permeability)/(dist));
        if(enableVolumeIntegral)
            entries[matrix] -= volume * faceArea / perimeter * (lambdaW * gV_w + lambdaN * gV_n)
                * ((unitDistVec*permeability)/(dist));     // = boundary integral - area integral

        //calculate right hand side
        entries[rhs] = faceArea   * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);
        entries[rhs] *= fabs(unitOuterNormal * permeability);
        if(enableVolumeIntegral)
            entries[rhs] -= volume * faceArea / perimeter * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n)
                            * (unitDistVec * permeability);
        entries[rhs] *= (gravity_ * unitDistVec);

        // calculate capillary pressure gradient
        Scalar pCGradient = (cellDataI.capillaryPressure() - cellDataJ.capillaryPressure()) / dist;

        // include capillary pressure fluxes
        switch (pressureType)
        {
        case pw:
        {
            //add capillary pressure term to right hand side
            entries[rhs] += lambdaN * dV_n * fabs(permeability * unitOuterNormal) * pCGradient * faceArea;
            if(enableVolumeIntegral)
                entries[rhs]-= lambdaN * gV_n * (permeability * unitDistVec) * pCGradient * volume * faceArea / perimeter;
            break;
        }
        case pn:
        {
            //add capillary pressure term to right hand side
            entries[rhs] -= lambdaW * dV_w * fabs(permeability * unitOuterNormal) * pCGradient * faceArea;
            if(enableVolumeIntegral)
                entries[rhs]+= lambdaW * gV_w * (permeability * unitDistVec) * pCGradient * volume * faceArea / perimeter;
            break;
        }
        }
    }   // end !first
}

/*!
 * \brief Get flux on Boundary
 *
 * for first == true, the flux is calculated in traditional fractional-flow forn as in FVPressure2P.
 * for first == false, the flux thorugh \f$ \gamma \f$  is calculated via a volume balance formulation
 *  \f[ - A_{\gamma} \mathbf{n}^T_{\gamma} \mathbf{K} \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \mathbf{d}_{ij}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}^T \mathbf{d}_{ij} \right)
    \sum_{\kappa} \frac{\partial v_{t}}{\partial C^{\kappa}} X^{\kappa}_{\alpha} \;, \f]
 * where we skip the volume integral assuming  \f$ \frac{\partial v_{t,i}}{\partial C^{\kappa}_i} \f$
 * to be constant at the boundary.
 * Here, \f$ \mathbf{d}_{ij} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma} \f$
 * represents the normal of the face \f$ \gamma \f$.
 *
 * If a Neumann BC is set, the given (mass-)flux is directly multiplied by the volume derivative and inserted.
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
    auto elementI = intersection.inside();
    const GlobalPosition& globalPos = elementI.geometry().center();

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

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

    /**********         Dirichlet Boundary        *************/
    if (bcType.isDirichlet(Indices::pressureEqIdx))
    {
        // get absolute permeability
        DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(elementI));

        if(regulateBoundaryPermeability)
        {
            int axis = intersection.indexInInside() / 2;
            if(permeabilityI[axis][axis] < minimalBoundaryPermeability)
                permeabilityI[axis][axis] = minimalBoundaryPermeability;
        }
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
            entries[matrix] += lambda * faceArea * fabs(permeability * unitOuterNormal) / (dist);
            Scalar pressBoundary = primaryVariablesOnBoundary[Indices::pressureEqIdx];
            entries[rhs] += lambda * faceArea * pressBoundary * fabs(permeability * unitOuterNormal) / (dist);
            Scalar rightentry = (fractionalWI * cellDataI.density(wPhaseIdx)
                                 + fractionalNWI * cellDataI.density(nPhaseIdx))
                                 * lambda * faceArea * fabs(unitOuterNormal * permeability)
                                 * ( unitDistVec * gravity_);
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
            if(getPropValue<TypeTag, Properties::BoundaryMobility>() == Indices::satDependent)
            {
                lambdaWBound = BCfluidState.saturation(wPhaseIdx)
                        / viscosityWBound;
                lambdaNWBound = BCfluidState.saturation(nPhaseIdx)
                        / viscosityNWBound;
            }
            else if(getPropValue<TypeTag, Properties::BoundaryMobility>() == Indices::permDependent)
            {
                lambdaWBound
                    = fluidMatrixInteraction.krw(BCfluidState.saturation(wPhaseIdx)) / viscosityWBound;
                lambdaNWBound
                    = fluidMatrixInteraction.krn(BCfluidState.saturation(wPhaseIdx)) / viscosityNWBound;
            }
            // get average density
            Scalar rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + densityWBound);
            Scalar rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + densityNWBound);

            Scalar potentialW = 0;
            Scalar potentialNW = 0;
            if (!first)
            {
//                            potentialW = problem().variables().potentialWetting(eIdxGlobalI, isIndex);
//                            potentialNW = problem().variables().potentialNonwetting(eIdxGlobalI, isIndex);
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
                densityW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialW, 0.0, 1.0e-30)) ? rhoMeanW : cellDataI.density(wPhaseIdx);
                dV_w = (dv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                                   + dv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
                dV_w *= densityW;
                lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialW, 0.0, 1.0e-30)) ? 0.5 * (cellDataI.mobility(wPhaseIdx) + lambdaWBound)
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
                densityNW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialNW, 0.0, 1.0e-30)) ? rhoMeanNW : cellDataI.density(nPhaseIdx);
                dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                        + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
                dV_n *= densityNW;
                lambdaNW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialNW, 0.0, 1.0e-30)) ? 0.5 * (cellDataI.mobility(nPhaseIdx) + lambdaNWBound)
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
            Scalar entry = (lambdaW * dV_w + lambdaNW * dV_n)
                    * (fabs(unitOuterNormal * permeability) / dist) * faceArea;

            //calculate right hand side
            Scalar rightEntry = (lambdaW * densityW * dV_w + lambdaNW * densityNW * dV_n)
                    * fabs(unitOuterNormal * permeability) * (gravity_ * unitDistVec) * faceArea ;

            // include capillary pressure fluxes
            // calculate capillary pressure gradient
            Scalar pCGradient = (cellDataI.capillaryPressure() - pcBound) / dist;
            switch (pressureType)
            {
            case pw:
                {
                    //add capillary pressure term to right hand side
                    rightEntry += lambdaNW * dV_n * pCGradient * fabs(unitOuterNormal * permeability) * faceArea;
                    break;
                }
            case pn:
                {
                    //add capillary pressure term to right hand side
                    rightEntry -= lambdaW * dV_w * pCGradient * fabs(unitOuterNormal * permeability) * faceArea;
                    break;
                }
            }


            // set diagonal entry and right hand side entry
            entries[matrix] += entry;
            entries[rhs] += entry * primaryVariablesOnBoundary[Indices::pressureEqIdx];
            entries[rhs] -= rightEntry;
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

/*!
 * \brief Updates secondary variables of one cell
 *
 * For each element, the secondary variables are updated according to the
 * primary variables. In case the method is called after the Transport,
 * i.e. at the end / post time step, CellData2p2c.reset() resets the volume
 * derivatives for the next time step.
 * \param element The element
 * \param postTimeStep Flag indicating if we have just completed a time step
 */
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLawsInElement(const Element& element, bool postTimeStep)
{
    // get global coordinate of cell center
    GlobalPosition globalPos = element.geometry().center();

    // cell Index and cell data
    int eIdxGlobal = problem().variables().index(element);
    CellData& cellData = problem().variables().cellData(eIdxGlobal);

    // acess the fluid state and prepare for manipulation
    FluidState& fluidState = cellData.manipulateFluidState();

    Scalar temperature_ = problem().temperatureAtPos(globalPos);

    // reset to calculate new timeStep if we are at the end of a time step
    if(postTimeStep)
        cellData.reset();

    // get the overall mass of first component Z0 = C^0 / (C^0+C^1) [-]
    Scalar Z0 = cellData.massConcentration(wCompIdx)
            / (cellData.massConcentration(wCompIdx)
                    + cellData.massConcentration(nCompIdx));

    // make sure only physical quantities enter flash calculation
    if(Z0 < 0. || Z0 > 1.)
    {
        Dune::dgrave << "Feed mass fraction unphysical: Z0 = " << Z0
               << " at global Idx " << eIdxGlobal
               << " , because totalConcentration(wCompIdx) = "
               << cellData.totalConcentration(wCompIdx)
               << " and totalConcentration(nCompIdx) = "
               << cellData.totalConcentration(nCompIdx)<< std::endl;
        if(Z0 < 0.)
            {
            Z0 = 0.;
            cellData.setTotalConcentration(wCompIdx, 0.);
            problem().transportModel().totalConcentration(wCompIdx, eIdxGlobal) = 0.;
            Dune::dgrave << "Regularize totalConcentration(wCompIdx) = "
                << cellData.totalConcentration(wCompIdx)<< std::endl;
            }
        else
            {
            Z0 = 1.;
            cellData.setTotalConcentration(nCompIdx, 0.);
            problem().transportModel().totalConcentration(nCompIdx,eIdxGlobal) = 0.;
            Dune::dgrave << "Regularize totalConcentration(eIdxGlobal, nCompIdx) = "
                << cellData.totalConcentration(nCompIdx)<< std::endl;
            }
    }

    PhaseVector pressure;
    CompositionalFlash<Scalar, FluidSystem> flashSolver;

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), element);


    if(getPropValue<TypeTag, Properties::EnableCapillarity>()) // iterate capillary pressure and saturation
    {
        unsigned int maxiter = 6;
        Scalar pc = cellData.capillaryPressure(); // initial guess for pc from last TS
        // start iteration loop
        for (unsigned int iter = 0; iter < maxiter; iter++)
        {
            switch (pressureType)
            {
                case pw:
                {
                    pressure[wPhaseIdx] = asImp_().pressure(eIdxGlobal);
                    pressure[nPhaseIdx] = asImp_().pressure(eIdxGlobal) + pc;
                    break;
                }
                case pn:
                {
                    pressure[wPhaseIdx] = asImp_().pressure(eIdxGlobal) - pc;
                    pressure[nPhaseIdx] = asImp_().pressure(eIdxGlobal);
                    break;
                }
            }

            // complete fluid state
            flashSolver.concentrationFlash2p2c(fluidState, Z0, pressure, temperature_);

            // calculate new pc
            Scalar oldPc = pc;
            pc = fluidMatrixInteraction.pc(fluidState.saturation(wPhaseIdx));

            if (fabs(oldPc-pc)<10 && iter != 0)
                break;

            if (iter++ == maxiter)
                Dune::dinfo << iter << "times iteration of pc was applied at Idx " << eIdxGlobal
                << ", pc delta still " << fabs(oldPc-pc) << std::endl;
        }
    }
    else  // capillary pressure neglected
    {
        pressure[wPhaseIdx] = pressure[nPhaseIdx] = asImp_().pressure()[eIdxGlobal];
        flashSolver.concentrationFlash2p2c(fluidState, Z0, pressure, temperature_);
    }

    // initialize mobilities
    cellData.setMobility(wPhaseIdx, fluidMatrixInteraction.krw(fluidState.saturation(wPhaseIdx))
                / cellData.viscosity(wPhaseIdx));
    cellData.setMobility(nPhaseIdx, fluidMatrixInteraction.krn(fluidState.saturation(wPhaseIdx))
                / cellData.viscosity(nPhaseIdx));

    // determine volume mismatch between actual fluid volume and pore volume
    Scalar sumConc = (cellData.totalConcentration(wCompIdx)
            + cellData.totalConcentration(nCompIdx));
    Scalar massw = sumConc * fluidState.phaseMassFraction(wPhaseIdx);
    Scalar massn = sumConc * fluidState.phaseMassFraction(nPhaseIdx);

    if (Dune::FloatCmp::eq<Scalar>((cellData.density(wPhaseIdx)*cellData.density(nPhaseIdx)), 0))
        DUNE_THROW(Dune::MathError, "Sequential2p2c::postProcessUpdate: try to divide by 0 density");
    Scalar vol = massw / cellData.density(wPhaseIdx) + massn / cellData.density(nPhaseIdx);
    if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(problem().timeManager().timeStepSize(), 0.0, 1.0e-30))
    {
        cellData.volumeError()=(vol - problem().spatialParams().porosity(element));

        using std::isnan;
        if (isnan(cellData.volumeError()))
        {
            DUNE_THROW(Dune::MathError, "Sequential2p2c::postProcessUpdate:\n"
                    << "volErr[" << eIdxGlobal << "] isnan: vol = " << vol
                    << ", massw = " << massw << ", rho_l = " << cellData.density(wPhaseIdx)
                    << ", massn = " << massn << ", rho_g = " << cellData.density(nPhaseIdx)
                    << ", poro = " << problem().spatialParams().porosity(element)
                    << ", dt = " << problem().timeManager().timeStepSize());
        }
    }
    else
        cellData.volumeError()=0.;

    return;
}

}//end namespace Dumux
#endif
