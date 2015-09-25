// $Id: fvpressure2p.hh 3357 2010-03-25 13:02:05Z lauser $
/*****************************************************************************
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
#ifndef DUMUX_FVPRESSURE2P2CMULTIPHYSICS_HH
#define DUMUX_FVPRESSURE2P2CMULTIPHYSICS_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include "dumux/common/pardiso.hh"
#include "dumux/common/math.hh"
#include <dumux/decoupled/2p2c/2p2cproperties.hh>
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
 *  solution of the problem.  Capillary forces and diffusion are
 *  neglected. Isothermal conditions and local thermodynamic
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
template<class TypeTag> class FVPressure2P2CMultiPhysics
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

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
        Sn = Indices::saturationNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<Scalar, 2> PhaseVector;

    // the typenames used for the stiffness matrix and solution vector
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) RHSVector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //function which assembles the system of equations to be solved
    void assemble(bool first);

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

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
    //the variables object is initialized, non-compositional before and compositional after first pressure calculation
    void initialMaterialLaws(bool compositional);

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    //initialization routine to prepare first timestep
    void initialize(bool solveTwice = false);

    //pressure solution routine: update estimate for secants, assemble, solve.
    void pressure(bool solveTwice = true)
    {
        //pre-transport to estimate update vector
        Scalar dt_estimate = 0.;
        Dune::dinfo << "secant guess"<< std::endl;
        problem_.transportModel().update(-1, dt_estimate, problem_.variables().updateEstimate(), false);
        //last argument false in update() makes shure that this is estimate and no "real" transport step
        problem_.variables().updateEstimate() *= problem_.timeManager().timeStepSize();

        assemble(false);           Dune::dinfo << "pressure calculation"<< std::endl;
        solve();

        return;
    }

    void calculateVelocity()
    {
        return;
    }

    //numerical volume derivatives wrt changes in mass, pressure
    void volumeDerivatives(GlobalPosition globalPos, ElementPointer ep, Scalar& dv_dC1, Scalar& dv_dC2, Scalar& dv_dp);

    /*! \name general methods for serialization, output */
    //@{
    // serialization methods
    template<class Restarter>
    void serialize(Restarter &res)
    {
        return;
    }

    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \brief Write data files
     /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem().variables().addOutputVtkFields(writer);

        // add multiphysics stuff
        Dune::BlockVector<Dune::FieldVector<int,1> > *subdomainPtr = writer.template createField<int, 1> (dv_dp.size());
        *subdomainPtr = problem_.variables().subdomain();
        writer.addCellData(subdomainPtr, "subdomain");

#if DUNE_MINIMAL_DEBUG_LEVEL <= 3
        // add debug stuff
        Dune::BlockVector<Dune::FieldVector<double,1> > *errorCorrPtr = writer.template createField<double, 1> (dv_dp.size());
        *errorCorrPtr = errorCorrection;
        writer.addCellData(errorCorrPtr, "Error Correction");
        // add debug stuff
        Dune::BlockVector<Dune::FieldVector<double,1> > *dv_dpPtr = writer.template createField<double, 1> (dv_dp.size());
        *dv_dpPtr = dv_dp;
        writer.addCellData(dv_dpPtr, "dv_dP");
                Dune::BlockVector<Dune::FieldVector<double,1> > *dV_dC1Ptr = writer.template createField<double, 1> (dv_dp.size());
        Dune::BlockVector<Dune::FieldVector<double,1> > *dV_dC2Ptr = writer.template createField<double, 1> (dv_dp.size());
        *dV_dC1Ptr = dV_[0];
        *dV_dC2Ptr = dV_[1];
        writer.addCellData(dV_dC1Ptr, "dV_dC1");
        writer.addCellData(dV_dC2Ptr, "dV_dC2");
#endif

        return;
    }

    //! Constructs a FVPressure2P2C object
    /**
     * \param problem a problem class object
     */
    FVPressure2P2CMultiPhysics(Problem& problem) :
        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (2 * dim + 1)
                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()),
                gravity(problem.gravity())
    {
        if (pressureType != pw && pressureType != pn && pressureType != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType != Sw && saturationType != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        //rezise block vectors
        dV_[wPhaseIdx].resize(problem_.variables().gridSize());
        dV_[nPhaseIdx].resize(problem_.variables().gridSize());
        dv_dp.resize(problem_.variables().gridSize());
        problem_.variables().subdomain().resize(problem_.variables().gridSize());
        nextSubdomain.resize(problem_.variables().gridSize());

        errorCorrection.resize(problem_.variables().gridSize());

        //prepare stiffness Matrix and Vectors
        initializeMatrix();
    }

private:
    Problem& problem_;
    Matrix A_;
    RHSVector f_;
    std::string solverName_;
    std::string preconditionerName_;

    //vectors for partial derivatives
    typename SolutionTypes::PhaseProperty dV_;
    typename SolutionTypes::ScalarSolution dv_dp;

    // subdomain map
    Dune::BlockVector<Dune::FieldVector<int,1> > nextSubdomain;  //! vector holding next subdomain

    // debug
    typename SolutionTypes::ScalarSolution errorCorrection;  //debug output

protected:
    const Dune::FieldVector<Scalar, dimWorld>& gravity; //!< vector including the gravity constant
    static const Scalar cFLFactor_ = GET_PROP_VALUE(TypeTag, PTAG(CFLFactor)); //!< determines the CFLfactor
    static const int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static const int saturationType = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$ 0 = S_w \f$, \f$ 1 = S_n \f$)
};

//! initializes the matrix to store the system of equations
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
                rowSize++;
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem_.variables().index(*outside);

                // add off diagonal index
                A_.addindex(globalIdxI, globalIdxJ);
            }
    }
    A_.endindices();

    return;
}

//! initializes the simulation run
/*!
 * Initializes the simulation to gain the initial pressure field.
 *
 * \param solveTwice flag to determine possible iterations of the initialization process
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::initialize(bool solveTwice)
{
    // assign whole domain to most complex subdomain, => 2p
    problem_.variables().subdomain() = 2;
    // initialguess: set saturations, determine visco and mobility for initial pressure equation
    // at this moment, the pressure is unknown. Hence, dont regard compositional effects.
    initialMaterialLaws(false);     Dune::dinfo << "first saturation guess"<<std::endl; //=J: initialGuess()

    assemble(true);                 Dune::dinfo << "first pressure guess"<<std::endl;
    solve();

    // update the compositional variables (hence true)
    initialMaterialLaws(true);      Dune::dinfo << "first guess for mass fractions"<<std::endl;//= J: transportInitial()

    // perform concentration update to determine secants
    Scalar dt_estimate = 0.;
    problem_.transportModel().update(0., dt_estimate, problem_.variables().updateEstimate(), false);   Dune::dinfo << "secant guess"<< std::endl;
    dt_estimate = std::min ( problem_.timeManager().timeStepSize(), dt_estimate);
    problem_.variables().updateEstimate() *= dt_estimate;

    // pressure calculation
    assemble(false);                 Dune::dinfo << "second pressure guess"<<std::endl;
    solve();
    // update the compositional variables
    initialMaterialLaws(true);

    if (solveTwice)
    {
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureOld(this->problem().variables().pressure());
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureDiff;
        Scalar pressureNorm = 1.;   //dummy initialization to perform at least 1 iteration
        int numIter = 1;

        while (pressureNorm > 1e-5 && numIter < 10)
        {
            Scalar dt_dummy=0.;    // without this dummy, it will never converge!
            // update for secants
            problem_.transportModel().update(-1, dt_dummy, problem_.variables().updateEstimate(), false);   Dune::dinfo << "secant guess"<< std::endl;
            problem_.variables().updateEstimate() *= dt_estimate;
            // pressure calculation
            assemble(false);                 Dune::dinfo << "pressure guess number "<< numIter <<std::endl;
            solve();
            // update the compositional variables
            initialMaterialLaws(true);

            pressureDiff = pressureOld;
            pressureDiff -= this->problem().variables().pressure();
            pressureNorm = pressureDiff.infinity_norm();
            pressureOld = this->problem().variables().pressure();
            pressureNorm /= pressureOld.infinity_norm();

            numIter++;
        }
    }
    return;
}

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
    A_ = 0;
    f_ = 0;

    // initialization: set the fluid volume derivatives to zero
    Scalar dv_dC1(0.), dv_dC2(0.);
    for (int i = 0; i < problem_.variables().gridSize(); i++)
    {
        dV_[wPhaseIdx][i] = 0;     // dv_dC1 = dV/dm1
        dV_[nPhaseIdx][i] = 0;     // dv / dC2
        dv_dp[i] = 0;      // dv / dp
    }

    // determine maximum error to scale error-term
    Scalar timestep_ = problem_.timeManager().timeStepSize();
    Scalar maxErr = fabs(problem_.variables().volErr().infinity_norm()) / timestep_;

    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // number of Faces of current element
        typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
        int numberOfFaces = ReferenceElements::general(gt).size(1);


        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().center();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

        /****************implement sources and sinks************************/
        Dune::FieldVector<Scalar,2> source(problem_.source(globalPos, *eIt));
        if(first || problem_.variables().subdomain(globalIdxI) == 1)
        {
                source[wPhaseIdx] /= densityWI;
                source[nPhaseIdx] /= densityNWI;
        }
        else
        {
                // derivatives of the fluid volume with respect to concentration of components, or pressure
                if (dV_[0][globalIdxI] == 0)
                    volumeDerivatives(globalPos, *eIt, dV_[wPhaseIdx][globalIdxI][0], dV_[nPhaseIdx][globalIdxI][0], dv_dp[globalIdxI][0]);

                source[wPhaseIdx] *= dV_[wPhaseIdx][globalIdxI];        // note: dV_[i][1] = dv_dC1 = dV/dm1
                source[nPhaseIdx] *= dV_[nPhaseIdx][globalIdxI];
        }
        f_[globalIdxI] = volume * (source[wPhaseIdx] + source[nPhaseIdx]);
        /********************************************************************/

        // get absolute permeability
        FieldMatrix permeabilityI(problem_.spatialParameters().intrinsicPermeability(globalPos, *eIt));

        // get mobilities and fractional flow factors
        Scalar lambdaWI = problem_.variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem_.variables().mobilityNonwetting(globalIdxI);
        Scalar fractionalWI=0, fractionalNWI=0;
        if (first)
        {
            fractionalWI = lambdaWI / (lambdaWI+ lambdaNWI);
            fractionalNWI = lambdaNWI / (lambdaWI+ lambdaNWI);
        }

        Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);

        // prepare meatrix entries
        Scalar entry=0.;
        Scalar rightEntry=0.;

        // iterate over all faces of the cell
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            typedef Dune::GenericReferenceElements<Scalar, dim - 1> FaceReferenceElements;
            const Dune::FieldVector<Scalar, dim - 1>& faceLocal =
                FaceReferenceElements::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition localPosFace(0);

            // get normal vector
            Dune::FieldVector<Scalar, dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                // gemotry info of neighbor
                Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector<Scalar, dimWorld> distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                FieldMatrix permeabilityJ
                    = problem_.spatialParameters().intrinsicPermeability(globalPosNeighbor,
                                                                            *neighborPointer);

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);
                Dumux::harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

                Dune::FieldVector<Scalar, dim> permeability(0);
                meanPermeability.mv(unitDistVec, permeability);

                // get mobilities in neighbor
                Scalar lambdaWJ = problem_.variables().mobilityWetting(globalIdxJ);
                Scalar lambdaNWJ = problem_.variables().mobilityNonwetting(globalIdxJ);

                // phase densities in cell in neighbor
                Scalar densityWJ = problem_.variables().densityWetting(globalIdxJ);
                Scalar densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);

                Scalar pcJ = problem_.variables().capillaryPressure(globalIdxJ);

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
                else if (problem_.variables().subdomain(globalIdxI) == 1) //the current cell in the 1p domain
                {
                    // 1p => no pC => only 1 pressure, potential
                    Scalar potential = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI]
                                 - problem_.variables().pressure()[globalIdxJ]) / (dist*dist);
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
                    if (problem_.variables().subdomain(globalIdxJ)== 2)
                    {
                        // determine volume derivatives
                        if (dV_[0][globalIdxJ] == 0)
                            volumeDerivatives(globalPosNeighbor, *neighborPointer, dV_[wPhaseIdx][globalIdxJ][0], dV_[nPhaseIdx][globalIdxJ][0], dv_dp[globalIdxJ][0]);
                        dv_dC1 = (dV_[wPhaseIdx][globalIdxI] + dV_[wPhaseIdx][globalIdxJ]) / 2; // dV/dm1= dV/dC^1
                        dv_dC2 = (dV_[nPhaseIdx][globalIdxI] + dV_[nPhaseIdx][globalIdxJ]) / 2;

                        graddv_dC1 = (dV_[wPhaseIdx][globalIdxJ] - dV_[wPhaseIdx][globalIdxI]) / dist;
                        graddv_dC2 = (dV_[nPhaseIdx][globalIdxJ] - dV_[nPhaseIdx][globalIdxI]) / dist;
                    }
                    else
                    {
                        dv_dC1 = dV_[wPhaseIdx][globalIdxI]; //only regard volume changes in 2p area => no averageing
                        dv_dC2 = dV_[nPhaseIdx][globalIdxI];

                    }

//                    potentialW = problem_.variables().potentialWetting(globalIdxI, isIndex);
//                    potentialNW = problem_.variables().potentialNonwetting(globalIdxI, isIndex);
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
                        potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI]
                                - problem_.variables().pressure()[globalIdxJ]) / (dist*dist);
                        potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI]
                                - problem_.variables().pressure()[globalIdxJ] + pcI - pcJ) / (dist*dist);
                        break;
                    }
                    case pn:
                    {
                        potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI]
                                - problem_.variables().pressure()[globalIdxJ] - pcI + pcJ) / (dist*dist);
                        potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI]
                                - problem_.variables().pressure()[globalIdxJ]) / (dist*dist);
                        break;
                    }
                    case pglobal:
                    {
                        DUNE_THROW(Dune::NotImplemented, "Global pressure not yet implemented for 2p2c");
//                        potentialW = (problem_.variables().pressure()[globalIdxI]
//                                - problem_.variables().pressure()[globalIdxJ] - fMeanNW * (pcI - pcJ)) / dist;
//                        potentialNW = (problem_.variables().pressure()[globalIdxI]
//                                - problem_.variables().pressure()[globalIdxJ] + fMeanW * (pcI - pcJ)) / dist;
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
                        dV_w = (dv_dC1 * problem_.variables().wet_X1(globalIdxI) + dv_dC2 * (1. - problem_.variables().wet_X1(globalIdxI)));
                        dV_w *= densityWI;
                        lambdaW = problem_.variables().mobilityWetting(globalIdxI);

                        if (problem_.variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_w = (graddv_dC1 * problem_.variables().wet_X1(globalIdxI) + graddv_dC2 * (1. - problem_.variables().wet_X1(globalIdxI)));
                            gV_w *= densityWI;
                        }
                    }
                    else
                    {
                        dV_w = (dv_dC1 * problem_.variables().wet_X1(globalIdxJ) + dv_dC2 * (1. - problem_.variables().wet_X1(globalIdxJ)));
                        dV_w *= densityWJ;
                        lambdaW = problem_.variables().mobilityWetting(globalIdxJ);
                        if (problem_.variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_w = (graddv_dC1 * problem_.variables().wet_X1(globalIdxJ) + graddv_dC2 * (1. - problem_.variables().wet_X1(globalIdxJ)));
                            gV_w *= densityWJ;
                        }
                    }
                    if (potentialNW >= 0.)
                    {
                        dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxI) + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
                        dV_n *= densityNWI;
                        lambdaN = problem_.variables().mobilityNonwetting(globalIdxI);
                        if (problem_.variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_n = (graddv_dC1 * problem_.variables().nonwet_X1(globalIdxI) + graddv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
                            gV_n *= densityNWI;
                        }
                    }
                    else
                    {
                        dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxJ) + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxJ)));
                        dV_n *= densityNWJ;
                        lambdaN = problem_.variables().mobilityNonwetting(globalIdxJ);
                        if (problem_.variables().subdomain(globalIdxJ)== 2)
                        {
                            gV_n = (graddv_dC1 * problem_.variables().nonwet_X1(globalIdxJ) + graddv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxJ)));
                            gV_n *= densityNWJ;
                        }
                    }

                    //calculate current matrix entry
                    entry = faceArea * (lambdaW * dV_w + lambdaN * dV_n);

                    //calculate right hand side
                    rightEntry = faceArea  * (unitOuterNormal * unitDistVec) * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);
                    if (problem_.variables().subdomain(globalIdxJ)!= 1) // complex 2p subdomain
                    {
                        //subtract area integral
                        entry -= volume / numberOfFaces * (lambdaW * gV_w + lambdaN * gV_n);
                        rightEntry -= volume / numberOfFaces * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n);
                    }
                    entry *= fabs((permeability*unitOuterNormal)/(dist));
                    rightEntry *= (permeability * gravity);
                }   // end !first

                //set right hand side
                f_[globalIdxI] -= rightEntry;
                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entry;
                // set off-diagonal entry
                A_[globalIdxI][globalIdxJ] = -entry;
            }   // end neighbor
            /****************************************************/


            /************* boundary face ************************/
            else
            {
                // get volume derivatives inside the cell
                dv_dC1 = dV_[wPhaseIdx][globalIdxI];
                dv_dC2 = dV_[nPhaseIdx][globalIdxI];

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().global(faceLocal);

                // geometrical information
                Dune::FieldVector<Scalar, dimWorld> distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = problem_.bcTypePress(globalPosFace, *isIt);

                /**********         Dirichlet Boundary        *************/
                if (bctype == BoundaryConditions::dirichlet)
                {
                    //permeability vector at boundary
                    Dune::FieldVector<Scalar, dim> permeability(0);
                    permeabilityI.mv(unitDistVec, permeability);

                    // create a fluid state for the boundary
                    FluidState BCfluidState;

                    Scalar temperatureBC = problem_.temperature(globalPosFace, *eIt);

                    //get dirichlet pressure boundary condition
                    Scalar pressBound = problem_.dirichletPress(globalPosFace, *isIt);

                    Scalar pressBC = 0.;
                    Scalar pcBound = 0.;
                    if(pressureType==pw)
                    {
                        pressBC = pressBound;
                    }
                    else if(pressureType==pn)
                    {
                        pcBound = problem_.variables().capillaryPressure(globalIdxI);
                        pressBC = pressBound - pcBound;
                    }

                    if (first)
                    {
                        Scalar lambda = lambdaWI+lambdaNWI;
                        A_[globalIdxI][globalIdxI] += lambda * faceArea * (permeability * unitOuterNormal) / (dist);
                        pressBC = problem_.dirichletPress(globalPosFace, *isIt);
                        f_[globalIdxI] += lambda * faceArea * pressBC * (permeability * unitOuterNormal) / (dist);
                        rightEntry = (fractionalWI * densityWI
                                             + fractionalNWI * densityNWI)
                                             * lambda * faceArea * (unitOuterNormal * unitDistVec) *(permeability*gravity);
                        f_[globalIdxI] -= rightEntry;
                    }
                    else    //not first
                    {
                        //get boundary condition type for compositional transport
                        BoundaryConditions2p2c::Flags bcform = problem_.bcFormulation(globalPosFace, *isIt);
                        if (bcform == BoundaryConditions2p2c::saturation) // saturation given
                        {
                            Scalar satBound = problem_.dirichletTransport(globalPosFace, *isIt);
                            BCfluidState.satFlash(satBound, pressBC, problem_.spatialParameters().porosity(globalPos, *eIt), temperatureBC);
                        }
                        else if (bcform == Dumux::BoundaryConditions2p2c::concentration) // mass fraction given
                        {
                            Scalar Z1Bound = problem_.dirichletTransport(globalPosFace, *isIt);
                            BCfluidState.update(Z1Bound, pressBC, problem_.spatialParameters().porosity(globalPos, *eIt), temperatureBC);
                        }
                        else    // nothing declared at boundary
                        {
                            Scalar satBound = problem_.variables().saturation()[globalIdxI];
                               BCfluidState.satFlash(satBound, pressBC, problem_.spatialParameters().porosity(globalPos, *eIt), temperatureBC);
                            Dune::dwarn << "no boundary saturation/concentration specified on boundary pos " << globalPosFace << std::endl;
                        }


                        // determine fluid properties at the boundary
                        Scalar lambdaWBound = 0.;
                        Scalar lambdaNWBound = 0.;

                        Scalar densityWBound =
                            FluidSystem::phaseDensity(wPhaseIdx, temperatureBC, pressBC, BCfluidState);
                        Scalar densityNWBound =
                            FluidSystem::phaseDensity(nPhaseIdx, temperatureBC, pressBC+pcBound, BCfluidState);
                        Scalar viscosityWBound =
                            FluidSystem::phaseViscosity(wPhaseIdx, temperatureBC, pressBC, BCfluidState);
                        Scalar viscosityNWBound =
                            FluidSystem::phaseViscosity(nPhaseIdx, temperatureBC, pressBC+pcBound, BCfluidState);

                        // mobility at the boundary
                        switch (GET_PROP_VALUE(TypeTag, PTAG(BoundaryMobility)))
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
                                    problem_.spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityWBound;
                            lambdaNWBound = MaterialLaw::krn(
                                    problem_.spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityNWBound;
                            }
                        }

                        Scalar rhoMeanW = 0.5 * (densityWI + densityWBound);
                        Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWBound);

                        Scalar potentialW = 0, potentialNW = 0;

//                            potentialW = problem_.variables().potentialWetting(globalIdxI, isIndex);
//                            potentialNW = problem_.variables().potentialNonwetting(globalIdxI, isIndex);
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
                                potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pressBound) / (dist * dist);
                                potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] + pcI - pressBound - pcBound)
                                        / (dist * dist);
                                break;
                            }
                            case pn:
                            {
                                potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pcI - pressBound + pcBound)
                                        / (dist * dist);
                                potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pressBound) / (dist * dist);
                                break;
                            }
                            case pglobal:
                            {
                                Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNWBound);
                                Scalar fractionalNWBound = lambdaNWBound / (lambdaWBound + lambdaNWBound);
                                Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
                                Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWBound);

                                potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pressBound - fMeanNW * (pcI
                                        - pcBound)) / (dist * dist);
                                potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pressBound + fMeanW * (pcI
                                        - pcBound)) / (dist * dist);
                                break;
                            }
                        }

                        potentialW += densityW * (unitDistVec * gravity);
                        potentialNW += densityNW * (unitDistVec * gravity);

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW, lambdaNW;
                        Scalar dV_w, dV_n;     // gV_a weglassen, da dV/dc am Rand ortsunabhängig angenommen -> am rand nicht bestimmbar -> nur Randintegral ohne Gebietsintegral

                        if(problem_.variables().subdomain(globalIdxI)==1)    // easy 1p subdomain
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
                                dV_w = (dv_dC1 * problem_.variables().wet_X1(globalIdxI)
                                         + dv_dC2 * (1. - problem_.variables().wet_X1(globalIdxI)));
                                dV_w *= densityW;
                                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaWI;
                            }
                            else
                            {
                                densityW = densityWBound;
                                dV_w = (dv_dC1 * BCfluidState.massFrac(wPhaseIdx, wCompIdx) + dv_dC2 * BCfluidState.massFrac(wPhaseIdx, nCompIdx));
                                dV_w *= densityW;
                                lambdaW = lambdaWBound;
                            }
                            if (potentialNW >= 0.)
                            {
                                densityNW = (potentialNW == 0) ? rhoMeanNW : densityNWI;
                                dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxI)
                                        + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
                                dV_n *= densityNW;
                                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNWI;
                            }
                            else
                            {
                                densityNW = densityNWBound;
                                dV_n = (dv_dC1 * BCfluidState.massFrac(nPhaseIdx, wCompIdx) + dv_dC2 * BCfluidState.massFrac(nPhaseIdx, nCompIdx));
                                dV_n *= densityNW;
                                lambdaNW = lambdaNWBound;
                            }

                            //calculate current matrix entry
                            entry = (lambdaW * dV_w + lambdaNW * dV_n) * ((permeability * unitDistVec) / dist) * faceArea
                                    * (unitOuterNormal * unitDistVec);

                            //calculate right hand side
                            rightEntry = (lambdaW * densityW * dV_w + lambdaNW * densityNW * dV_n) * (permeability * gravity)
                                    * faceArea ;
                        } //end 2p subdomain


                        // set diagonal entry and right hand side entry
                        A_[globalIdxI][globalIdxI] += entry;
                        f_[globalIdxI] += entry * pressBound;
                        f_[globalIdxI] -= rightEntry * (unitOuterNormal * unitDistVec);
                    }    //end of if(first) ... else{...
                }   // end dirichlet

                /**********************************
                 * set neumann boundary condition
                 **********************************/
                else
                {
                    Dune::FieldVector<Scalar,2> J = problem_.neumann(globalPosFace, *isIt);
                    if (first || problem_.variables().subdomain(globalIdxI)==1)
                    {
                        J[wPhaseIdx] /= densityWI;
                        J[nPhaseIdx] /= densityNWI;
                    }
                    else
                    {
                        J[wPhaseIdx] *= dv_dC1;
                        J[nPhaseIdx] *= dv_dC2;
                    }

                    f_[globalIdxI] -= (J[wPhaseIdx] + J[nPhaseIdx]) * faceArea;
                }
                /*************************************************/
            }
        } // end all intersections

        // compressibility term
        if (!first && timestep_ != 0.)
        {
            if (dv_dp[globalIdxI] == 0.)
            {
                // if incompressible fluids are used, the following has to be changed
                assert(problem_.variables().subdomain(globalIdxI)==1);

                // numerical derivative of fluid volume with respect to pressure
                Scalar p_ = problem_.variables().pressure()[globalIdxI] + 1e-2;
                Scalar Z1 = problem_.variables().totalConcentration(globalIdxI, wCompIdx)
                        / (problem_.variables().totalConcentration(globalIdxI, wCompIdx)
                                + problem_.variables().totalConcentration(globalIdxI, nCompIdx));
                PseudoOnePTwoCFluidState<TypeTag> pseudoFluidState;
                pseudoFluidState.update(Z1, p_, problem_.variables().saturation(globalIdxI));
                if(problem_.variables().saturation(globalIdxI)==1)  //only w-phase
                {
                    Scalar v_w_ = 1. / FluidSystem::phaseDensity(wPhaseIdx,
                                                        problem_.temperature(globalPos, *eIt),
                                                        p_, pseudoFluidState);
                    dv_dp[globalIdxI] = (problem_.variables().totalConcentration(globalIdxI, wCompIdx)
                            + problem_.variables().totalConcentration(globalIdxI, nCompIdx))
                            * ( v_w_  - 1./densityWI) / 1e-2;
                }
                if(problem_.variables().saturation(globalIdxI)==0.)  //only nw-phase
                {
                    Scalar v_n_ = 1. / FluidSystem::phaseDensity(nPhaseIdx,
                                                        problem_.temperature(globalPos, *eIt),
                                                        p_, pseudoFluidState);
                    dv_dp[globalIdxI] = (problem_.variables().totalConcentration(globalIdxI, wCompIdx)
                            + problem_.variables().totalConcentration(globalIdxI, nCompIdx))
                            * ( v_n_  - 1./densityNWI) / 1e-2;
                }
            }    //end calculation of 1p compress_term
            Scalar compress_term = dv_dp[globalIdxI] / timestep_;

            A_[globalIdxI][globalIdxI] -= compress_term*volume;
            f_[globalIdxI] -= problem_.variables().pressure()[globalIdxI] * compress_term * volume;

            if (isnan(compress_term) || isinf(compress_term))
                DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << globalIdxI);

            if(!GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)))
                DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
        }

        // error reduction routine: volumetric error is damped and inserted to right hand side
        // if damping is not done, the solution method gets unstable!
        problem_.variables().volErr()[globalIdxI] /= timestep_;
        Scalar erri = fabs(problem_.variables().volErr()[globalIdxI]);
        Scalar x_lo = 0.6;
        Scalar x_mi = 0.9;
        Scalar fac  = 0.05;
        Scalar lofac = 0.;
        Scalar hifac = 0.;
        hifac /= fac;

        if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxErr))
        {
            if (erri <= x_mi * maxErr)
                f_[globalIdxI] += errorCorrection[globalIdxI] = fac* (1-x_mi*(lofac/fac-1)/(x_lo-x_mi) + (lofac/fac-1)/(x_lo-x_mi)*erri/maxErr)
                                    * problem_.variables().volErr()[globalIdxI] * volume;
            else
                f_[globalIdxI] += errorCorrection[globalIdxI] = fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr)
                                    * problem_.variables().volErr()[globalIdxI] * volume;
        }
    } // end grid traversal
    return;
}

//! solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::solve()
{
    typedef typename GET_PROP(TypeTag, PTAG(SolverParameters)) SolverParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressurePreconditioner)) Preconditioner;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureSolver)) Solver;

//    printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//    printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);

    Dune::MatrixAdapter<Matrix, RHSVector, RHSVector> op(A_); // make linear operator from A_
    Dune::InverseOperatorResult result;
    Scalar reduction = SolverParameters::reductionSolver;
    int maxItSolver = SolverParameters::maxIterationNumberSolver;
    int iterPreconditioner = SolverParameters::iterationNumberPreconditioner;
    int verboseLevelSolver = SolverParameters::verboseLevelSolver;
    Scalar relaxation = SolverParameters::relaxationPreconditioner;

    if (verboseLevelSolver)
        std::cout << "FVPressure2P: solve for pressure" << std::endl;

    Preconditioner preconditioner(A_, iterPreconditioner, relaxation);
    Solver solver(op, preconditioner, reduction, maxItSolver, verboseLevelSolver);
    solver.apply(problem_.variables().pressure(), f_, result);

    return;
}

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
void FVPressure2P2CMultiPhysics<TypeTag>::initialMaterialLaws(bool compositional)
{
    // initialize the fluid system
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        // assign an Index for convenience
        int globalIdx = problem_.variables().index(*eIt);

        // get the temperature
        Scalar temperature_ = problem_.temperature(globalPos, *eIt);

        // initial conditions
        problem_.variables().capillaryPressure(globalIdx) = 0.;
        Scalar pressW = 0;
        Scalar pressNW = 0;
        Scalar sat_0=0.;

        BoundaryConditions2p2c::Flags ictype = problem_.initFormulation(globalPos, *eIt);            // get type of initial condition

        if(!compositional) //means that we do the first approximate guess without compositions
        {
            // phase pressures are unknown, so start with an exemplary
            Scalar exemplaryPressure = problem_.referencePressure(globalPos, *eIt);
            pressW = pressNW = problem_.variables().pressure()[globalIdx] = exemplaryPressure;

            if (ictype == BoundaryConditions2p2c::saturation)  // saturation initial condition
            {
                sat_0 = problem_.initSat(globalPos, *eIt);
                fluidState.satFlash(sat_0, pressW, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
            else if (ictype == BoundaryConditions2p2c::concentration) // concentration initial condition
            {
                Scalar Z1_0 = problem_.initConcentration(globalPos, *eIt);
                fluidState.update(Z1_0, pressW, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
        }
        else if(compositional)    //means we regard compositional effects since we know an estimate pressure field
        {
            //determine phase pressures from primary pressure variable
            switch (pressureType)
            {
                case pw:
                {
                    pressW = problem_.variables().pressure()[globalIdx];
                    pressNW = problem_.variables().pressure()[globalIdx] + problem_.variables().capillaryPressure(globalIdx);
                    break;
                }
                case pn:
                {
                    pressW = problem_.variables().pressure()[globalIdx] - problem_.variables().capillaryPressure(globalIdx);
                    pressNW = problem_.variables().pressure()[globalIdx];
                    break;
                }
            }

            if (ictype == Dumux::BoundaryConditions2p2c::saturation)  // saturation initial condition
            {
                sat_0 = problem_.initSat(globalPos, *eIt);
                fluidState.satFlash(sat_0, pressW, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
            else if (ictype == Dumux::BoundaryConditions2p2c::concentration) // concentration initial condition
            {
                Scalar Z1_0 = problem_.initConcentration(globalPos, *eIt);
                fluidState.update(Z1_0, pressW, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
            }
        }

        // initialize densities
        problem_.variables().densityWetting(globalIdx) = FluidSystem::phaseDensity(wPhaseIdx, temperature_, pressW, fluidState);
        problem_.variables().densityNonwetting(globalIdx) = FluidSystem::phaseDensity(nPhaseIdx, temperature_, pressNW, fluidState);

        // initialize mass fractions
        problem_.variables().wet_X1(globalIdx) = fluidState.massFrac(wPhaseIdx, wCompIdx);
        problem_.variables().nonwet_X1(globalIdx) = fluidState.massFrac(nPhaseIdx, wCompIdx);


        // initialize viscosities
        problem_.variables().viscosityWetting(globalIdx) = FluidSystem::phaseViscosity(wPhaseIdx, temperature_, pressW, fluidState);
        problem_.variables().viscosityNonwetting(globalIdx) = FluidSystem::phaseViscosity(nPhaseIdx, temperature_, pressNW, fluidState);

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx) = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem_.variables().viscosityWetting(globalIdx);
        problem_.variables().mobilityNonwetting(globalIdx) = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem_.variables().viscosityNonwetting(globalIdx);

        // initialize cell concentration
        problem_.variables().totalConcentration(globalIdx, wCompIdx) = fluidState.massConcentration(wCompIdx);
        problem_.variables().totalConcentration(globalIdx, nCompIdx) = fluidState.massConcentration(nCompIdx);
        problem_.variables().saturation()[globalIdx] = fluidState.saturation(wPhaseIdx);
    }
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
    // this method only completes the variables: = old postprocessupdate()

    // instantiate standard 2p2c and pseudo1p fluid state objects
    FluidState fluidState;
    PseudoOnePTwoCFluidState<TypeTag> pseudoFluidState;

    //get timestep for error term
    Scalar dt = problem_.timeManager().timeStepSize();

    // next subdomain map
    if (problem_.timeManager().time() == 0.)
        nextSubdomain = 2;  // start with complicated sub in initialization
    else
        nextSubdomain = 1;  // reduce complexity after first TS

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        int globalIdx = problem_.variables().index(*eIt);

        Scalar temperature_ = problem_.temperature(globalPos, *eIt);
        // reset volErr
        problem_.variables().volErr()[globalIdx] = 0;

        // get the overall mass of component 1:  Z1 = C^k / (C^1+C^2) [-]
        Scalar Z1 = problem_.variables().totalConcentration(globalIdx, wCompIdx)
                / (problem_.variables().totalConcentration(globalIdx, wCompIdx)
                        + problem_.variables().totalConcentration(globalIdx, nCompIdx));
        // make shure only physical quantities enter flash calculation
        #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
        if(Z1<0. || Z1 > 1.)
        {
            std::cout << "Feed mass fraction unphysical: Z1 = " << Z1
                   << " at global Idx " << globalIdx
                   << " , because totalConcentration(globalIdx, wCompIdx) = "
                   << problem_.variables().totalConcentration(globalIdx, wCompIdx)
                   << " and totalConcentration(globalIdx, nCompIdx) = "
                   << problem_.variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            if(Z1<0.)
            {
            Z1 = 0.;
            // add this error to volume error term for correction in next TS
            problem_.variables().volErr()[globalIdx] +=
                    problem_.variables().totalConcentration(globalIdx, wCompIdx)
                    / problem_.variables().densityWetting(globalIdx);
            //regul!
            problem_.variables().totalConcentration(globalIdx, wCompIdx) = 0.;
            Dune::dgrave << "Regularize totalConcentration(globalIdx, wCompIdx) = "
                << problem_.variables().totalConcentration(globalIdx, wCompIdx)<< std::endl;
            }
        else
            {
            Z1 = 1.;
            // add this error to volume error term for correction in next TS
            problem_.variables().volErr()[globalIdx] +=
                    problem_.variables().totalConcentration(globalIdx, nCompIdx)
                    / problem_.variables().densityNonwetting(globalIdx);
            //regul!
            problem_.variables().totalConcentration(globalIdx, nCompIdx) = 0.;
            Dune::dgrave << "Regularize totalConcentration(globalIdx, nCompIdx) = "
                << problem_.variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            }
        }
        #endif

        if (problem_.variables().subdomain(globalIdx)==2)   //=> 2p domain
        {
        //determine phase pressures from primary pressure variable
        Scalar pressW(0.), pressNW(0.);
        pressW =pressNW= problem_.variables().pressure()[globalIdx];

            //complete fluid state
            fluidState.update(Z1, pressW, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);

            /******** update variables in variableclass **********/
            // initialize saturation
            problem_.variables().saturation(globalIdx) = fluidState.saturation(wPhaseIdx);

            // the following makes shure pC is neglected
            problem_.variables().capillaryPressure(globalIdx) = 0.0;

            // initialize viscosities
            problem_.variables().viscosityWetting(globalIdx)
                    = FluidSystem::phaseViscosity(wPhaseIdx, temperature_, pressW, fluidState);
            problem_.variables().viscosityNonwetting(globalIdx)
                    = FluidSystem::phaseViscosity(nPhaseIdx, temperature_, pressNW, fluidState);

            // initialize mobilities
            problem_.variables().mobilityWetting(globalIdx) =
                    MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                        / problem_.variables().viscosityWetting(globalIdx);
            problem_.variables().mobilityNonwetting(globalIdx) =
                    MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                        / problem_.variables().viscosityNonwetting(globalIdx);

            // initialize mass fractions
            problem_.variables().wet_X1(globalIdx) = fluidState.massFrac(wPhaseIdx, wCompIdx);
            problem_.variables().nonwet_X1(globalIdx) = fluidState.massFrac(nPhaseIdx, wCompIdx);

            // initialize densities
            problem_.variables().densityWetting(globalIdx)
                    = FluidSystem::phaseDensity(wPhaseIdx, temperature_, pressW, fluidState);
            problem_.variables().densityNonwetting(globalIdx)
                    = FluidSystem::phaseDensity(nPhaseIdx, temperature_, pressNW, fluidState);

            problem_.spatialParameters().update(fluidState.saturation(wPhaseIdx), *eIt);

            // determine volume mismatch between actual fluid volume and pore volume
            Scalar sumConc = (problem_.variables().totalConcentration(globalIdx, wCompIdx)
                    + problem_.variables().totalConcentration(globalIdx, nCompIdx));
            Scalar massw = problem_.variables().numericalDensity(globalIdx, wPhaseIdx) = sumConc * fluidState.phaseMassFraction(wPhaseIdx);
            Scalar massn = problem_.variables().numericalDensity(globalIdx, nPhaseIdx) =sumConc * fluidState.phaseMassFraction(nPhaseIdx);
            Scalar vol = massw / problem_.variables().densityWetting(globalIdx)
                       + massn / problem_.variables().densityNonwetting(globalIdx);
            if (dt != 0)
            {
                problem_.variables().volErr()[globalIdx] += (vol - problem_.spatialParameters().porosity(globalPos, *eIt));

                Scalar volErrI = problem_.variables().volErr(globalIdx);
                if (std::isnan(volErrI))
                {
                    DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate:\n"
                            << "volErr[" << globalIdx << "] isnan: vol = " << vol
                            << ", massw = " << massw << ", rho_l = " << problem_.variables().densityWetting(globalIdx)
                            << ", massn = " << massn << ", rho_g = " << problem_.variables().densityNonwetting(globalIdx)
                            << ", poro = " << problem_.spatialParameters().porosity(globalPos, *eIt) << ", dt = " << dt);
                }
            }
            else
            {
                problem_.variables().volErr()[globalIdx] = 0;
            }

            // check subdomain consistency
            if (problem_.variables().saturation(globalIdx) != (1. || 0.)) // cell still 2p
            {
                // mark this element
                nextSubdomain[globalIdx] = 2;

                // mark neighbors
                IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt!=isItEnd; ++isIt)
                {
                    if (isIt->neighbor())
                    {
                        int globalIdxJ = problem_.variables().index(*(isIt->outside()));
                        // mark neighbor Element
                        nextSubdomain[globalIdxJ] = 2;
                    }
                }
                }
            // end subdomain check
        }
        else    // simple
        {
            // Todo: only variables of present phase need to be updated
            Scalar press1p = problem_.variables().pressure()[globalIdx];
            pseudoFluidState.update(Z1, press1p, problem_.variables().saturation(globalIdx));

            // initialize viscosities
            problem_.variables().viscosityWetting(globalIdx)
                    = FluidSystem::phaseViscosity(wPhaseIdx, temperature_, press1p, pseudoFluidState);
            problem_.variables().viscosityNonwetting(globalIdx)
                    = FluidSystem::phaseViscosity(nPhaseIdx, temperature_, press1p, pseudoFluidState);

            // initialize mobilities
            problem_.variables().mobilityWetting(globalIdx) =
                    MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), pseudoFluidState.saturation(wPhaseIdx))
                        / problem_.variables().viscosityWetting(globalIdx);
            problem_.variables().mobilityNonwetting(globalIdx) =
                    MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), pseudoFluidState.saturation(wPhaseIdx))
                        / problem_.variables().viscosityNonwetting(globalIdx);

            // initialize mass fractions
            problem_.variables().wet_X1(globalIdx) = pseudoFluidState.massFrac(wPhaseIdx, wCompIdx);
            problem_.variables().nonwet_X1(globalIdx) = pseudoFluidState.massFrac(nPhaseIdx, wCompIdx);

            // initialize densities
            problem_.variables().densityWetting(globalIdx)
                    = FluidSystem::phaseDensity(wPhaseIdx, temperature_, press1p, pseudoFluidState);
            problem_.variables().densityNonwetting(globalIdx)
                    = FluidSystem::phaseDensity(nPhaseIdx, temperature_, press1p, pseudoFluidState);

            // error term handling
            Scalar sumConc = (problem_.variables().totalConcentration(globalIdx, wCompIdx)
                    + problem_.variables().totalConcentration(globalIdx, nCompIdx));
            Scalar vol(0.);
            if(problem_.variables().saturation(globalIdx) == 1.)    //only w-phase
                vol = sumConc / problem_.variables().densityWetting(globalIdx);
            else    //only nw-phase
                vol = sumConc / problem_.variables().densityNonwetting(globalIdx);

            if (dt != 0)
                problem_.variables().volErr()[globalIdx] = (vol - problem_.spatialParameters().porosity(globalPos, *eIt));


        }
    }// end grid traversal

    problem_.variables().subdomain() = nextSubdomain;

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
 * \param[out] dv_dC1 partial derivative of fluid volume w.r.t. mass of component 1 \f$\mathrm{[m^3/kg]}\f$
 * \param[out] dv_dC2 partial derivative of fluid volume w.r.t. mass of component 2 \f$\mathrm{[m^3/kg]}\f$
 * \param[out] dv_dp partial derivative of fluid volume w.r.t. pressure \f$\mathrm{[1/Pa]}\f$
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::volumeDerivatives(GlobalPosition globalPos, ElementPointer ep, Scalar& dv_dC1, Scalar& dv_dC2, Scalar& dv_dp)
{
    // cell index
    int globalIdx = problem_.variables().index(*ep);

    // get cell temperature
    Scalar temperature_ = problem_.temperature(globalPos, *ep);

    // initialize an Fluid state for the update
    FluidState updFluidState;

    /**********************************
     * a) get necessary variables
     **********************************/
    //determine phase pressures from primary pressure variable
    Scalar pressW=0.;
    switch (pressureType)
    {
    case pw:
    {
        pressW = problem_.variables().pressure()[globalIdx];
        break;
    }
    case pn:
    {
        pressW = problem_.variables().pressure()[globalIdx] - problem_.variables().capillaryPressure(globalIdx);
        break;
    }
    }

    Scalar v_w = 1. / problem_.variables().densityWetting(globalIdx);
    Scalar v_g = 1. / problem_.variables().densityNonwetting(globalIdx);
    Scalar sati = problem_.variables().saturation()[globalIdx];
    // mass of components inside the cell
    Scalar m1 = problem_.variables().totalConcentration(globalIdx, wCompIdx);
    Scalar m2 = problem_.variables().totalConcentration(globalIdx, nCompIdx);
    // mass fraction of wetting phase
    Scalar nuw1 =  sati/ v_w / (sati/v_w + (1-sati)/v_g);
    // actual fluid volume
    Scalar volalt = (m1+m2) * (nuw1 * v_w + (1-nuw1) * v_g);

    /**********************************
     * b) define increments
     **********************************/
    // increments for numerical derivatives
    Scalar inc1 = (fabs(problem_.variables().updateEstimate(globalIdx, wCompIdx)) > 1e-8 / v_w) ?  problem_.variables().updateEstimate(globalIdx,wCompIdx) : 1e-8/v_w;
    Scalar inc2 =(fabs(problem_.variables().updateEstimate(globalIdx, nCompIdx)) > 1e-8 / v_g) ?  problem_.variables().updateEstimate(globalIdx,nCompIdx) : 1e-8 / v_g;
    Scalar incp = 1e-2;


    /**********************************
     * c) Secant method for derivatives
     **********************************/

    // numerical derivative of fluid volume with respect to pressure
    Scalar p_ = pressW + incp;
    Scalar Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1,
            p_, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
    Scalar v_w_ = 1. / FluidSystem::phaseDensity(wPhaseIdx, temperature_, p_, updFluidState);
    Scalar v_g_ = 1. / FluidSystem::phaseDensity(nPhaseIdx, temperature_, p_, updFluidState);
    dv_dp = ((m1+m2) * (nuw1 * v_w_ + (1-nuw1) * v_g_) - volalt) /incp;

    // numerical derivative of fluid volume with respect to mass of component 1
    m1 +=  inc1;
    Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1, pressW, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
    Scalar satt = updFluidState.saturation(wPhaseIdx);
    Scalar nuw = satt / v_w / (satt/v_w + (1-satt)/v_g);
    dv_dC1 = ((m1+m2) * (nuw * v_w + (1-nuw) * v_g) - volalt) /inc1;
    m1 -= inc1;

    // numerical derivative of fluid volume with respect to mass of component 2
    m2 += inc2;
    Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1, pressW, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
    satt = updFluidState.saturation(wPhaseIdx);
    nuw = satt / v_w / (satt/v_w + (1-satt)/v_g);
    dv_dC2 = ((m1+m2) * (nuw * v_w + (1-nuw) * v_g) - volalt)/ inc2;
    m2 -= inc2;
}


}//end namespace Dumux
#endif
