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
#include "dumux/common/pardiso.hh"
#include "dumux/common/math.hh"
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
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2P2C
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

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

protected:
    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

    Problem& problem()
    {
        return problem_;
    }
    const Problem& problem() const
    {
        return problem_;
    }

public:

    //function which assembles the system of equations to be solved
    void assemble(bool first);

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

        problem_.variables().communicateUpdateEstimate();

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

#if DUNE_MINIMAL_DEBUG_LEVEL <= 2
        int size_ = problem_.variables().gridSize();
        // add debug stuff
        typename SolutionTypes::ScalarSolution *numdensityW = writer.allocateManagedBuffer (size_);
        typename SolutionTypes::ScalarSolution *numdensityNW = writer.allocateManagedBuffer (size_);
        *numdensityW = problem_.variables().numericalDensity(wPhaseIdx);
        *numdensityNW = problem_.variables().numericalDensity(nPhaseIdx);
        writer.attachCellData(*numdensityW, "numerical density (mass/volume) w_phase");
        writer.attachCellData(*numdensityNW, "numerical density (mass/volume) nw_phase");

        typename SolutionTypes::ScalarSolution *errorCorrPtr = writer.allocateManagedBuffer (size_);
        *errorCorrPtr = problem_.variables().errorCorrection();
//        int size = subdomainPtr.size();
        writer.attachCellData(*errorCorrPtr, "Error Correction");

        typename SolutionTypes::ScalarSolution *dv_dpPtr = writer.allocateManagedBuffer (size_);
        *dv_dpPtr = problem_.variables().dv_dp();
        writer.attachCellData(*dv_dpPtr, "dv_dp");

        typename SolutionTypes::ScalarSolution *dV_dC1Ptr = writer.allocateManagedBuffer (size_);
        typename SolutionTypes::ScalarSolution *dV_dC2Ptr = writer.allocateManagedBuffer (size_);
        *dV_dC1Ptr = problem_.variables().dv()[0];
        *dV_dC2Ptr = problem_.variables().dv()[1];
        writer.attachCellData(*dV_dC1Ptr, "dV_dC1");
        writer.attachCellData(*dV_dC2Ptr, "dV_dC2");

        Dune::BlockVector<Dune::FieldVector<double,1> > *updEstimate1 = writer.allocateManagedBuffer (size_);
        Dune::BlockVector<Dune::FieldVector<double,1> > *updEstimate2 = writer.allocateManagedBuffer (size_);
        *updEstimate1 = problem_.variables().updateEstimate()[0];
        *updEstimate2 = problem_.variables().updateEstimate()[1];
        writer.attachCellData(*updEstimate1, "updEstimate comp 1");
        writer.attachCellData(*updEstimate2, "updEstimate comp 2");
#endif
        return;
    }

    //! \brief Write additional debug info in a special writer
    // used via pseudoTS thorugh the initialization procedure.
    void debugOutput(double pseudoTS = 0.)
    {
        std::cout << "Writing debug for current time step\n";
        debugWriter_.beginWrite(problem_.timeManager().time() + pseudoTS);
        addOutputVtkFields(debugWriter_);

#if DUNE_MINIMAL_DEBUG_LEVEL <= 2
        int size_ = problem_.variables().gridSize();
        // output porosity, permeability
        Dune::BlockVector<Dune::FieldVector<double,1> > *poroPtr = debugWriter_.allocateManagedBuffer (size_);
        Dune::BlockVector<Dune::FieldVector<double,1> > *permPtr = debugWriter_.allocateManagedBuffer (size_);

        Dune::BlockVector<Dune::FieldVector<double,1> > poro_(0.), perm_(0.);
        poro_.resize(size_); perm_.resize(size_);
        // iterate over all elements of domain
        for (ElementIterator eIt = problem_.gridView().template begin<0> ();
                eIt != problem_.gridView().template end<0>(); ++eIt)
        {
            // get position, index
            GlobalPosition globalPos = eIt->geometry().center();
            int globalIdx = problem_.variables().index(*eIt);
            poro_[globalIdx] = problem().spatialParameters().porosity(globalPos, *eIt);
            perm_[globalIdx] = problem().spatialParameters().intrinsicPermeability(globalPos, *eIt)[0][0];
        }
        *poroPtr = poro_;
        *permPtr = perm_;
        debugWriter_.attachCellData(*poroPtr, "porosity");
        debugWriter_.attachCellData(*permPtr, "permeability");
#endif

        debugWriter_.endWrite();
        return;
    }
    //@}

    //! Constructs a FVPressure2P2C object
    /**
     * \param problem a problem class object
     */
    FVPressure2P2C(Problem& problem) :
        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (2 * dim + 1)
                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()),
                debugWriter_(problem.gridView(),"debugOutput2p2c")
    {
        cFLFactor_ = GET_PARAM(TypeTag, Scalar, CFLFactor);

        if (pressureType != pw && pressureType != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
    }

protected:
    Problem& problem_;
    Matrix A_;
    RHSVector f_;
    std::string solverName_;
    std::string preconditionerName_;

    // debug
    Dumux::VtkMultiWriter<GridView> debugWriter_;

    Scalar cFLFactor_; //!< determines the CFLfactor
    static constexpr int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
};

//! initializes the matrix to store the system of equations
template<class TypeTag>
void FVPressure2P2C<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // set perimeter to zero
        problem().variables().perimeter(globalIdxI) = 0;

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            problem().variables().perimeter(globalIdxI)
                    += isIt->geometry().volume();
            if (isIt->neighbor())
                rowSize++;
        }
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
void FVPressure2P2C<TypeTag>::initialize(bool solveTwice)
{
    //prepare stiffness Matrix and Vectors
    problem_.pressureModel().initializeMatrix();

    // initialguess: set saturations, determine visco and mobility for initial pressure equation
    // at this moment, the pressure is unknown. Hence, dont regard compositional effects.
    Dune::dinfo << "first saturation guess"<<std::endl; //=J: initialGuess()
        problem_.pressureModel().initialMaterialLaws(false);
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                debugOutput();
            #endif
    Dune::dinfo << "first pressure guess"<<std::endl;
        problem_.pressureModel().assemble(true);
        problem_.pressureModel().solve();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                debugOutput(1e-6);
            #endif
    // update the compositional variables (hence true)
    Dune::dinfo << "first guess for mass fractions"<<std::endl;//= J: transportInitial()
        problem_.pressureModel().initialMaterialLaws(true);

    // perform concentration update to determine secants
    Dune::dinfo << "secant guess"<< std::endl;
        Scalar dt_estimate = 0.;
        problem_.transportModel().update(0., dt_estimate, problem_.variables().updateEstimate(), false);
        dt_estimate = std::min ( problem_.timeManager().timeStepSize(), dt_estimate);
        //make sure the right time-step is used by all processes in the parallel case
        #if HAVE_MPI
        dt_estimate = problem_.gridView().comm().min(dt_estimate);
        #endif

        problem_.variables().updateEstimate() *= dt_estimate;
        //communicate in parallel case
        problem_.variables().communicateUpdateEstimate();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                debugOutput(2e-6);
            #endif
    // pressure calculation
    Dune::dinfo << "second pressure guess"<<std::endl;
        problem_.pressureModel().assemble(false);
        problem_.pressureModel().solve();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                debugOutput(3e-6);
            #endif

        // update the compositional variables
        problem_.pressureModel().initialMaterialLaws(true);


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
            problem_.transportModel().update(0., dt_dummy, problem_.variables().updateEstimate(), false);   Dune::dinfo << "secant guess"<< std::endl;
            problem_.variables().updateEstimate() *= dt_estimate;

            //communicate in parallel case
            problem_.variables().communicateUpdateEstimate();

            // pressure calculation
            problem_.pressureModel().assemble(false);                 Dune::dinfo << "pressure guess number "<< numIter <<std::endl;
            problem_.pressureModel().solve();
            // update the compositional variables
            problem_.pressureModel().initialMaterialLaws(true);

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
void FVPressure2P2C<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    // initialization: set the fluid volume derivatives to zero
    Scalar dv_dC1(0.), dv_dC2(0.);
    for (int i = 0; i < problem_.variables().gridSize(); i++)
    {
        problem_.variables().dv(i, wPhaseIdx) = 0;     // dv_dC1 = dV/dm1
        problem_.variables().dv(i, nPhaseIdx) = 0;     // dv / dC2
        problem_.variables().dv_dp(i) = 0;      // dv / dp
    }

    // determine maximum error to scale error-term
    Scalar timestep_ = problem_.timeManager().timeStepSize();
    Scalar maxErr = fabs(problem_.variables().volErr().infinity_norm()) / timestep_;

    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().center();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // cell volume & perimeter, assume linear map here
        Scalar volume = eIt->geometry().volume();
        Scalar perimeter = problem().variables().perimeter(globalIdxI);

        const GlobalPosition& gravity_ = problem().gravity();

        Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

        /****************implement sources and sinks************************/
        Dune::FieldVector<Scalar,2> source(problem_.source(globalPos, *eIt));
        if(first)
        {
                source[wPhaseIdx] /= densityWI;
                source[nPhaseIdx] /= densityNWI;
        }
        else
        {
		        // derivatives of the fluid volume with respect to concentration of components, or pressure
				if (problem_.variables().dv(globalIdxI, wCompIdx) == 0)
					volumeDerivatives(globalPos, *eIt,
							problem_.variables().dv(globalIdxI, wPhaseIdx),
							problem_.variables().dv(globalIdxI, nPhaseIdx),
							problem_.variables().dv_dp(globalIdxI));

				source[wPhaseIdx] *= problem_.variables().dv(globalIdxI, wPhaseIdx);		// note: dV_[i][1] = dv_dC1 = dV/dm1
				source[nPhaseIdx] *= problem_.variables().dv(globalIdxI, nPhaseIdx);
        }
		f_[globalIdxI] = volume * (source[wPhaseIdx] + source[nPhaseIdx]);
        /***********************************/

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
            // get normal vector
            const GlobalPosition& unitOuterNormal = isIt->centerUnitOuterNormal();

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                // gemotry info of neighbor
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                GlobalPosition unitDistVec(distVec);
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
                    rightEntry = factor * lambda * faceArea * (permeability * gravity_) * (unitOuterNormal * unitDistVec);
                }
                else
                {
                	// determine volume derivatives
                    if (problem_.variables().dv(globalIdxJ, wCompIdx) == 0)
                    	volumeDerivatives(globalPosNeighbor, *neighborPointer,
                    			problem_.variables().dv(globalIdxJ, wPhaseIdx),
                    			problem_.variables().dv(globalIdxJ, nPhaseIdx),
                    			problem_.variables().dv_dp(globalIdxJ));
                    dv_dC1 = (problem_.variables().dv(globalIdxJ, wPhaseIdx)
                    			+ problem_.variables().dv(globalIdxI, wPhaseIdx)) / 2; // dV/dm1= dv/dC^1
                    dv_dC2 = (problem_.variables().dv(globalIdxJ, nPhaseIdx)
                    			+ problem_.variables().dv(globalIdxI, nPhaseIdx)) / 2;

                    Scalar graddv_dC1 = (problem_.variables().dv(globalIdxJ, wPhaseIdx)
                    						- problem_.variables().dv(globalIdxI, wPhaseIdx)) / dist;
                    Scalar graddv_dC2 = (problem_.variables().dv(globalIdxJ, nPhaseIdx)
                    						- problem_.variables().dv(globalIdxI, nPhaseIdx)) / dist;


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

                    switch (pressureType)
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
                    }

                    potentialW += densityW * (unitDistVec * gravity_);
                    potentialNW += densityNW * (unitDistVec * gravity_);

					// initialize convenience shortcuts
					Scalar lambdaW, lambdaN;
					Scalar dV_w(0.), dV_n(0.);		// dV_a = \sum_k \rho_a * dv/dC^k * X^k_a
					Scalar gV_w(0.), gV_n(0.);		// multipaper eq(3.3) line 3 analogon dV_w


					//do the upwinding of the mobility depending on the phase potentials
					if (potentialW >= 0.)
					{
						dV_w = (dv_dC1 * problem_.variables().wet_X1(globalIdxI) + dv_dC2 * (1. - problem_.variables().wet_X1(globalIdxI)));
						lambdaW = problem_.variables().mobilityWetting(globalIdxI);
						gV_w = (graddv_dC1 * problem_.variables().wet_X1(globalIdxI) + graddv_dC2 * (1. - problem_.variables().wet_X1(globalIdxI)));
						dV_w *= densityWI; gV_w *= densityWI;
					}
					else
					{
						dV_w = (dv_dC1 * problem_.variables().wet_X1(globalIdxJ) + dv_dC2 * (1. - problem_.variables().wet_X1(globalIdxJ)));
						lambdaW = problem_.variables().mobilityWetting(globalIdxJ);
						gV_w = (graddv_dC1 * problem_.variables().wet_X1(globalIdxJ) + graddv_dC2 * (1. - problem_.variables().wet_X1(globalIdxJ)));
						dV_w *= densityWJ; gV_w *= densityWJ;
					}
					if (potentialNW >= 0.)
					{
						dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxI) + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
						lambdaN = problem_.variables().mobilityNonwetting(globalIdxI);
						gV_n = (graddv_dC1 * problem_.variables().nonwet_X1(globalIdxI) + graddv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
						dV_n *= densityNWI; gV_n *= densityNWI;
					}
					else
					{
						dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxJ) + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxJ)));
						lambdaN = problem_.variables().mobilityNonwetting(globalIdxJ);
						gV_n = (graddv_dC1 * problem_.variables().nonwet_X1(globalIdxJ) + graddv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxJ)));
						dV_n *= densityNWJ; gV_n *= densityNWJ;
					}

	                //calculate current matrix entry
                    entry = faceArea * (lambdaW * dV_w + lambdaN * dV_n)
                            - volume * faceArea / perimeter * (lambdaW * gV_w + lambdaN * gV_n); 	// = boundary integral - area integral
                    entry *= fabs((permeability*unitOuterNormal)/(dist));

                    //calculate right hand side
                    rightEntry = faceArea  * (unitOuterNormal * unitDistVec) * (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n);
                	rightEntry -= volume * faceArea / perimeter * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n);
                    rightEntry *= (permeability * gravity_); 		// = multipaper eq(3.3) line 2+3

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
                dv_dC1 = problem_.variables().dv(globalIdxI, wCompIdx);
                dv_dC2 = problem_.variables().dv(globalIdxI, nCompIdx);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                // geometrical information
                GlobalPosition distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = problem_.bcTypePress(globalPosFace, *isIt);

                // prepare pressure boundary condition
                PhaseVector pressBC(0.);
                Scalar pcBound (0.);

                /**********         Dirichlet Boundary        *************/
                if (bctype == BoundaryConditions::dirichlet)
                {
                    //permeability vector at boundary
                    Dune::FieldVector<Scalar, dim> permeability(0);
                    permeabilityI.mv(unitDistVec, permeability);

                	// create a fluid state for the boundary
                	FluidState BCfluidState;

                    Scalar temperatureBC = problem_.temperature(globalPosFace, *eIt);

                    if (first)
                    {
                    	Scalar lambda = lambdaWI+lambdaNWI;
                        A_[globalIdxI][globalIdxI] += lambda * faceArea * (permeability * unitOuterNormal) / (dist);
                        Scalar pressBC = problem_.dirichletPress(globalPosFace, *isIt);
                        f_[globalIdxI] += lambda * faceArea * pressBC * (permeability * unitOuterNormal) / (dist);
                        Scalar rightentry = (fractionalWI * densityWI
                                             + fractionalNWI * densityNWI)
                                             * lambda * faceArea * (unitOuterNormal * unitDistVec) *(permeability*gravity_);
                        f_[globalIdxI] -= rightentry;
                    }
                    else    //not first
                    {
                        // read boundary values
                        problem_.transportModel().evalBoundary(globalPosFace,
                                                                    isIt,
                                                                    eIt,
                                                                    BCfluidState,
                                                                    pressBC);
                        pcBound = pressBC[nPhaseIdx] - pressBC[wPhaseIdx];

                        // determine fluid properties at the boundary
                        Scalar lambdaWBound = 0.;
                        Scalar lambdaNWBound = 0.;

                        Scalar densityWBound =
                            FluidSystem::phaseDensity(wPhaseIdx, temperatureBC, pressBC[wPhaseIdx], BCfluidState);
                        Scalar densityNWBound =
                            FluidSystem::phaseDensity(nPhaseIdx, temperatureBC, pressBC[nPhaseIdx], BCfluidState);
                        Scalar viscosityWBound =
                            FluidSystem::phaseViscosity(wPhaseIdx, temperatureBC, pressBC[wPhaseIdx], BCfluidState);
                        Scalar viscosityNWBound =
                            FluidSystem::phaseViscosity(nPhaseIdx, temperatureBC, pressBC[nPhaseIdx], BCfluidState);

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
//                            potentialW = problem_.variables().potentialWetting(globalIdxI, isIndex);
//                            potentialNW = problem_.variables().potentialNonwetting(globalIdxI, isIndex);
//
//                            // do potential upwinding according to last potGradient vs Jochen: central weighting
//                            densityW = (potentialW > 0.) ? densityWI : densityWBound;
//                            densityNW = (potentialNW > 0.) ? densityNWI : densityNWBound;
//
//                            densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                            densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                            densityW=rhoMeanW; densityNW=rhoMeanNW;

                            //calculate potential gradient
                            switch (pressureType)
                            {
                                case pw:
                                {
                                    potentialW = (unitOuterNormal * distVec) *
                                            (problem_.variables().pressure()[globalIdxI] - problem_.dirichletPress(globalPosFace, *isIt)) / (dist * dist);
                                    potentialNW = (unitOuterNormal * distVec) *
                                            (problem_.variables().pressure()[globalIdxI] + pcI
                                                    - problem_.dirichletPress(globalPosFace, *isIt) - pcBound) / (dist * dist);
                                    break;
                                }
                                case pn:
                                {
                                    potentialW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] - pcI -
                                            problem_.dirichletPress(globalPosFace, *isIt) + pcBound)
                                            / (dist * dist);
                                    potentialNW = (unitOuterNormal * distVec) * (problem_.variables().pressure()[globalIdxI] -
                                            problem_.dirichletPress(globalPosFace, *isIt)) / (dist * dist);
                                    break;
                                }
                            }
                            potentialW += densityW * (unitDistVec * gravity_);
                            potentialNW += densityNW * (unitDistVec * gravity_);
                       }   //end !first

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW, lambdaNW;
                        Scalar dV_w, dV_n; 	// gV_a weglassen, da dV/dc am Rand ortsunabhängig angenommen -> am rand nicht bestimmbar -> nur Randintegral ohne Gebietsintegral

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
                            dV_w = (dv_dC1 * BCfluidState.massFrac(wPhaseIdx, wCompIdx)
                            		 + dv_dC2 * BCfluidState.massFrac(wPhaseIdx, nCompIdx));
                            dV_w *= densityW;
                            lambdaW = lambdaWBound;
                        }
                        if (potentialNW >= 0.)
                        {
                            densityNW = (potentialNW == 0) ? rhoMeanNW : densityNWI;
                            dV_n = (dv_dC1 * problem_.variables().nonwet_X1(globalIdxI) + dv_dC2 * (1. - problem_.variables().nonwet_X1(globalIdxI)));
                            dV_n *= densityNW;
                            lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNWI;
                        }
                        else
                        {
                            densityNW = densityNWBound;
                            dV_n = (dv_dC1 * BCfluidState.massFrac(nPhaseIdx, wCompIdx)
                            		+ dv_dC2 * BCfluidState.massFrac(nPhaseIdx, nCompIdx));
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


                        // set diagonal entry and right hand side entry
                        A_[globalIdxI][globalIdxI] += entry;
                        f_[globalIdxI] += entry * problem_.dirichletPress(globalPosFace, *isIt);
                        f_[globalIdxI] -= rightEntry * (unitOuterNormal * unitDistVec);
                    }	//end of if(first) ... else{...
                }   // end dirichlet

                /**********************************
                 * set neumann boundary condition
                 **********************************/
                else
                {
                	Dune::FieldVector<Scalar,2> J = problem_.neumann(globalPosFace, *isIt);
					if (first)
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
            }
        } // end all intersections
//        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//        printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
        // compressibility term
        if (!first && timestep_ != 0)
        {
            Scalar compress_term = problem_.variables().dv_dp(globalIdxI) / timestep_;

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
        Scalar x_lo = 0.2;
        Scalar x_mi = 0.9;
        Scalar fac  = 0.5;
        if (pressureType == pw)
            fac = 0.05;
        Scalar lofac = 0.;
        Scalar hifac = 0.;

        if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxErr) && (!problem_.timeManager().willBeFinished()))
        {
            if (erri <= x_mi * maxErr)
                f_[globalIdxI] +=
                		problem_.variables().errorCorrection(globalIdxI) =
                				fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/maxErr)
                                    * problem_.variables().volErr()[globalIdxI] * volume;
            else
                f_[globalIdxI] +=
                		problem_.variables().errorCorrection(globalIdxI) =
                				fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr)
                                    * problem_.variables().volErr()[globalIdxI] * volume;
        }
//        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
    } // end grid traversal
//    printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
    return;
}

//! solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure2P2C<TypeTag>::solve()
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LinearSolver)) Solver;

    int verboseLevelSolver = GET_PARAM(TypeTag, int, LSVerbosity);

    if (verboseLevelSolver)
        std::cout << "compositional 2p2c: solve for pressure" << std::endl;

    Solver solver(problem_);
    solver.solve(A_, problem_.variables().pressure(), f_);
    //                printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
    //                printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
    //                printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);

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
void FVPressure2P2C<TypeTag>::initialMaterialLaws(bool compositional)
{
    problem().variables().communicateTransportedQuantity();
    problem().variables().communicatePressure();

	// initialize the fluid system
    FluidState fluidState;

	// iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    ElementIterator eIt = problem_.gridView().template begin<0>();
    for (; eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        // assign an Index for convenience
        int globalIdx = problem_.variables().index(*eIt);

        // get the temperature
        Scalar temperature_ = problem_.temperature(globalPos, *eIt);

        // initial conditions
		PhaseVector pressure(0.);
		Scalar sat_0=0.;

        BoundaryConditions2p2c::Flags ictype = problem_.initFormulation(globalPos, *eIt);            // get type of initial condition

        if(!compositional) //means that we do the first approximate guess without compositions
        {
			// phase pressures are unknown, so start with an exemplary
        	Scalar exemplaryPressure = problem_.referencePressure(globalPos, *eIt);
			pressure[wPhaseIdx] = pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx] = exemplaryPressure;
			problem_.variables().capillaryPressure(globalIdx) = 0.;
			if (ictype == BoundaryConditions2p2c::saturation)  // saturation initial condition
			{
				sat_0 = problem_.initSat(globalPos, *eIt);
				fluidState.satFlash(sat_0, pressure, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
			}
			else if (ictype == BoundaryConditions2p2c::concentration) // concentration initial condition
			{
				Scalar Z1_0 = problem_.initConcentration(globalPos, *eIt);
				fluidState.update(Z1_0, pressure, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
			}
        }
        else if(compositional)	//means we regard compositional effects since we know an estimate pressure field
        {
			if (ictype == Dumux::BoundaryConditions2p2c::saturation)  // saturation initial condition
			{
			    //get saturation, determine pc
				sat_0 = problem_.initSat(globalPos, *eIt);
		        if(GET_PROP_VALUE(TypeTag, PTAG(EnableCapillarity)))
		        {
                    problem_.variables().capillaryPressure(globalIdx)
                            = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
                                    sat_0);
		        }
		        else
		            problem_.variables().capillaryPressure(globalIdx) = 0.;

		        //determine phase pressures from primary pressure variable
                switch (pressureType)
                {
                    case pw:
                    {
                        pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx];
                        pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx] + problem_.variables().capillaryPressure(globalIdx);
                        break;
                    }
                    case pn:
                    {
                        pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx] - problem_.variables().capillaryPressure(globalIdx);
                        pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx];
                        break;
                    }
                }

				fluidState.satFlash(sat_0, pressure, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);
			}
			else if (ictype == Dumux::BoundaryConditions2p2c::concentration) // concentration initial condition
			{
				Scalar Z1_0 = problem_.initConcentration(globalPos, *eIt);
				// If total concentrations are given at the boundary, saturation is unknown.
				// This may affect pc and hence p_alpha and hence again saturation -> iteration.

		        // iterations in case of enabled capillary pressure
		        if(GET_PROP_VALUE(TypeTag, PTAG(EnableCapillarity)))
		        {
		            //start with pc from last TS
		            Scalar pc(problem_.variables().capillaryPressure(globalIdx));

		            int maxiter = 3;
		            //start iteration loop
		            for(int iter=0; iter < maxiter; iter++)
		            {
		                //determine phase pressures from primary pressure variable
		                switch (pressureType)
		                {
		                    case pw:
		                    {
		                        pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx];
		                        pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx] + pc;
		                        break;
		                    }
		                    case pn:
		                    {
		                        pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx] - pc;
		                        pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx];
		                        break;
		                    }
		                }

		                //store old pc
		                Scalar oldPc = pc;
		                //update with better pressures
		                fluidState.update(Z1_0, pressure, problem_.spatialParameters().porosity(globalPos, *eIt),
		                                    problem_.temperature(globalPos, *eIt));
		                pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
		                                    fluidState.saturation(wPhaseIdx));
		                // TODO: get right criterion, do output for evaluation
		                //converge criterion
		                if (abs(oldPc-pc)<10)
		                    iter = maxiter;

                        pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
                                fluidState.saturation(wPhaseIdx));
		            }
		            problem_.variables().capillaryPressure(globalIdx) = pc; //complete iteration procedure
		        }
		        else  // capillary pressure neglected
		        {
		            problem_.variables().capillaryPressure(globalIdx) = 0.;
		            pressure[wPhaseIdx] = pressure[nPhaseIdx]
		                = problem_.variables().pressure()[globalIdx];
		            fluidState.update(Z1_0, pressure, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);

		        }
			} //end conc initial condition
        } //end compositional

        // initialize densities
        problem_.variables().densityWetting(globalIdx) = FluidSystem::phaseDensity(wPhaseIdx, temperature_, pressure[wPhaseIdx], fluidState);
        problem_.variables().densityNonwetting(globalIdx) = FluidSystem::phaseDensity(nPhaseIdx, temperature_, pressure[nPhaseIdx], fluidState);

        // initialize mass fractions
        problem_.variables().wet_X1(globalIdx) = fluidState.massFrac(wPhaseIdx, wCompIdx);
        problem_.variables().nonwet_X1(globalIdx) = fluidState.massFrac(nPhaseIdx, wCompIdx);


        // initialize viscosities
        problem_.variables().viscosityWetting(globalIdx) = FluidSystem::phaseViscosity(wPhaseIdx, temperature_, pressure[wPhaseIdx], fluidState);
        problem_.variables().viscosityNonwetting(globalIdx) = FluidSystem::phaseViscosity(nPhaseIdx, temperature_, pressure[nPhaseIdx], fluidState);

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx) = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem_.variables().viscosityWetting(globalIdx);
        problem_.variables().mobilityNonwetting(globalIdx) = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                / problem_.variables().viscosityNonwetting(globalIdx);

        // initialize cell concentration
        problem_.variables().totalConcentration(globalIdx, wCompIdx) = fluidState.massConcentration(wCompIdx);
        problem_.variables().totalConcentration(globalIdx, nCompIdx) = fluidState.massConcentration(nCompIdx);
        problem_.variables().saturation()[globalIdx] = fluidState.saturation(wPhaseIdx);

        // to prevent output errors, set derivatives 0
        problem_.variables().dv(globalIdx, wCompIdx) = 0.;     // dv_dC1 = dV/dm1
        problem_.variables().dv(globalIdx, nCompIdx) = 0.;     // dv / dC2
        problem_.variables().dv_dp(globalIdx) = 0.;      // dv / dp
    }
    return;
}

//! constitutive functions are updated once if new concentrations are calculated and stored in the variables container
template<class TypeTag>
void FVPressure2P2C<TypeTag>::updateMaterialLaws()
{
    problem_.variables().communicateTransportedQuantity();
    problem_.variables().communicatePressure();

	// instantiate a brandnew fluid state object
    FluidState fluidState;

    //get timestep for error term
    Scalar dt = problem_.timeManager().timeStepSize();

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        int globalIdx = problem_.variables().index(*eIt);

        Scalar temperature_ = problem_.temperature(globalPos, *eIt);
        // reset volume error
        problem_.variables().volErr()[globalIdx] = 0;


        // get the overall mass of component 1 Z1 = C^k / (C^1+C^2) [-]
        Scalar Z1 = problem_.variables().totalConcentration(globalIdx, wCompIdx)
                / (problem_.variables().totalConcentration(globalIdx, wCompIdx)
                        + problem_.variables().totalConcentration(globalIdx, nCompIdx));

        // make shure only physical quantities enter flash calculation
        #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
        if(Z1<0. || Z1 > 1.)
        {
            Dune::dgrave << "Feed mass fraction unphysical: Z1 = " << Z1
                   << " at global Idx " << globalIdx
                   << " , because totalConcentration(globalIdx, wCompIdx) = "
                   << problem_.variables().totalConcentration(globalIdx, wCompIdx)
                   << " and totalConcentration(globalIdx, nCompIdx) = "
                   << problem_.variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
            if(Z1<0.)
                {
                Z1 = 0.;
                problem_.variables().totalConcentration(globalIdx, wCompIdx) = 0.;
                Dune::dgrave << "Regularize totalConcentration(globalIdx, wCompIdx) = "
                    << problem_.variables().totalConcentration(globalIdx, wCompIdx)<< std::endl;
                }
            else
                {
                Z1 = 1.;
                problem_.variables().totalConcentration(globalIdx, nCompIdx) = 0.;
                Dune::dgrave << "Regularize totalConcentration(globalIdx, nCompIdx) = "
                    << problem_.variables().totalConcentration(globalIdx, nCompIdx)<< std::endl;
                }
        }
        #endif

        //determine phase pressures from primary pressure variable and pc of last TS
        PhaseVector pressure(0.);
        switch (pressureType)
        {
        case pw:
        {
            pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx];
            pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx]
                      + problem_.variables().capillaryPressure(globalIdx);
            break;
        }
        case pn:
        {
            pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx]
                     - problem_.variables().capillaryPressure(globalIdx);
            pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx];
            break;
        }
        }

        //complete fluid state
        fluidState.update(Z1, pressure, problem_.spatialParameters().porosity(globalPos, *eIt), temperature_);

        // iterations part in case of enabled capillary pressure
        Scalar pc(0.), oldPc(0.);
        if(GET_PROP_VALUE(TypeTag, PTAG(EnableCapillarity)))
        {
            pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
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
                fluidState.update(Z1, pressure, problem_.spatialParameters().porosity(globalPos, *eIt),
                                    problem_.temperature(globalPos, *eIt));
                pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
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
		problem_.variables().saturation(globalIdx) = fluidState.saturation(wPhaseIdx);
		problem_.variables().capillaryPressure(globalIdx) = pc; // pc=0 if EnableCapillarity=false

        // initialize viscosities
		problem_.variables().viscosityWetting(globalIdx)
        		= FluidSystem::phaseViscosity(wPhaseIdx, temperature_, pressure[wPhaseIdx], fluidState);
		problem_.variables().viscosityNonwetting(globalIdx)
        		= FluidSystem::phaseViscosity(nPhaseIdx, temperature_, pressure[nPhaseIdx], fluidState);

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
        		= FluidSystem::phaseDensity(wPhaseIdx, temperature_, pressure[wPhaseIdx], fluidState);
        problem_.variables().densityNonwetting(globalIdx)
        		= FluidSystem::phaseDensity(nPhaseIdx, temperature_, pressure[nPhaseIdx], fluidState);

        // determine volume mismatch between actual fluid volume and pore volume
        Scalar sumConc = (problem_.variables().totalConcentration(globalIdx, wCompIdx)
                + problem_.variables().totalConcentration(globalIdx, nCompIdx));
        Scalar massw = problem_.variables().numericalDensity(globalIdx, wPhaseIdx) = sumConc * fluidState.phaseMassFraction(wPhaseIdx);
        Scalar massn = problem_.variables().numericalDensity(globalIdx, nPhaseIdx) = sumConc * fluidState.phaseMassFraction(nPhaseIdx);

        if ((problem_.variables().densityWetting(globalIdx)*problem_.variables().densityNonwetting(globalIdx)) == 0)
        	DUNE_THROW(Dune::MathError, "Decoupled2p2c::postProcessUpdate: try to divide by 0 density");
        Scalar vol = massw / problem_.variables().densityWetting(globalIdx) + massn / problem_.variables().densityNonwetting(globalIdx);
        if (dt != 0)
        {
            problem_.variables().volErr()[globalIdx] = (vol - problem_.spatialParameters().porosity(globalPos, *eIt));

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
            problem_.variables().volErr()[globalIdx] = 0;
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
    int globalIdx = problem_.variables().index(*ep);

    // get cell temperature
    Scalar temperature_ = problem_.temperature(globalPos, *ep);

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
            pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx];
            pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx] + problem_.variables().capillaryPressure(globalIdx);
            break;
        }
        case pn:
        {
            pressure[wPhaseIdx] = problem_.variables().pressure()[globalIdx] - problem_.variables().capillaryPressure(globalIdx);
            pressure[nPhaseIdx] = problem_.variables().pressure()[globalIdx];
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
    PhaseVector p_(incp);
    p_ += pressure;
    Scalar Z1 = m1 / (m1 + m2);
    updFluidState.update(Z1,
            p_, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
    Scalar v_w_ = 1. / FluidSystem::phaseDensity(wPhaseIdx, temperature_, p_[wPhaseIdx], updFluidState);
    Scalar v_g_ = 1. / FluidSystem::phaseDensity(nPhaseIdx, temperature_, p_[nPhaseIdx], updFluidState);
    dv_dp = ((m1+m2) * (nuw1 * v_w_ + (1-nuw1) * v_g_) - volalt) /incp;
    if (dv_dp>0)
    {
        // dV_dp > 0 is unphysical: Try inverse increment for secant
        Dune::dinfo << "dv_dp larger 0 at Idx " << globalIdx << " , try and invert secant"<< std::endl;

        p_ -= 2*incp;
        updFluidState.update(Z1,
                    p_, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
        Scalar v_w_ = 1. / FluidSystem::phaseDensity(wPhaseIdx, temperature_, p_[wPhaseIdx], updFluidState);
        Scalar v_g_ = 1. / FluidSystem::phaseDensity(nPhaseIdx, temperature_, p_[nPhaseIdx], updFluidState);
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
	updFluidState.update(Z1, pressure, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
	Scalar satt = updFluidState.saturation(wPhaseIdx);
	Scalar nuw = satt / v_w / (satt/v_w + (1-satt)/v_g);
    dv_dC1 = ((m1+m2) * (nuw * v_w + (1-nuw) * v_g) - volalt) /inc1;
    m1 -= inc1;

    // numerical derivative of fluid volume with respect to mass of component 2
    m2 += inc2;
    Z1 = m1 / (m1 + m2);
	updFluidState.update(Z1, pressure, problem_.spatialParameters().porosity(globalPos, *ep), temperature_);
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


}//end namespace Dumux
#endif
