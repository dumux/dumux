// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef DUMUX_FVPRESSURE2P_ADAPTIVE_HH
#define DUMUX_FVPRESSURE2P_ADAPTIVE_HH

// dumux environment
#include <dumux/decoupled/2p/2pproperties.hh>


/**
 * @file
 * @brief  Finite Volume Discretization of a pressure equation.
 * @author Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{

//! \ingroup FV2p
//! \brief Finite Volume discretization of the pressure equation of the sequential IMPES Model.
/*! Provides a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The definition of the total velocity \f$\boldsymbol{v}_total\f$ depends on the kind of pressure chosen. This could be a wetting (w) phase pressure leading to
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * a non-wetting (n) phase pressure yielding
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K}  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * or a global pressure leading to
 * \f[ - \text{div}\, \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q.\f]
 *  Here, \f$p\f$ denotes a pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation,\f$f\f$ the fractional flow function of a phase, \f$\rho\f$ a phase density, \f$g\f$ the gravity constant and \f$q\f$ the source term.
 * For all cases, \f$p = p_D\f$ on \f$\Gamma_{Neumann}\f$, and \f$\boldsymbol{v}_{total}  = q_N\f$
 * on \f$\Gamma_{Dirichlet}\f$.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2Padaptive
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

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
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;

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
    //! updates and stores constitutive relations
    void updateMaterialLaws();
    //! \copydoc Dumux::FVPressure1P::initialize()
    void initialize(bool solveTwice = true)
    {
        updateMaterialLaws();

        assemble(true);
        solve();
        if (solveTwice)
        {
            Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureOld(this->problem().variables().pressure());

            assemble(false);
            solve();

            Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureDiff(pressureOld);
            pressureDiff -= this->problem().variables().pressure();
            pressureOld = this->problem().variables().pressure();
            Scalar pressureNorm = pressureDiff.infinity_norm();
            pressureNorm /= pressureOld.infinity_norm();
            int numIter = 0;
            while (pressureNorm > 1e-5 && numIter < 10)
            {
                updateMaterialLaws();
                assemble(false);
                solve();

                pressureDiff = pressureOld;
                pressureDiff -= this->problem().variables().pressure();
                pressureNorm = pressureDiff.infinity_norm();
                pressureOld = this->problem().variables().pressure();
                pressureNorm /= pressureOld.infinity_norm();

                numIter++;
            }
            //            std::cout<<"Pressure defect = "<<pressureNorm<<"; "<<numIter<<" Iterations needed for initial pressure field"<<std::endl;
        }
        return;
    }
    //! \copydoc Dumux::FVPressure1P::pressure()
    void pressure(bool solveTwice = true)
    {
        assemble(false);
        solve();

        return;
    }
    //! updates the pressure field (analog to update function in Dumux::IMPET)
    void update()
    {
        updateMaterialLaws();

        pressure(false);

        return;
    }

    // serialization methods
    //! \copydoc Dumux::FVPressure1P::serialize(Restarter &res)
    template<class Restarter>
    void serialize(Restarter &res)
    {
        return;
    }

    //! \copydoc Dumux::FVPressure1P::deserialize(Restarter &res)
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \copydoc Dumux::FVPressure1P::addOutputVtkFields(MultiWriter &writer)
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        typename Variables::ScalarSolutionType *pressure = writer.allocateManagedBuffer (problem_.gridView().size(0));

        *pressure = problem_.variables().pressure();

        if (pressureType == pw)
        {
        writer.attachCellData(*pressure, "wetting pressure");
        }

        if (pressureType == pn)
        {
        writer.attachCellData(*pressure, "nonwetting pressure");
        }

        if (pressureType == pglobal)
        {
        writer.attachCellData(*pressure, "global pressure");
        }

        // output  phase-dependent stuff
        typename Variables::ScalarSolutionType *pC = writer.allocateManagedBuffer (problem_.gridView().size(0));
        *pC = problem_.variables().capillaryPressure();
        writer.attachCellData(*pC, "capillary pressure");

        typename Variables::ScalarSolutionType *densityWetting = writer.allocateManagedBuffer (problem_.gridView().size(0));
        *densityWetting = problem_.variables().densityWetting();
        writer.attachCellData(*densityWetting, "wetting density");

        typename Variables::ScalarSolutionType *densityNonwetting = writer.allocateManagedBuffer (problem_.gridView().size(0));
        *densityNonwetting = problem_.variables().densityNonwetting();
        writer.attachCellData(*densityNonwetting, "nonwetting density");

        typename Variables::ScalarSolutionType *viscosityWetting = writer.allocateManagedBuffer (problem_.gridView().size(0));
        *viscosityWetting = problem_.variables().viscosityWetting();
        writer.attachCellData(*viscosityWetting, "wetting viscosity");

        typename Variables::ScalarSolutionType *viscosityNonwetting = writer.allocateManagedBuffer (problem_.gridView().size(0));
        *viscosityNonwetting = problem_.variables().viscosityNonwetting();
        writer.attachCellData(*viscosityNonwetting, "nonwetting viscosity");

//        typename Variables::ScalarSolutionType *saturation = writer.allocateManagedBuffer (problem_.gridView().size(0));
//
//        *saturation = problem_.variables().saturation();
//
//        writer.attachCellData(*saturation, "wetting saturation");

        return;
    }

    //! Constructs a FVPressure2Padaptive object
    /**
     * \param problem a problem class object
     */
    FVPressure2Padaptive(Problem& problem) :
        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (4 * dim + 1)
                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()),
                gravity(problem.gravity())
//        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (4 * dim + 1)
//                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()),
//                gravity(problem.gravity())
    {// Comment: the term (4 * dim + 1) only works in the 2D case (or less). For 3D applications, it has to be (8*dim+1)
        if (pressureType != pw && pressureType != pn && pressureType != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType != Sw && saturationType != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        initializeMatrix();
    }

private:
    Problem& problem_;
    Matrix A_;
    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > f_;
protected:
    const GlobalPosition& gravity; //!< vector including the gravity constant
    static const bool compressibility = GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    static const int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

//initializes the matrix to store the system of equations
template<class TypeTag>
void FVPressure2Padaptive<TypeTag>::initializeMatrix()
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
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
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
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
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

//!function which assembles the system of equations to be solved. Only works for 2D.
template<class TypeTag>
void FVPressure2Padaptive<TypeTag>::assemble(bool first)
{
	updateMaterialLaws();
	// Adapt size of matrix A_ and Vector f_
	{int size = problem_.variables().gridSize();
	A_.setSize(size,size);
	f_.resize(size);}
	initializeMatrix(); // define non-zero entries in matrix

    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    int nodecount = 0;

    BoundaryTypes bcType;

    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().center();


        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

        // set sources
        PrimaryVariables source(0.0);
        problem().source(source, *eIt);
        if (!compressibility)
        {
            source[wPhaseIdx] /= densityWI;
            source[nPhaseIdx] /= densityNWI;
        }
        f_[globalIdxI] += volume * (source[wPhaseIdx] + source[nPhaseIdx]);

        Scalar porosity = problem_.spatialParameters().porosity(*eIt);

        // get mobilities and fractional flow factors
        Scalar lambdaWI = problem_.variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem_.variables().mobilityNonwetting(globalIdxI);
        Scalar fractionalWI = problem_.variables().fracFlowFuncWetting(globalIdxI);
        Scalar fractionalNWI = problem_.variables().fracFlowFuncNonwetting(globalIdxI);
        Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);

        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            int isIndex = isIt->indexInInside();

            // get normal vector
            const GlobalPosition& unitOuterNormal = isIt->centerUnitOuterNormal();

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                // "normal" case: neighbor has same level
                if (neighborPointer->level()==eIt.level())
                {
					// neighbor cell center in global coordinates
					const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

					// distance vector between barycenters
					GlobalPosition distVec = globalPosNeighbor - globalPos;

					// compute distance between cell centers
//					Scalar dist = distVec.two_norm();
					Scalar dist = distVec*unitOuterNormal;

					// compute vectorized permeabilities
					FieldMatrix meanPermeability(0);

					problem_.spatialParameters().meanK(meanPermeability,
							problem_.spatialParameters().intrinsicPermeability(*eIt),
							problem_.spatialParameters().intrinsicPermeability(*neighborPointer));

					Dune::FieldVector<Scalar, dim> permeability(0);
					meanPermeability.mv(unitOuterNormal, permeability);

					// get mobilities and fractional flow factors
					Scalar lambdaWJ = problem_.variables().mobilityWetting(globalIdxJ);
					Scalar lambdaNWJ = problem_.variables().mobilityNonwetting(globalIdxJ);
					Scalar fractionalWJ = problem_.variables().fracFlowFuncWetting(globalIdxJ);
					Scalar fractionalNWJ = problem_.variables().fracFlowFuncNonwetting(globalIdxJ);
					Scalar densityWJ = problem_.variables().densityWetting(globalIdxJ);
					Scalar densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);

					Scalar pcJ = problem_.variables().capillaryPressure(globalIdxJ);

					Scalar rhoMeanW = 0.5 * (densityWI + densityWJ);
					Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWJ);
					Scalar fMeanW = 0.5 * (fractionalWI + fractionalWJ);
					Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWJ);

					// update diagonal entry
					Scalar entry;

					//calculate potential gradients
					Scalar potentialW = 0;
					Scalar potentialNW = 0;

					Scalar densityW = 0;
					Scalar densityNW = 0;

					//if we are at the very first iteration we can't calculate phase potentials
					if (!first)
					{
						// After grid adaption, the potentials from the previous time step
						// can not be used for upwinding. The mean values is assumed instead.
						densityW = rhoMeanW;
						densityNW = rhoMeanNW;

						switch (pressureType)
						{
						case pw:
						{
							potentialW = (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ]);
							potentialNW = (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ] + pcI
									- pcJ);
							break;
						}
						case pn:
						{
							potentialW
									= (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ] - pcI + pcJ);
							potentialNW = (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ]);
							break;
						}
						case pglobal:
						{
							potentialW = (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ] - fMeanNW
									* (pcI - pcJ));
							potentialNW = (problem_.variables().pressure()[globalIdxI] - problem_.variables().pressure()[globalIdxJ] + fMeanW
									* (pcI - pcJ));
							break;
						}
						}

						potentialW += densityW * (distVec * gravity);
						potentialNW += densityNW * (distVec * gravity);

						//store potentials for further calculations (velocity, saturation, ...)
						problem_.variables().potentialWetting(globalIdxI, isIndex) = potentialW;
						problem_.variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;
					}

					//do the upwinding of the mobility depending on the phase potentials
					Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
					lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
					Scalar lambdaNW = (potentialNW > 0) ? lambdaNWI : lambdaNWJ;
					lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;

					densityW = (potentialW > 0.) ? densityWI : densityWJ;
					densityNW = (potentialNW > 0.) ? densityNWI : densityNWJ;

					densityW = (potentialW == 0) ? rhoMeanW : densityW;
					densityNW = (potentialNW == 0) ? rhoMeanNW : densityNW;

					//calculate current matrix entry
					entry = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

					//calculate right hand side
					Scalar rightEntry = (lambdaW * densityW + lambdaNW * densityNW) * (permeability * gravity) * faceArea;

					switch (pressureType)
					{
					case pw:
					{
						// calculate capillary pressure gradient
						Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
						pCGradient *= (pcI - pcJ) / dist;

						//add capillary pressure term to right hand side
						rightEntry += 0.5 * (lambdaNWI + lambdaNWJ) * (permeability * pCGradient) * faceArea;
						break;
					}
					case pn:
					{
						// calculate capillary pressure gradient
						Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
						pCGradient *= (pcI - pcJ) / dist;

						//add capillary pressure term to right hand side
						rightEntry -= 0.5 * (lambdaWI + lambdaWJ) * (permeability * pCGradient) * faceArea;
						break;
					}
					}

					//set right hand side
					f_[globalIdxI] -= rightEntry;

					// set diagonal entry
					A_[globalIdxI][globalIdxI] += entry;

					// set off-diagonal entry
					A_[globalIdxI][globalIdxJ] -= entry;

                } // end of case "neighbor has same level"

                // hanging node situation: neighbor has higher level
                if (neighborPointer->level()==eIt.level()+1)
                {
                	// Count number of hanging nodes
					// nodecount counts each hanging node twice, we have to divide by 2 in the end
					// not really necessary
					nodecount++;

					int globalIdxK = 0;
					ElementPointer thirdCellPointer = isIt->outside();
					bool foundK=false;
					bool foundIJ=false;
					// We are looking for two things:
					// IsIndexJ, the index of the interface from the neighbor-cell point of view
					// GlobalIdxK, the index of the third cell
					// for efficienty this is done in one IntersectionIterator-Loop

					// Intersectioniterator around cell J
					IntersectionIterator isItEndJ = this->problem().gridView().iend(*neighborPointer);

					for (IntersectionIterator isItJ = this->problem().gridView().ibegin(*neighborPointer); isItJ != isItEndJ; ++isItJ)
					{
						if (isItJ->neighbor())
						{
							ElementPointer neighborPointer2 = isItJ->outside();

							// Neighbor of neighbor is Cell I?
							if (this->problem().variables().index(*neighborPointer2)==globalIdxI)
							{
								foundIJ=true;
							}
							else
							{
								if (neighborPointer2->level()==eIt.level()+1)
								{
									// To verify, we found the correct Cell K, we check
									// - is level(K)=level(J)?
									// - is distance(IJ)=distance(IK)?
									// - K is neighbor of J.
									const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();
									const GlobalPosition& globalPosThirdCell = neighborPointer2->geometry().center();
									Dune::FieldVector<Scalar, dimWorld> distVecIJ = globalPosNeighbor - globalPos;
									Dune::FieldVector<Scalar, dimWorld> distVecIK = globalPosThirdCell - globalPos;
									if (((distVecIK.two_norm()-distVecIJ.two_norm()))/distVecIJ.two_norm() < 0.001)
									{
										globalIdxK= this->problem().variables().index(*neighborPointer2);
										thirdCellPointer = neighborPointer2;
										foundK=true;
									}

								}
							}
						}
						if (foundIJ && foundK) break;
					}

                    int isIndexJ = isIt->indexInOutside();

					// neighbor cell center in global coordinates
					const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();
//					const GlobalPosition& globalPosThirdCell = thirdCellPointer->geometry().center();
					const GlobalPosition& globalPosInterface = isIt->geometry().center();

					Dune::FieldVector<Scalar,dimWorld> distVec = globalPosInterface - globalPos;
					Scalar lI= distVec*unitOuterNormal;
					distVec = globalPosNeighbor - globalPosInterface;
					Scalar lJ= distVec*unitOuterNormal;
					Scalar l=lI+lJ;

					FieldMatrix permeabilityI(0);
					FieldMatrix permeabilityJ(0);
					FieldMatrix permeabilityK(0);

					problem_.spatialParameters().meanK(permeabilityI,
							problem_.spatialParameters().intrinsicPermeability(*eIt));
					problem_.spatialParameters().meanK(permeabilityJ,
							problem_.spatialParameters().intrinsicPermeability(*neighborPointer));
					problem_.spatialParameters().meanK(permeabilityK,
							problem_.spatialParameters().intrinsicPermeability(*thirdCellPointer));

					// Calculate permeablity component normal to interface
					Scalar kI, kJ, kK, kMean, ng;
					Dune::FieldVector<Scalar, dim> permI(0);
					Dune::FieldVector<Scalar, dim> permJ(0);
					Dune::FieldVector<Scalar, dim> permK(0);

					permeabilityI.mv(unitOuterNormal, permI);
					permeabilityJ.mv(unitOuterNormal, permJ);
					permeabilityK.mv(unitOuterNormal, permK);

					// kI,kJ,kK=(n^T)Kn
					kI=unitOuterNormal*permI;
					kJ=unitOuterNormal*permJ;
					kK=unitOuterNormal*permK;

					// See Diplomarbeit Michael Sinsbeck
					kMean=kI*kJ*kK*l/(kJ*kK*lI+kI*(kJ+kK)/2*lJ);

					ng=this->gravity*unitOuterNormal;

					// get mobilities
					// fractional flow function is not evaluated, since we do not use p_global
					Scalar lambdaWJ = problem_.variables().mobilityWetting(globalIdxJ);
					Scalar lambdaNWJ = problem_.variables().mobilityNonwetting(globalIdxJ);
					Scalar densityWJ = problem_.variables().densityWetting(globalIdxJ);
					Scalar densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);
					Scalar fractionalWJ = problem_.variables().fracFlowFuncWetting(globalIdxJ);
					Scalar fractionalNWJ = problem_.variables().fracFlowFuncNonwetting(globalIdxJ);

//					Scalar lambdaWK = problem_.variables().mobilityWetting(globalIdxK);
//					Scalar lambdaNWK = problem_.variables().mobilityNonwetting(globalIdxK);
					Scalar densityWK = problem_.variables().densityWetting(globalIdxK);
					Scalar densityNWK = problem_.variables().densityNonwetting(globalIdxK);
					Scalar fractionalWK = problem_.variables().fracFlowFuncWetting(globalIdxK);
					Scalar fractionalNWK = problem_.variables().fracFlowFuncNonwetting(globalIdxK);

					Scalar pcJ = problem_.variables().capillaryPressure(globalIdxJ);
					Scalar pcK = problem_.variables().capillaryPressure(globalIdxK);
					Scalar pcJK=(pcJ+pcK)/2;

					// Potentials from previous time step are not available.
					// Instead, calculate mean density, then find potentials,
					// then upwind density.
					// pressure from previous time step might also be incorrect.

					Scalar rhoMeanWIJ = (lJ*densityWI + lI*densityWJ)/l;
					Scalar rhoMeanNWIJ = (lJ*densityNWI + lI*densityNWJ)/l;
					Scalar rhoMeanWIK = (lJ*densityWI + lI*densityWK)/l;
					Scalar rhoMeanNWIK = (lJ*densityNWI + lI*densityNWK)/l;

					Scalar densityWIJ = 0;
					Scalar densityNWIJ = 0;
					Scalar densityWIK = 0;
					Scalar densityNWIK = 0;

					Scalar fMeanWIJ = (lJ*fractionalWI+lI*fractionalWJ)/l;
					Scalar fMeanNWIJ = (lJ*fractionalNWI+lI*fractionalNWJ)/l;
					Scalar fMeanWIK = (lJ*fractionalWI+lI*fractionalWK)/l;
					Scalar fMeanNWIK = (lJ*fractionalNWI+lI*fractionalNWK)/l;

					Scalar potentialWIJ = 0;
					Scalar potentialNWIJ = 0;
					Scalar potentialWIK = 0;
					Scalar potentialNWIK = 0;

					//if we are at the very first iteration we can't calculate phase potentials
					if (!first)
					{
						// potentials from previous time step no available.
						densityWIJ = rhoMeanWIJ;
						densityNWIJ = rhoMeanNWIJ;
						densityWIK = rhoMeanWIK;
						densityNWIK = rhoMeanNWIK;

						// Old pressure field
						Scalar pressI=problem_.variables().pressure()[globalIdxI];
						Scalar pressJ=problem_.variables().pressure()[globalIdxJ];
						Scalar pressK=problem_.variables().pressure()[globalIdxK];
						Scalar pressJK=(pressJ+pressK)/2;
						Scalar pcJK=(pcJ+pcK)/2;

						switch (pressureType)
						{
							case pw:
							{
								potentialWIJ = (pressI-pressJK)/l+
										(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
								potentialNWIJ = (pressI+pcI-(pressJK+pcJK))/l+
										(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
								potentialWIK = (pressI-pressJK)/l+
										(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
								potentialNWIK = (pressI+pcI-(pressJK+pcJK))/l+
										(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
								break;
							}
							case pn:
							{
								potentialWIJ = (pressI-pcI-(pressJK-pcJK))/l+
										(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
								potentialNWIJ = (pressI-pressJK)/l+
										(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
								potentialWIK = (pressI-pcI-(pressJK-pcJK))/l+
										(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
								potentialNWIK = (pressI-pressJK)/l+
										(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
								break;
							}
							case pglobal:
							{
								potentialWIJ = (pressI-fMeanNWIJ*pcI-(pressJK-fMeanNWIJ*pcJK))/l+
										(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
								potentialNWIJ = (pressI+fMeanWIJ*pcI-(pressJK+fMeanWIJ*pcJK))/l+
										(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
								potentialWIK = (pressI-fMeanNWIK*pcI-(pressJK-fMeanNWIK*pcJK))/l+
										(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
								potentialNWIK = (pressI+fMeanWIK*pcI-(pressJK+fMeanWIK*pcJK))/l+
										(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
								break;
							}
						}

						//store potentials for further calculations (velocity, saturation, ...)

						problem_.variables().potentialWetting(globalIdxI, isIndex) = potentialWIJ;
						problem_.variables().potentialNonwetting(globalIdxI, isIndex) = potentialNWIJ;
						problem_.variables().potentialWetting(globalIdxJ, isIndexJ) = -potentialWIJ;
						problem_.variables().potentialNonwetting(globalIdxJ, isIndexJ) = -potentialNWIJ;
					}

					//do the upwinding of the mobility depending on the phase potentials
					Scalar lambdaWIJ = (potentialWIJ > 0.) ? lambdaWI : lambdaWJ;
					lambdaWIJ = (potentialWIJ == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
					Scalar lambdaNWIJ = (potentialNWIJ > 0) ? lambdaNWI : lambdaNWJ;
					lambdaNWIJ = (potentialNWIJ == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNWIJ;
					densityWIJ = (potentialWIJ > 0.) ? densityWI : densityWJ;
					densityNWIJ = (potentialNWIJ > 0.) ? densityNWI : densityNWJ;
					densityWIJ = (potentialWIJ == 0) ? rhoMeanWIJ : densityWIJ;
					densityNWIJ = (potentialNWIJ == 0) ? rhoMeanNWIJ : densityNWIJ;

					densityWIK = (potentialWIK > 0.) ? densityWI : densityWK;
					densityNWIK = (potentialNWIK > 0.) ? densityNWI : densityNWK;
					densityWIK = (potentialWIK == 0) ? rhoMeanWIK : densityWIK;
					densityNWIK = (potentialNWIK == 0) ? rhoMeanNWIK : densityNWIK;


					// update diagonal entry and right hand side

					Scalar entry;
					Scalar rightEntry;

					entry = (lambdaWIJ + lambdaNWIJ)*kMean/l*faceArea;
					rightEntry  = faceArea*lambdaWIJ*kMean*ng*((densityWIJ)-(lJ/l)*(kI+kK)/kI*(densityWIK-densityWIJ)/2);
					rightEntry += faceArea*lambdaNWIJ*kMean*ng*((densityNWIJ)-(lJ/l)*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2);

					switch (pressureType)
					{
						case pw:
						{
							rightEntry += faceArea*lambdaNWIJ*kMean*(pcJK-pcI)/l;
							break;
						}
						case pn:
						{
							rightEntry -= faceArea*lambdaWIJ*kMean*(pcJK-pcI)/l;
							break;
						}
					}

					f_[globalIdxI] -= rightEntry;
					f_[globalIdxJ] += rightEntry;

					// set diagonal entry
					A_[globalIdxI][globalIdxI] += entry;
					A_[globalIdxI][globalIdxJ] -= 0.5*entry;
					A_[globalIdxI][globalIdxK] -= 0.5*entry;

					// set off-diagonal entry
					A_[globalIdxJ][globalIdxI] -= entry;
					A_[globalIdxJ][globalIdxJ] += 0.5*entry;
					A_[globalIdxJ][globalIdxK] += 0.5*entry;

                }
            }

            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                problem().boundaryTypes(bcType, *isIt);
                PrimaryVariables boundValues(0.0);

                GlobalPosition distVec(globalPosFace - globalPos);
//                Scalar dist = distVec.two_norm();
                Scalar dist = distVec * unitOuterNormal;

                if (bcType.isDirichlet(eqIdxPress))
                {
                    problem().dirichlet(boundValues, *isIt);

                    //permeability vector at boundary
                    // compute vectorized permeabilities
                    FieldMatrix meanPermeability(0);

                    problem_.spatialParameters().meanK(meanPermeability,
                            problem_.spatialParameters().intrinsicPermeability(*eIt));

                    Dune::FieldVector<Scalar, dim> permeability(0);
                    meanPermeability.mv(unitOuterNormal, permeability);

                    //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                    Scalar satBound;
                    if (bcType.isDirichlet(eqIdxSat))
                    {
                        satBound = boundValues[saturationIdx];
                    }
                    else
                    {
                        satBound = problem_.variables().saturation()[globalIdxI];
                    }
                    Scalar temperature = problem_.temperature(*eIt);

                    //get dirichlet pressure boundary condition
                    Scalar pressBound = boundValues[pressureIdx];

                    //calculate consitutive relations depending on the kind of saturation used
                    //determine phase saturations from primary saturation variable
                    Scalar satW = 0;
                    //Scalar satNW = 0;
                    switch (saturationType)
                    {
                    case Sw:
                    {
                        satW = satBound;
                        //satNW = 1 - satBound;
                        break;
                    }
                    case Sn:
                    {
                        satW = 1 - satBound;
                        //satNW = satBound;
                        break;
                    }
                    }

                    Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);
                    Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

                    //determine phase pressures from primary pressure variable
                    Scalar pressW = 0;
                    Scalar pressNW = 0;
                    switch (pressureType)
                    {
                    case pw:
                    {
                        pressW = pressBound;
                        pressNW = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressW = pressBound - pcBound;
                        pressNW = pressBound;
                        break;
                    }
                    }

                    Scalar densityWBound = 0;
                    Scalar densityNWBound = 0;
                    Scalar lambdaWBound = 0;
                    Scalar lambdaNWBound = 0;

                    if (compressibility)
                    {
                        FluidState fluidState;
                        fluidState.update(satW, pressW, pressNW, temperature);
                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, pressW, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, pressNW, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressNW, fluidState);
                        lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosityWBound * densityWBound;
                        lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosityNWBound * densityNWBound;
                    }
                    else
                    {
                        Scalar referencePressure = problem_.referencePressure(*eIt);
                        FluidState fluidState;
                        fluidState.update(satW, referencePressure, referencePressure, temperature);

                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
                        lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosityWBound;
                        lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosityNWBound;
                    }
                    Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNWBound);
                    Scalar fractionalNWBound = lambdaNWBound / (lambdaWBound + lambdaNWBound);

                    Scalar rhoMeanW = 0.5 * (densityWI + densityWBound);
                    Scalar rhoMeanNW = 0.5 * (densityNWI + densityNWBound);
                    Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
                    Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWBound);

                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    Scalar densityW = 0;
                    Scalar densityNW = 0;

                    if (!first)
                    {
                    	// Comment: In the non-adaptive case, one can upwind the density using the
                    	// old potential. Here, we use the mean density instead.
                    	// For the incompressible case, it does not make a difference, anyway.
						densityW = rhoMeanW;
						densityNW = rhoMeanNW;

                        //calculate potential gradient
                        switch (pressureType)
                        {
                        case pw:
                        {
                            potentialW = (problem_.variables().pressure()[globalIdxI] - pressBound);
                            potentialNW = (problem_.variables().pressure()[globalIdxI] + pcI - pressBound - pcBound);
                            break;
                        }
                        case pn:
                        {
                            potentialW = (problem_.variables().pressure()[globalIdxI] - pcI - pressBound + pcBound);
                            potentialNW = (problem_.variables().pressure()[globalIdxI] - pressBound);
                            break;
                        }
                        case pglobal:
                        {
                            potentialW = (problem_.variables().pressure()[globalIdxI] - pressBound - fMeanNW * (pcI - pcBound));
                            potentialNW = (problem_.variables().pressure()[globalIdxI] - pressBound + fMeanW * (pcI - pcBound));
                            break;
                        }
                        }

                        potentialW += densityW * (distVec * gravity);
                        potentialNW += densityNW * (distVec * gravity);

                        //store potential gradients for further calculations
                        problem_.variables().potentialWetting(globalIdxI, isIndex) = potentialW;
                        problem_.variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;
                    }

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
                    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
                    Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
                    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;

                    densityW = (potentialW > 0.) ? densityWI : densityWBound;
                    densityW = (potentialW == 0) ? rhoMeanW : densityW;
                    densityNW = (potentialNW > 0.) ? densityNWI : densityNWBound;
                    densityNW = (potentialNW == 0) ? rhoMeanNW : densityNW;

                    //calculate current matrix entry
                    Scalar entry = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

                    //calculate right hand side
                    Scalar rightEntry = (lambdaW * densityW + lambdaNW * densityNW) * (permeability * gravity) * faceArea;

                    switch (pressureType)
                    {
                    case pw:
                    {
                        // calculate capillary pressure gradient
                        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                        pCGradient *= (pcI - pcBound) / dist;

                        //add capillary pressure term to right hand side
                        rightEntry += 0.5 * (lambdaNWI + lambdaNWBound) * (permeability * pCGradient) * faceArea;
                        break;
                    }
                    case pn:
                    {
                        // calculate capillary pressure gradient
                        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                        pCGradient *= (pcI - pcBound) / dist;

                        //add capillary pressure term to right hand side
                        rightEntry -= 0.5 * (lambdaWI + lambdaWBound) * (permeability * pCGradient) * faceArea;
                        break;
                    }
                    }

                    // set diagonal entry and right hand side entry
                    A_[globalIdxI][globalIdxI] += entry;
                    f_[globalIdxI] += entry * pressBound;
                    f_[globalIdxI] -= rightEntry;
                }
                //set neumann boundary condition

                else if (bcType.isNeumann(eqIdxPress))
                {
                    problem().neumann(boundValues, *isIt);

                    if (!compressibility)
                    {
                        boundValues[wPhaseIdx] /= densityWI;
                        boundValues[nPhaseIdx] /= densityNWI;
                    }
                    f_[globalIdxI] -= (boundValues[wPhaseIdx] + boundValues[nPhaseIdx]) * faceArea;

                    //Assumes that the phases flow in the same direction at the neumann boundary, which is the direction of the total flux!!!
                    //needed to determine the upwind direction in the saturation equation
                    problem_.variables().potentialWetting(globalIdxI, isIndex) = boundValues[wPhaseIdx];
                    problem_.variables().potentialNonwetting(globalIdxI, isIndex) = boundValues[nPhaseIdx];
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
                }
            }
        } // end all intersections

        //volume correction due to density differences
        if (compressibility)
        {
            switch (saturationType)
            {
            case Sw:
            {
                f_[globalIdxI] -= problem_.variables().volumecorrection(globalIdxI) * porosity * volume * (densityWI - densityNWI);
                break;
            }
            case Sn:
            {
                f_[globalIdxI] -= problem_.variables().volumecorrection(globalIdxI) * porosity * volume * (densityNWI - densityWI);
                break;
            }
            }
        }

//        Dune::printmatrix(std::cout,A_,"title","-");
    } // end grid traversal
    return;
}

//!solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure2Padaptive<TypeTag>::solve()
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LinearSolver)) Solver;

    int verboseLevelSolver = GET_PARAM(TypeTag, int, LinearSolver, Verbosity);

    if (verboseLevelSolver)
    std::cout << "FVPressure2Padaptive: solve for pressure" << std::endl;

//    printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//    printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);

    Solver solver(problem_);
    solver.solve(A_, problem_.variables().pressure(), f_);

//                    printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);

    //for parallel use
//    problem_.variables().communicatePressure();

    return;
}
//!constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class TypeTag>
void FVPressure2Padaptive<TypeTag>::updateMaterialLaws()
{
//    //for parallel use
////    printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);
////    problem_.variables().communicatePressure();
////    printvector(std::cout, (problem_.variables().pressure()), "pressureComm", "row", 200, 1, 3);
//
////    printvector(std::cout, (problem_.variables().saturation()), "sat", "row", 200, 1, 3);
//    problem_.variables().communicateTransportedQuantity();
////    printvector(std::cout, (problem_.variables().saturation()), "satComm", "row", 200, 1, 3);

    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        Scalar temperature = problem_.temperature(*eIt);

        //determine phase saturations from primary saturation variable
        Scalar satW = 0;
        //Scalar satNW = 0;
        switch (saturationType)
        {
        case Sw:
        {
            satW = problem_.variables().saturation()[globalIdx];
            //satNW = 1 - problem_.variables().saturation()[globalIdx];
            break;
        }
        case Sn:
        {
            satW = 1 - problem_.variables().saturation()[globalIdx];
            //satNW = problem_.variables().saturation()[globalIdx];
            break;
        }
        }

        problem_.variables().capillaryPressure(globalIdx) = MaterialLaw::pC(
            problem_.spatialParameters().materialLawParams(*eIt), satW);

        //determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNW = 0;
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

        Scalar densityW = 0;
        Scalar densityNW = 0;
        Scalar viscosityW = 0;
        Scalar viscosityNW = 0;

        if (compressibility)
        {
            fluidState.update(satW, pressW, pressNW, temperature);
        }
        else
        {
            pressW = problem_.referencePressure(*eIt);
            pressNW = problem_.referencePressure(*eIt);
            fluidState.update(satW, pressW, pressNW, temperature);
        }

        densityW = FluidSystem::phaseDensity(wPhaseIdx, temperature, pressW, fluidState);
        densityNW = FluidSystem::phaseDensity(nPhaseIdx, temperature, pressNW, fluidState);

        viscosityW = FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
        viscosityNW = FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressNW, fluidState);

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityW;
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityNW;

        if (compressibility)
        {
            mobilityW *= densityW;
            mobilityNW *= densityNW;
        }

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx) = mobilityW;
        problem_.variables().mobilityNonwetting(globalIdx) = mobilityNW;

        // initialize densities
        problem_.variables().densityWetting(globalIdx) = densityW;
        problem_.variables().densityNonwetting(globalIdx) = densityNW;

        // initialize viscosities
        problem_.variables().viscosityWetting(globalIdx) = viscosityW;
        problem_.variables().viscosityNonwetting(globalIdx) = viscosityNW;

        //initialize fractional flow functions
        problem_.variables().fracFlowFuncWetting(globalIdx) = mobilityW / (mobilityW + mobilityNW);
        problem_.variables().fracFlowFuncNonwetting(globalIdx) = mobilityNW / (mobilityW + mobilityNW);
    }
    return;
}

}
#endif
