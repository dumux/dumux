// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
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
/*!
 * \file
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */
#ifndef DUMUX_2PDFM_VOLUME_VARIABLES_HH
#define DUMUX_2PDFM_VOLUME_VARIABLES_HH

#include "2pdfmproperties.hh"

//#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/boxmodels/2p/2pvolumevariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPDFMVolumeVariables : public TwoPVolumeVariables <TypeTag>
//	public BoxVolumeVariables<TypeTag>
{
    typedef TwoPVolumeVariables<TypeTag> ParentType;
//    typedef BoxVolumeVariables<TypeTag> ParentType;
public:

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pwSn = Indices::pwSn,
        pnSw = Indices::pnSw,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::UGGrid<2> GridType; //TODO
    typedef typename GridType::ctype DT;

    enum {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
        typedef Dune::FieldVector<Scalar, dim> LocalPosition;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The local primary variable vector
     * \param problem The problem object
     * \param element The current element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        completeFluidState(priVars, problem, element, elemGeom, scvIdx, fluidState_);

        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(wPhaseIdx);

        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(nPhaseIdx);

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(priVars, problem, element, elemGeom, scvIdx, isOldSol);
        asImp_().updateFracture(priVars, problem, element, elemGeom, scvIdx, isOldSol);
    }

    void updateFracture(const PrimaryVariables &priVars,
						const Problem &problem,
						const Element &element,
						const FVElementGeometry &elemGeom,
						int scvIdx,
						bool isOldSol)
    {
    	PrimaryVariables varsFracture;
    	 const MaterialLawParams &materialParamsMatrix =
    	            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);
		Scalar pressure[numPhases];
		Scalar pMatrix[numPhases];
		Scalar pFract[numPhases];

		// coordinates of the vertex
		const GlobalPosition &global = element.geometry().corner(scvIdx);

		satNMatrix_  = priVars[saturationIdx];
		satWMatrix_  = 1.0 - satNMatrix_;
		satN_ = satNMatrix_;
		satW_ = satWMatrix_;

//        Valgrind::CheckDefined(satWMatrix);
//        Valgrind::CheckDefined(materialParamsMatrix);

		pCMatrix_ = MaterialLaw::pC(materialParamsMatrix, satWMatrix_);
		if (problem.one_phase_problem){
			pCMatrix_ = 0.0;
		}
		pC_ = pCMatrix_;
		//pressures
		pMatrix[wPhaseIdx] = priVars[pressureIdx];
		pMatrix[nPhaseIdx] = pMatrix[wPhaseIdx] + pCMatrix_;
		//Initialize pFract with the same values as the ones in the matrix
		pFract[wPhaseIdx] = pMatrix[wPhaseIdx];
		pFract[nPhaseIdx] = satNMatrix_;

		varsFracture[pressureIdx] = pFract[wPhaseIdx];
		varsFracture[saturationIdx] = pFract[wPhaseIdx];
		completeFluidState(varsFracture, problem, element, elemGeom, scvIdx, fluidStateFracture_);

		//Checks if the node is on a fracture
		isNodeOnFracture_ = problem.spatialParameters().isVertexFracture(element, scvIdx);

		///////////////////////////////////////////////////////////////////////////////
		if (isNodeOnFracture_){
			const MaterialLawParams &materialParamsFracture =
					problem.spatialParameters().materialLawParamsFracture(element, elemGeom, scvIdx);

			satNFracture_ = priVars[saturationIdx];
			satWFracture_ = 1 - satNFracture_;
			pCFracture_ = MaterialLaw::pC(materialParamsFracture, satWFracture_);
			if (problem.one_phase_problem){
				pCFracture_ = 0.0;
			}
			pFract[wPhaseIdx] = priVars[pressureIdx];
			pFract[nPhaseIdx] = pFract[wPhaseIdx] + pCFracture_;
			pEntryMatrix_ = MaterialLaw::pC(materialParamsMatrix, 1);

			//use interface condition - extended capillary pressure inteface condition
			if (problem.use_interface_condition){
				interfaceCondition(materialParamsMatrix);
			}
			pC_ = pCFracture_;
			satW_ = satWFracture_; //for plotting we are interested in the saturations of the fracture
			satN_ = satNFracture_;
			mobilityFracture_[wPhaseIdx] =
		            MaterialLaw::krw(materialParamsFracture, fluidStateFracture_.saturation(wPhaseIdx))
		            / fluidStateFracture_.viscosity(wPhaseIdx);

			mobilityFracture_[nPhaseIdx] =
					MaterialLaw::krn(materialParamsFracture, fluidStateFracture_.saturation(wPhaseIdx))
	            		/ fluidStateFracture_.viscosity(nPhaseIdx);

			dSM_dSF_ = (1 - problem.spatialParameters().SwrM_) / (1 - problem.spatialParameters().SwrF_)
					* pow((problem.spatialParameters().pdM_/ problem.spatialParameters().pdF_),problem.spatialParameters().lambdaM_)
					* (problem.spatialParameters().lambdaM_ / problem.spatialParameters().lambdaF_)
					* pow((satWFracture_ - problem.spatialParameters().SwrF_ ) / (1 - problem.spatialParameters().SwrF_),
							(problem.spatialParameters().lambdaM_ / problem.spatialParameters().lambdaF_) - 1);
		}// end if (node)
		///////////////////////////////////////////////////////////////////////////////
		else {
			satNFracture_ = -1;
			satWFracture_ = -1;
			pCFracture_ = -1e100;
			pFract[wPhaseIdx] = -1e100;
			pFract[nPhaseIdx] = -1e100;
			pEntryMatrix_ = -1e100;
			mobilityFracture_[wPhaseIdx] = 0.0;
			mobilityFracture_[nPhaseIdx] = 0.0;
		}
		///////////////////////////////////////////////////////////////////////////////
		pressure[wPhaseIdx] = priVars[pressureIdx];
		pressure[nPhaseIdx] = pressure[wPhaseIdx] + pC_;


//		fluidState_.update(satN_, pressure[wPhaseIdx], pressure[nPhaseIdx], temperature());
//
//		mobility_[wPhaseIdx] =
//			 MaterialLaw::krw(materialParamsMatrix, 1 - satNMatrix_)
//			 /
//			 FluidSystem::phaseViscosity(wPhaseIdx,
//					 	 	 	 	 	 temperature(),
//										 pMatrix[wPhaseIdx],
//										 fluidState_);
//		mobility_[nPhaseIdx] =
//			MaterialLaw::krn(materialParamsMatrix, 1 - satNMatrix_)
//			/
//			FluidSystem::phaseViscosity(nPhaseIdx,
//										temperature(),
//										pMatrix[nPhaseIdx],
//										fluidState_);

		porosityFracture_ = problem.spatialParameters().porosityFracture(element,
																 elemGeom,
																 scvIdx);


		Valgrind::CheckDefined(mobilityFracture_[wPhaseIdx]);
		Valgrind::CheckDefined(mobilityFracture_[nPhaseIdx]);

		Valgrind::CheckDefined(fluidState_);
		Valgrind::CheckDefined(porosity_);
		Valgrind::CheckDefined(porosityFracture_);
		Valgrind::CheckDefined(temperature_);

		Valgrind::CheckDefined(mobility_);
		Valgrind::CheckDefined(mobilityFracture_);

		Valgrind::CheckDefined(satW);
		Valgrind::CheckDefined(satWFracture);
		Valgrind::CheckDefined(satN);
		Valgrind::CheckDefined(satNFracture);
		Valgrind::CheckDefined(pC);
		Valgrind::CheckDefined(pCFracture);
		Valgrind::CheckDefined(pEntryMatrix);

		Valgrind::CheckDefined(isNodeOnFracture);
    }


    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& primaryVariables,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& elementGeometry,
                                   int scvIdx,
                                   FluidState& fluidState)
    {
        Scalar t = Implementation::temperature_(primaryVariables, problem, element,
                                                elementGeometry, scvIdx);
        fluidState.setTemperature(t);

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const typename MaterialLaw::Params &materialParams =
            problem.spatialParameters().materialLawParams(element, elementGeometry, scvIdx);


        if (int(formulation) == pwSn) {
            Scalar Sn = primaryVariables[saturationIdx];
            fluidState.setSaturation(nPhaseIdx, Sn);
            fluidState.setSaturation(wPhaseIdx, 1 - Sn);

            Scalar pW = primaryVariables[pressureIdx];
            fluidState.setPressure(wPhaseIdx, pW);
            fluidState.setPressure(nPhaseIdx,
                                   pW + MaterialLaw::pC(materialParams, 1 - Sn));
        }
        else if (int(formulation) == pnSw) {
            Scalar Sw = primaryVariables[saturationIdx];
            fluidState.setSaturation(wPhaseIdx, Sw);
            fluidState.setSaturation(nPhaseIdx, 1 - Sw);

            Scalar pN = primaryVariables[pressureIdx];
            fluidState.setPressure(nPhaseIdx, pN);
            fluidState.setPressure(wPhaseIdx,
                                   pN - MaterialLaw::pC(materialParams, Sw));
        }

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // compute and set the enthalpy
            Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }


    /*
     * \brief Extended capillary pressure saturation interface condition //TODO implement here
     */
     void interfaceCondition(const MaterialLawParams &materialParamsMatrix)
     {
         /*2nd condition Niessner, J., R. Helmig, H. Jakobs, and J.E. Roberts. 2005, eq.10
          * if the capillary pressure in the fracture is smaller than the entry pressure
          * in the matrix than in the matrix
          * */
         	if (pCFracture_ <= pEntryMatrix_)
         	{
         		satWMatrix_ = 1.0; ///TODO
         		satNMatrix_ = 1 - satWMatrix_;
         	}
         	//3rd condition Niessner, J., R. Helmig, H. Jakobs, and J.E. Roberts. 2005, eq.10
         	else
         	{
         		/*
         		 * Inverse capillary pressure function SwM = pcM^(-1)(pcF(SwF))
         		 */
         		satWMatrix_ = MaterialLaw::Sw(materialParamsMatrix, pCFracture_);
         		satNMatrix_ = 1 - satWMatrix_;
         	}
     }


    /*
     * \brief Calculates the volume of the fracture inside the SCV
     */
    Scalar calculateSCVFractureVolume ( const Problem &problem,
    		const Element &element,
            const FVElementGeometry &elemGeom,
            int scvIdx)
    {
    	Scalar volSCVFracture;
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::GenericReferenceElementContainer<DT,dim>::value_type&
        	refElem = Dune::GenericReferenceElements<DT,dim>::general(gt);

        for (int faceIdx=0; faceIdx<refElem.size(1); faceIdx++)
        {
    		SCVFace face = elemGeom.subContVolFace[faceIdx];
        	int i=face.i;
			int j=face.j;

        	if (problem.spatialParameters().isEdgeFracture(element, faceIdx)
        			&& (i==scvIdx || j==scvIdx))
        	{
                Scalar fracture_width = problem.spatialParameters().fractureWidth();

                const GlobalPosition global_i = element.geometry().corner(i);
                const GlobalPosition global_j = element.geometry().corner(j);
				GlobalPosition diff_ij = global_j;
				diff_ij -=global_i;
				//fracture length in the subcontrol volume is half d_ij
				Scalar fracture_length = 0.5*diff_ij.two_norm();

//    				std::cout<<"I["<<global<<"] fracture_length"<<fracture_length <<"\n"; //TODO delete
//    				std::cout<<"fracture_width"<<fracture_width<<"\n"; //TODO delete
//    				//half is taken when looping on the other element
				volSCVFracture += 0.5 * fracture_length * fracture_width;
        	}
        }
        return volSCVFracture;
    }



    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    Scalar saturationFracture(int phaseIdx) const
    {
   	 if (phaseIdx == wPhaseIdx)
   		 return satWFracture_;
   	 else
   		 return satNFracture_;
    }
    Scalar saturationMatrix(int phaseIdx) const
    {
		 if (phaseIdx == wPhaseIdx)
			 return satWMatrix_;
		 else
			 return satNMatrix_;
    }


    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume [Pa].
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }


    Scalar mobilityFracture(int phaseIdx) const
         { return mobilityFracture_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }


    Scalar porosityFracture() const
    { return porosityFracture_; }

    /*
     *
     */
    Scalar dSM_dSF() const
    {
   	return dSM_dSF_;
    }


protected:
    friend class TwoPVolumeVariables<TypeTag>;
    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &elemGeom,
                               int scvIdx)
    {
    	return problem.boxTemperature(element, elemGeom, scvIdx);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return 0;
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int vertIdx,
                       bool isOldSol)
    { }

    FluidState fluidState_;
    FluidState fluidStateFracture_;
    Scalar porosity_;
    Scalar porosityFracture_;
    Scalar mobility_[numPhases];
    Scalar mobilityFracture_[numPhases];

    Scalar satW_;
    Scalar satWFracture_;
    Scalar satWMatrix_;
    Scalar satN_;
    Scalar satNFracture_;
    Scalar satNMatrix_;

    Scalar pC_;
    Scalar pCFracture_;
    Scalar pCMatrix_;
    Scalar pEntryMatrix_;
    Scalar dSM_dSF_;

    bool isNodeOnFracture_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
