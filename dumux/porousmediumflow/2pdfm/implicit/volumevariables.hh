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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase discrete fracture-matrix model.
 */
#ifndef DUMUX_MODELS_2PDFM_VOLUME_VARIABLES_HH
#define DUMUX_MODELS_2PDFM_VOLUME_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/porousmediumflow/2p/implicit/volumevariables.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPDFMModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase discrete fracture-matrix model.
 */
template <class TypeTag>
class TwoPDFMVolumeVariables : public TwoPVolumeVariables<TypeTag>
{
    typedef TwoPVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pwsn = Indices::pwsn,
        pnsw = Indices::pnsw,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
    typedef typename GridType::ctype DT;

    enum {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Dune::ReferenceElements<DT, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<DT, dim> ReferenceElement;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename Dumux::VertIdxToMinPcMapper<TypeTag> VertIdxToMinPcMapper;

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        // initially we set the matrix saturations to the actual solution
        satNMatrix_  = priVars[saturationIdx];
        satWMatrix_  = 1.0 - satNMatrix_;

        // fracture variables are set to zero for the case that no fracture is present in the scv
        satNFracture_ = 0.0;
        satWFracture_ = 0.0;
        porosityFracture_ = 0.0;
        permeabilityFracture_ = 0.0;
        mobilityFracture_[wPhaseIdx] = 0;
        mobilityFracture_[nPhaseIdx] = 0;

        // calculate the matrix mobilities
        this->completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);
        const MaterialLawParams &materialParams = problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        mobilityMatrix_[wPhaseIdx] = MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx)) / fluidState_.viscosity(wPhaseIdx);
        mobilityMatrix_[nPhaseIdx] = MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx)) / fluidState_.viscosity(nPhaseIdx);

        // update the fracture if fracture is present in the scv
        if (problem.spatialParams().isVertexFracture(element, scvIdx))
            asImp_().updateFracture(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \brief Construct the volume variables for all fracture vertices.
     *
     * \param priVars Primary variables.
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx Sub-control volume index
     * \param isOldSol Tells whether the model's previous or current solution should be used.
     */
    void updateFracture(const PrimaryVariables &priVars,
                        const Problem &problem,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        int scvIdx,
                        bool isOldSol)
    {
        // the scv is on a fracture, the actual solution corresponds to the solution in the fracture
        // the matrix saturation will be modified below if interface condition is used
        satNFracture_ = priVars[saturationIdx];
        satWFracture_ = 1.0 - satNFracture_;

        //update fracture saturation using interface condition
        const VertexMapper &vertexMapper = problem.vertexMapper();
        int globalIdx = vertexMapper.subIndex(element, scvIdx, dim);
        const MaterialLawParams &materialParamsFracture = problem.spatialParams().materialLawParamsFracture(element, fvGeometry, scvIdx);

        // this mapper provides access to the minimum pc for neighboring elements
        const VertIdxToMinPcMapper &vertIdxToMinPcMapper = problem.vertIdxToMinPcMapper();

        //calculates capillary pressure and entry pressure for current scv
        Scalar pc = MaterialLaw::pc(materialParamsFracture, satWFracture_);
        Scalar pe = MaterialLaw::pc(materialParamsFracture, 1-materialParamsFracture.snr());

        FVElementGeometry minPcElemFvGeometry;
        Element minPcElem = vertIdxToMinPcMapper.vertexElement(globalIdx);
        minPcElemFvGeometry.update(problem.gridView(), minPcElem);
        //choose materialLawParamsFracture as fracture has higher pe
        MaterialLawParams minPcElemMaterialParams = problem.spatialParams().materialLawParamsFracture(minPcElem, minPcElemFvGeometry, scvIdx);
        //calculate capillary pressure based on the pc-sw relation of the minPc element
        Scalar pcmin = MaterialLaw::pc(minPcElemMaterialParams, satWFracture_);
        // update fracture saturation
        if (std::abs(pc-pcmin) < 1e-6){
            }
        else if (pcmin < pe){
            satNFracture_ = std::min(materialParamsFracture.snr(), 0.0 /* SnInitial*/);
        }
        else{
            satNFracture_ = 1 - MaterialLaw::sw(materialParamsFracture, pcmin);
        }
        satWFracture_ = 1.0 - satNFracture_;
        //primary variables for update of fluid state for fracture
        PrimaryVariables priVarsFracture = priVars;
       if (int(formulation) == pwsn) {
            priVarsFracture[saturationIdx]= satNFracture_;
            }
       else if (int(formulation) == pnsw) {
            priVarsFracture[saturationIdx]= satWFracture_;}

        const MaterialLawParams &materialParamsMatrix =
                    problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        // the fluid state of the fracture
        this->completeFluidState(priVarsFracture, problem, element, fvGeometry, scvIdx, fluidStateFracture_);


        // use interface condition - extended capillary pressure inteface condition
        if (problem.useInterfaceCondition())
        {
            std::cout << " old satM: "<< satNMatrix_ << std::endl;
            // updated solution in the fracture is used for the interface condition
            interfaceCondition(priVarsFracture, materialParamsMatrix, materialParamsFracture);
            std::cout << " new satM: "<< satNMatrix_ << std::endl;
            // After modifying the matrix saturations we have to update the fluid state and the matrix mobilities
            PrimaryVariables updatedMatrixPriVars(priVarsFracture);
            updatedMatrixPriVars[saturationIdx] = satNMatrix_;
            this->completeFluidState(updatedMatrixPriVars, problem, element, fvGeometry, scvIdx, fluidState_);

            mobilityMatrix_[wPhaseIdx] = MaterialLaw::krw(materialParamsMatrix, satWMatrix_) / fluidState_.viscosity(wPhaseIdx);
            mobilityMatrix_[nPhaseIdx] = MaterialLaw::krn(materialParamsMatrix, satWMatrix_) / fluidState_.viscosity(nPhaseIdx);
        }

        // calculate the mobilities in the fracture using the material parameters in the fracture
        mobilityFracture_[wPhaseIdx] =
                MaterialLaw::krw(materialParamsFracture, fluidStateFracture_.saturation(wPhaseIdx))
                / fluidStateFracture_.viscosity(wPhaseIdx);

        mobilityFracture_[nPhaseIdx] =
                MaterialLaw::krn(materialParamsFracture, fluidStateFracture_.saturation(wPhaseIdx))
                    / fluidStateFracture_.viscosity(nPhaseIdx);

        // set the porosity and permeability
        porosityFracture_ = problem.spatialParams().porosityFracture(element, fvGeometry, scvIdx);
        permeabilityFracture_ = problem.spatialParams().intrinsicPermeabilityFracture(element, fvGeometry, scvIdx);
    }

    /*!
     * \brief Extended capillary pressure saturation interface condition
     *
     * \param priVars Primary variables
     * \param materialParamsMatrix the material law to calculate the sw as inverse of capillary pressure function
     * \param materialParamsFracture the material law to calculate the sw as inverse of capillary pressure function
     *
     * This method is called by updateFracture
     */
    void interfaceCondition(const PrimaryVariables &priVars, const MaterialLawParams &materialParamsMatrix, const MaterialLawParams &materialParamsFracture)
    {
        // calculate the capillary pressure in the fracture and the entry pressure in the matrix
        // the wetting saturation in the fracture is 1 - Sn
        Scalar pcFracture = MaterialLaw::pc(materialParamsFracture, 1 - priVars[saturationIdx]);
        Scalar pEntryMatrix = MaterialLaw::pc(materialParamsMatrix, 1.0);

        /*2nd condition Niessner, J., R. Helmig, H. Jakobs, and J.E. Roberts. 2005, eq.10
        * if the capillary pressure in the fracture is smaller than the entry pressure
        * in the matrix than in the matrix
        * */
        if (pcFracture <= pEntryMatrix)
        {
            satWMatrix_ = 1.0;
            satNMatrix_ = 1 - satWMatrix_;
        }
        //3rd condition Niessner, J., R. Helmig, H. Jakobs, and J.E. Roberts. 2005, eq.10
        else
        {
            /*
             * Inverse capillary pressure function SwM = pcM^(-1)(pcF(SwF))
             */
            satWMatrix_ = MaterialLaw::sw(materialParamsMatrix, pcFracture);
            satNMatrix_ = 1 - satWMatrix_;
        }
    }

    /*!
     * \brief Calculates the volume of the fracture inside the SCV
     */
    Scalar calculateSCVFractureVolume ( const Problem &problem,
                                        const Element &element,
                                        const FVElementGeometry &fvGeometry,
                                        int scvIdx)
    {
        Scalar volSCVFracture;
        const auto geometry = element.geometry();
        Dune::GeometryType geomType = geometry.type();
        const ReferenceElement &refElement = ReferenceElements::general(geomType);

        for (int fIdx=0; fIdx<refElement.size(1); fIdx++)
        {
            SCVFace face = fvGeometry.subContVolFace[fIdx];
            int i=face.i;
            int j=face.j;

            if (problem.spatialParams().isEdgeFracture(element, fIdx)
                && (i == scvIdx || j == scvIdx))
            {
                Scalar fracture_width = problem.spatialParams().fractureWidth();

                const GlobalPosition global_i = geometry.corner(i);
                const GlobalPosition global_j = geometry.corner(j);
                GlobalPosition diff_ij = global_j;
                diff_ij -=global_i;
                //fracture length in the subcontrol volume is half d_ij
                Scalar fracture_length = 0.5*diff_ij.two_norm();

                volSCVFracture += 0.5 * fracture_length * fracture_width;
            }
        }
        return volSCVFracture;
    }

    /*!
     * \brief Returns the effective saturation fracture of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
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
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobilityMatrix_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobilityFracture(int phaseIdx) const
    { return mobilityFracture_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the fracture.
     */
    Scalar porosityFracture() const
    { return porosityFracture_; }

    /*!
     * \brief Returns the average permeability within the fracture.
     */
    Scalar permeabilityFracture() const
    { return permeabilityFracture_; }

protected:
    FluidState fluidState_;
    FluidState fluidStateFracture_;

    Scalar porosityFracture_;
    Scalar permeabilityFracture_;

    Scalar mobilityMatrix_[numPhases];
    Scalar mobilityFracture_[numPhases];

    Scalar satWMatrix_;
    Scalar satNMatrix_;

    Scalar satWFracture_;
    Scalar satNFracture_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // end namespace

#endif // DUMUX_MODELS_2PDFM_VOLUME_VARIABLES_HH
