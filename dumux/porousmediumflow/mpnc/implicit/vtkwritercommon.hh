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
/*!
 * \file
 *
 * \brief VTK writer module for the common quantities of the MpNc
 *        model.
 */
#ifndef DUMUX_MPNC_VTK_WRITER_COMMON_HH
#define DUMUX_MPNC_VTK_WRITER_COMMON_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "vtkwritermodule.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 *
 * \brief VTK writer module for the common quantities of the MpNc
 *        model.
 */
template<class TypeTag>
class MPNCVtkWriterCommon : public MPNCVtkWriterModule<TypeTag>
{
    typedef MPNCVtkWriterModule<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename ParentType::ScalarVector ScalarVector;
    typedef typename ParentType::PhaseVector PhaseVector;
    typedef typename ParentType::ComponentVector ComponentVector;
    typedef typename ParentType::PhaseComponentMatrix PhaseComponentMatrix;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimWorldVector;
    typedef Dune::BlockVector<DimWorldVector> DimWorldField;
    typedef std::array<DimWorldField, numPhases> PhaseDimWorldField;
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    MPNCVtkWriterCommon(const Problem &problem)
    : ParentType(problem), velocityOutput_(problem)
    {
        porosityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPorosity);
        permeabilityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPermeability);
        boundaryTypesOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddBoundaryTypes);
        saturationOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddSaturations);
        pressureOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddPressures);
        densityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddDensities);
        mobilityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddMobilities);
        averageMolarMassOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddAverageMolarMass);
        massFracOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddMassFractions);
        moleFracOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddMoleFractions);
        molarityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddMolarities);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        if (porosityOutput_)
            this->resizeScalarBuffer_(porosity_, isBox);
        if (permeabilityOutput_)
            this->resizeScalarBuffer_(permeability_);
        if (boundaryTypesOutput_)
            this->resizeScalarBuffer_(boundaryTypes_, isBox);

        if (velocityOutput_.enableOutput()) {
            Scalar numDofs = this->problem_.model().numDofs();
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                velocity_[phaseIdx].resize(numDofs);
                velocity_[phaseIdx] = 0;
            }
        }

        if (saturationOutput_) this->resizePhaseBuffer_(saturation_, isBox);
        if (pressureOutput_) this->resizePhaseBuffer_(pressure_, isBox);
        if (densityOutput_) this->resizePhaseBuffer_(density_, isBox);
        if (mobilityOutput_) this->resizePhaseBuffer_(mobility_, isBox);
        if (averageMolarMassOutput_) this->resizePhaseBuffer_(averageMolarMass_, isBox);
        if (moleFracOutput_) this->resizePhaseComponentBuffer_(moleFrac_, isBox);
        if (massFracOutput_) this->resizePhaseComponentBuffer_(massFrac_, isBox);
        if (molarityOutput_) this->resizePhaseComponentBuffer_(molarity_, isBox);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     *
     *        \param element The finite element
     *        \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *        \param elemVolVars The volume variables of the current element
     *        \param elemBcTypes The types of the boundary conditions for all vertices of the element
     */
    void processElement(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
        {
            int dofIdxGlobal = this->problem_.model().dofMapper().subIndex(element, scvIdx, dofCodim);

            const VolumeVariables &volVars = elemVolVars[scvIdx];

            if (porosityOutput_) porosity_[dofIdxGlobal] = volVars.porosity();

            // works for scalar permeability in spatialparameters
            if (permeabilityOutput_) permeability_[dofIdxGlobal] = this->problem_.spatialParams().intrinsicPermeability(element,fvGeometry,scvIdx);

            // calculate a single value for the boundary type: use one
            // bit for each equation and set it to 1 if the equation
            // is used for a dirichlet condition
            int tmp = 0;
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (elemBcTypes[scvIdx].isDirichlet(eqIdx))
                    tmp += (1 << eqIdx);
            }
            if (boundaryTypesOutput_) boundaryTypes_[dofIdxGlobal] = tmp;

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (saturationOutput_) saturation_[phaseIdx][dofIdxGlobal] = volVars.fluidState().saturation(phaseIdx);
                if (pressureOutput_) pressure_[phaseIdx][dofIdxGlobal] = volVars.fluidState().pressure(phaseIdx);
                if (densityOutput_) density_[phaseIdx][dofIdxGlobal] = volVars.fluidState().density(phaseIdx);
                if (mobilityOutput_) mobility_[phaseIdx][dofIdxGlobal] = volVars.mobility(phaseIdx);
                if (averageMolarMassOutput_) averageMolarMass_[phaseIdx][dofIdxGlobal] = volVars.fluidState().averageMolarMass(phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (moleFracOutput_) moleFrac_[phaseIdx][compIdx][dofIdxGlobal] = volVars.fluidState().moleFraction(phaseIdx, compIdx);
                    if (massFracOutput_) massFrac_[phaseIdx][compIdx][dofIdxGlobal] = volVars.fluidState().massFraction(phaseIdx, compIdx);
                    if (molarityOutput_) molarity_[phaseIdx][compIdx][dofIdxGlobal] = volVars.fluidState().molarity(phaseIdx, compIdx);
                }
            }
        }

        // velocity output
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            velocityOutput_.calculateVelocity(velocity_[phaseIdx], elemVolVars, fvGeometry, element, phaseIdx);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if (saturationOutput_)
            this->commitPhaseBuffer_(writer, "S_%s", saturation_, isBox);

        if (pressureOutput_)
            this->commitPhaseBuffer_(writer, "p_%s", pressure_, isBox);

        if (densityOutput_)
            this->commitPhaseBuffer_(writer, "rho_%s", density_, isBox);

        if (averageMolarMassOutput_)
            this->commitPhaseBuffer_(writer, "M_%s", averageMolarMass_, isBox);

        if (mobilityOutput_)
            this->commitPhaseBuffer_(writer, "lambda_%s", mobility_, isBox);

        if (porosityOutput_)
            this->commitScalarBuffer_(writer, "porosity", porosity_, isBox);

        if (boundaryTypesOutput_)
            this->commitScalarBuffer_(writer, "boundary types", boundaryTypes_, isBox);

        if (moleFracOutput_)
            this->commitPhaseComponentBuffer_(writer, "x_%s^%s", moleFrac_, isBox);

        if (massFracOutput_)
            this->commitPhaseComponentBuffer_(writer, "X_%s^%s", massFrac_, isBox);

        if(molarityOutput_)
            this->commitPhaseComponentBuffer_(writer, "c_%s^%s", molarity_, isBox);

        if (velocityOutput_.enableOutput()) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                std::ostringstream oss;
                oss << "velocity_" << FluidSystem::phaseName(phaseIdx);
                writer.attachDofData(velocity_[phaseIdx],
                                     oss.str(), isBox, dim);
            }
        }
    }

private:
    bool porosityOutput_;
    bool permeabilityOutput_ ;
    bool boundaryTypesOutput_;
    bool saturationOutput_;
    bool pressureOutput_;
    bool densityOutput_;
    bool mobilityOutput_;
    bool massFracOutput_;
    bool moleFracOutput_;
    bool molarityOutput_;
    bool averageMolarMassOutput_;

    PhaseVector saturation_;
    PhaseVector pressure_;
    PhaseVector density_;
    PhaseVector mobility_;
    PhaseVector averageMolarMass_;

    PhaseDimWorldField velocity_;

    ScalarVector porosity_;
    ScalarVector permeability_;
    ScalarVector temperature_;
    ScalarVector boundaryTypes_;

    PhaseComponentMatrix moleFrac_;
    PhaseComponentMatrix massFrac_;
    PhaseComponentMatrix molarity_;

    ComponentVector fugacity_;

    ImplicitVelocityOutput<TypeTag> velocityOutput_;
};

}

#endif
