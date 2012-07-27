// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief VTK writer module for the common quantities of the MpNc
 *        model.
 */
#ifndef DUMUX_MPNC_VTK_WRITER_COMMON_HH
#define DUMUX_MPNC_VTK_WRITER_COMMON_HH

#include "mpncvtkwritermodule.hh"

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
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename ParentType::ScalarVector ScalarVector;
    typedef typename ParentType::PhaseVector PhaseVector;
    typedef typename ParentType::ComponentVector ComponentVector;
    typedef typename ParentType::PhaseComponentMatrix PhaseComponentMatrix;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::BlockVector<DimVector> DimField;
    typedef std::tr1::array<DimField, numPhases> PhaseDimField;

public:
    MPNCVtkWriterCommon(const Problem &problem)
        : ParentType(problem)
    {
        porosityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddPorosity);
        boundaryTypesOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddBoundaryTypes);
        saturationOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddSaturations);
        pressureOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddPressures);
        velocityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddVelocities);
        densityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddDensities);
        mobilityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddMobilities);
        averageMolarMassOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddAverageMolarMass);
        massFracOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddMassFractions);
        moleFracOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddMoleFractions);
        molarityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MPNC, VtkAddMolarities);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        if (porosityOutput_)
            this->resizeScalarBuffer_(porosity_);
        if (boundaryTypesOutput_)
            this->resizeScalarBuffer_(boundaryTypes_);

        if (velocityOutput_) {
            Scalar nVerts = this->problem_.gridView().size(dim);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                velocity_[phaseIdx].resize(nVerts);
                velocity_[phaseIdx] = 0;
            }
            this->resizeScalarBuffer_(boxSurface_);
        }

        if (saturationOutput_) this->resizePhaseBuffer_(saturation_);
        if (pressureOutput_) this->resizePhaseBuffer_(pressure_);
        if (densityOutput_) this->resizePhaseBuffer_(density_);
        if (mobilityOutput_) this->resizePhaseBuffer_(mobility_);
        if (averageMolarMassOutput_) this->resizePhaseBuffer_(averageMolarMass_);
        if (moleFracOutput_) this->resizePhaseComponentBuffer_(moleFrac_);
        if (massFracOutput_) this->resizePhaseComponentBuffer_(massFrac_);
        if (molarityOutput_) this->resizePhaseComponentBuffer_(molarity_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const Element &elem,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        int numLocalVertices = elem.geometry().corners();
        for (int localVertexIdx = 0; localVertexIdx < numLocalVertices; ++localVertexIdx) {
            int globalIdx = this->problem_.vertexMapper().map(elem, localVertexIdx, dim);
            const VolumeVariables &volVars = elemVolVars[localVertexIdx];

            if (porosityOutput_) porosity_[globalIdx] = volVars.porosity();

            // calculate a single value for the boundary type: use one
            // bit for each equation and set it to 1 if the equation
            // is used for a dirichlet condition
            int tmp = 0;
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (elemBcTypes[localVertexIdx].isDirichlet(eqIdx))
                    tmp += (1 << eqIdx);
            }
            if (boundaryTypesOutput_) boundaryTypes_[globalIdx] = tmp;

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (saturationOutput_) saturation_[phaseIdx][globalIdx] = volVars.fluidState().saturation(phaseIdx);
                if (pressureOutput_) pressure_[phaseIdx][globalIdx] = volVars.fluidState().pressure(phaseIdx);
                if (densityOutput_) density_[phaseIdx][globalIdx] = volVars.fluidState().density(phaseIdx);
                if (mobilityOutput_) mobility_[phaseIdx][globalIdx] = volVars.mobility(phaseIdx);
                if (averageMolarMassOutput_) averageMolarMass_[phaseIdx][globalIdx] = volVars.fluidState().averageMolarMass(phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (moleFracOutput_) moleFrac_[phaseIdx][compIdx][globalIdx] = volVars.fluidState().moleFraction(phaseIdx, compIdx);
                    if (massFracOutput_) massFrac_[phaseIdx][compIdx][globalIdx] = volVars.fluidState().massFraction(phaseIdx, compIdx);
                    if (molarityOutput_) molarity_[phaseIdx][compIdx][globalIdx] = volVars.fluidState().molarity(phaseIdx, compIdx);
                }
            }
        }

        // calculate velocities if requested by the problem
        if (velocityOutput_) {
            for (int faceIdx = 0; faceIdx < fvGeometry.numEdges; ++ faceIdx) {
                int i = fvGeometry.subContVolFace[faceIdx].i;
                int I = this->problem_.vertexMapper().map(elem, i, dim);

                int j = fvGeometry.subContVolFace[faceIdx].j;
                int J = this->problem_.vertexMapper().map(elem, j, dim);

                FluxVariables fluxVars(this->problem_,
                                       elem,
                                       fvGeometry,
                                       faceIdx,
                                       elemVolVars);

                Scalar scvfArea = fvGeometry.subContVolFace[faceIdx].normal.two_norm();
                scvfArea *= fluxVars.extrusionFactor();

                boxSurface_[I] += scvfArea;
                boxSurface_[J] += scvfArea;

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Dune::FieldVector<Scalar, dim> darcyVelocity;
                    darcyVelocity = fluxVars.velocity(phaseIdx);
                    darcyVelocity *= scvfArea;
                    velocity_[phaseIdx][I] += darcyVelocity;
                    velocity_[phaseIdx][J] += darcyVelocity;

                } // end for all phases
            } // end for all faces
        } // end if velocityOutput_
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if (saturationOutput_)
            this->commitPhaseBuffer_(writer, "S_%s", saturation_);

        if (pressureOutput_)
            this->commitPhaseBuffer_(writer, "p_%s", pressure_);

        if (densityOutput_)
            this->commitPhaseBuffer_(writer, "rho_%s", density_);

        if (averageMolarMassOutput_)
            this->commitPhaseBuffer_(writer, "M_%s", averageMolarMass_);

        if (mobilityOutput_)
            this->commitPhaseBuffer_(writer, "lambda_%s", mobility_);

        if (porosityOutput_)
            this->commitScalarBuffer_(writer, "porosity", porosity_);

        if (boundaryTypesOutput_)
            this->commitScalarBuffer_(writer, "boundary types", boundaryTypes_);

        if (moleFracOutput_)
            this->commitPhaseComponentBuffer_(writer, "x_%s^%s", moleFrac_);

        if (massFracOutput_)
            this->commitPhaseComponentBuffer_(writer, "X_%s^%s", massFrac_);

        if(molarityOutput_)
            this->commitPhaseComponentBuffer_(writer, "c_%s^%s", molarity_);

        if (velocityOutput_) {
            int nVerts = this->problem_.gridView().size(dim);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (int i = 0; i < nVerts; ++i)
                    velocity_[phaseIdx][i] /= boxSurface_[i];
                // commit the phase velocity
                std::ostringstream oss;
                oss << "velocity_" << FluidSystem::phaseName(phaseIdx);
                writer.attachVertexData(velocity_[phaseIdx],
                                        oss.str(),
                                        dim);
            }
        }
    }

private:
    bool porosityOutput_;
    bool boundaryTypesOutput_;
    bool saturationOutput_;
    bool pressureOutput_;
    bool velocityOutput_;
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

    PhaseDimField velocity_;
    ScalarVector boxSurface_;

    ScalarVector porosity_;
    ScalarVector temperature_;
    ScalarVector boundaryTypes_;

    PhaseComponentMatrix moleFrac_;
    PhaseComponentMatrix massFrac_;
    PhaseComponentMatrix molarity_;

    ComponentVector fugacity_;
};

}

#endif
