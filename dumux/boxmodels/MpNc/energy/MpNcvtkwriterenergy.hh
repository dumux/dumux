// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
/*!
 * \file
 *
 * \brief VTK writer module for the energy related quantities of the
 *        MpNc model.
 */
#ifndef DUMUX_MPNC_VTK_WRITER_ENERGY_HH
#define DUMUX_MPNC_VTK_WRITER_ENERGY_HH

#include "../MpNcvtkwritermodule.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 *
 * \brief VTK writer module for the energy related quantities of the
 *        MpNc model.
 *
 * This is the specialization for the case without energy.
 */
template<class TypeTag,
         bool enableEnergy /* = false */,
         bool enableKineticEnergy /* = false */>
class MPNCVtkWriterEnergy : public MPNCVtkWriterModule<TypeTag>
{
    static_assert(enableKineticEnergy == false,
                  "If you enable kinetic energy transfer between fluids, you"
                  "also have to enable the energy in general!");

    typedef MPNCVtkWriterModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;
    typedef typename ParentType::ComponentBuffer ComponentBuffer;
    typedef typename ParentType::PhaseComponentBuffer PhaseComponentBuffer;

public:
    MPNCVtkWriterEnergy(const Problem &problem)
        : ParentType(problem)
    {
        temperatureOutput_ = GET_PARAM(TypeTag, bool, MPNC, VtkAddTemperatures);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        if (temperatureOutput_) this->resizeScalarBuffer_(temperature_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const Element &elem,
                        const FVElementGeometry &fvElemGeom,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        int n = elem.template count<dim>();
        for (int i = 0; i < n; ++i) {
            int I = this->problem_.vertexMapper().map(elem, i, dim);
            const VolumeVariables &volVars = elemVolVars[i];

            if (temperatureOutput_)
                temperature_[I] = volVars.fluidState().temperature(0);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if (temperatureOutput_)
            this->commitScalarBuffer_(writer, "T", temperature_);
    }

private:
    bool temperatureOutput_;

    ScalarBuffer temperature_;
};

/*!
 * \ingroup MPNCModel
 *
 * \brief VTK writer module for the energy related quantities of the
 *        MpNc model.
 *
 * This is the specialization for the case with an energy equation but
 * local thermal equilibrium. (i.e. no kinetic energy transfer)
 */
template<class TypeTag>
class MPNCVtkWriterEnergy<TypeTag, /* enableEnergy = */ true, /* enableKineticEnergy = */ false >
    : public MPNCVtkWriterModule<TypeTag>
{
    typedef MPNCVtkWriterModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;
    typedef typename ParentType::ComponentBuffer ComponentBuffer;
    typedef typename ParentType::PhaseComponentBuffer PhaseComponentBuffer;


public:
    MPNCVtkWriterEnergy(const Problem &problem)
        : ParentType(problem)
    {
        temperatureOutput_ = GET_PARAM(TypeTag, bool, MPNC, VtkAddTemperatures);
        enthalpyOutput_ = GET_PARAM(TypeTag, bool, MPNC, VtkAddEnthalpies);
        internalEnergyOutput_ = GET_PARAM(TypeTag, bool, MPNC, VtkAddInternalEnergies);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        if (temperatureOutput_) this->resizeScalarBuffer_(temperature_);
        if (enthalpyOutput_) this->resizePhaseBuffer_(enthalpy_);
        if (internalEnergyOutput_) this->resizePhaseBuffer_(internalEnergy_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const Element &elem,
                        const FVElementGeometry &fvElemGeom,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        int n = elem.template count<dim>();
        for (int i = 0; i < n; ++i) {
            int I = this->problem_.vertexMapper().map(elem, i, dim);
            const VolumeVariables &volVars = elemVolVars[i];

            if (temperatureOutput_) temperature_[I] = volVars.fluidState().temperature(0);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (enthalpyOutput_)
                    enthalpy_[phaseIdx][I] = volVars.fluidState().temperature(phaseIdx);
                if (internalEnergyOutput_)
                    internalEnergy_[phaseIdx][I] = volVars.fluidState().internalEnergy(phaseIdx);
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if (temperatureOutput_)
            this->commitScalarBuffer_(writer, "T", temperature_);
        if (enthalpyOutput_)
            this->commitPhaseBuffer_(writer, "h_%s", enthalpy_);
        if (internalEnergyOutput_)
            this->commitPhaseBuffer_(writer, "u_%s", internalEnergy_);
    }

private:
    bool temperatureOutput_;
    bool enthalpyOutput_;
    bool internalEnergyOutput_;

    ScalarBuffer temperature_;
    PhaseBuffer enthalpy_;
    PhaseBuffer internalEnergy_;
};

}

#endif
