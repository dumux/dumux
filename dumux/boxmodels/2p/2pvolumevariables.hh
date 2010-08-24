// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Quantities required by the twophase box model defined on a vertex.
 */
#ifndef DUMUX_2P_VOLUME_VARIABLES_HH
#define DUMUX_2P_VOLUME_VARIABLES_HH

#include "2pproperties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum {
        pwSn = Indices::pwSn,
        pnSw = Indices::pnSw,

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };


    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef TwoPFluidState<TypeTag> FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        primaryVars_ = priVars;

        asImp().updateTemperature_(priVars,
                                   element,
                                   elemGeom,
                                   scvIdx,
                                   problem);

        // material law parameters
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);

        Scalar p[numPhases];
        Scalar Sn;
        if (int(formulation) == pwSn) {
            Sn = priVars[saturationIdx];
            p[wPhaseIdx] = priVars[pressureIdx];
            p[nPhaseIdx] =
                p[wPhaseIdx] +
                MaterialLaw::pC(materialParams, 1 - Sn);
        }
        else if (int(formulation) == pnSw) {
            Sn = 1 - priVars[saturationIdx];
            p[nPhaseIdx] = priVars[pressureIdx];
            p[wPhaseIdx] =
                p[nPhaseIdx] -
                MaterialLaw::pC(materialParams, 1 - Sn);
        }

        fluidState_.update(Sn, p[wPhaseIdx], p[nPhaseIdx], temperature_);

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, 1 - Sn)
            /
            FluidSystem::phaseViscosity(wPhaseIdx,
                                        temperature_,
                                        p[wPhaseIdx],
                                        fluidState_);
        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, 1 - Sn)
            /
            FluidSystem::phaseViscosity(nPhaseIdx,
                                        temperature_,
                                        p[nPhaseIdx],
                                        fluidState_);

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
    }

    void updateTemperature_(const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem)
    {
        temperature_ = problem.temperature(element, elemGeom, scvIdx);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    /*!
     * \brief Sets the evaluation point used in the by the local jacobian.
     */
    void setEvalPoint(const Implementation *ep)
    { }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.phasePressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

protected:
    PrimaryVariables primaryVars_;
    FluidState fluidState_;
    Scalar porosity_;
    Scalar temperature_;

    Scalar mobility_[numPhases];

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
