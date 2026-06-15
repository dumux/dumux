// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 * \copydoc Dumux::BedloadVolumeVariables
 */
#ifndef DUMUX_BEDLOAD_VOLUME_VARIABLES_HH
#define DUMUX_BEDLOAD_VOLUME_VARIABLES_HH

#include <dumux/material/sediment/bedloadtransport/secondarycurrents.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTransportModel
 * \brief Volume variables for the bedload transport model.
 *
 * \tparam Traits Class encapsulating types to be used by the vol vars
 */
template<class Traits>
class BedloadVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the type encapsulating primary variable indices
    using Indices = typename ModelTraits::Indices;
    static constexpr int numGrainClasses = ModelTraits::numGrainClasses();

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        priVars_ = elemSol[scv.localDofIndex()];
        static const int nGrainClasses = getParam<int>("Sediment.NumberGrainClasses");
        // It is important to update massFractions_ first
        massFractions_.assign(numGrainClasses, 0.0);
        Scalar temp = 0.0;
        for (int i=0; i<nGrainClasses; i++)
        {
            massFractions_[i] = sedimentMass(i)/sedimentMass();
            temp += massFractions_[i] / problem.spatialParams().grainDensity(i);
        }
        harmonicAverageGrainDensity_ = 1 / temp;
        waterDepth_ = problem.spatialParams().waterDepth(element, scv);
        velocity_[0] = problem.spatialParams().velocityX(element, scv);
        velocity_[1] = problem.spatialParams().velocityY(element, scv);
        bottomActiveLayer_ = problem.bottomActiveLayer(element, scv);
        layerThickness_ = sedimentMass() / (harmonicAverageGrainDensity_ * (1 - problem.spatialParams().porosity()) * scv.volume());
        tanAngleSecondaryCurrents_ = calculateTanAngleSecondaryCurrents(waterDepth(), problem.spatialParams().curvatureRadius(scv));
        // It is important to update bottomShearStress_ before bedloadDischarge_, since the bedloadDischarge_ uses bottomShearstress_
        bottomShearStress_ = Dune::FieldVector<Scalar, 2> {problem.spatialParams().bottomShearStressX(element, scv),
                                                           problem.spatialParams().bottomShearStressY(element, scv)};
        frictionParameter_ = problem.spatialParams().frictionParameter(scv);
        bedloadDischarge_.assign(numGrainClasses, Dune::FieldVector<Scalar, 2>{0.0, 0.0});
        for (int i=0; i<nGrainClasses; i++)
        {
            bedloadDischarge_[i] = problem.spatialParams().bedloadFormula(i).bedloadDischarge(*this);
        }
    }

    /*!
     * \brief Return the bedload discharge inside the sub-control volume for the specified grain class.
     *
     * \param grainClassIdx index of the grain class starting with 0
     */
    Dune::FieldVector<Scalar, 2> bedloadDischarge(int grainClassIdx) const
    {
        return bedloadDischarge_[grainClassIdx];
    }

    /*!
     * \brief Return the bed surface inside the sub-control volume.
     *
     */
    Scalar bedSurface() const
    {
        return this->bottomActiveLayer() + this->layerThickness();
    }

    /*!
     * \brief Return the bottom of the active layer.
     *
     */
    Scalar bottomActiveLayer() const
    {
        return bottomActiveLayer_;
    }

    /*!
     * \brief Return the bottom shear stress inside the sub-control volume.
     *
     */
    Dune::FieldVector<Scalar, 2> bottomShearStress() const
    {
        return bottomShearStress_;
    }

    /*!
     * \brief Return the extrusion factor (dummy variable).
     */
    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \brief Return the friction coefficient.
     *
     */
    Scalar frictionParameter() const
    {
        return frictionParameter_;
    }

    /*!
     * \brief Return the average grain density of the active layer inside the sub-control volume.
     *
     * The average grain density may change over time if the composition of the active layer changes.
     * Note that the harmonic average is calculated.
     *
     */
    Scalar harmonicAverageGrainDensity() const
    {
        return harmonicAverageGrainDensity_;
    }

    /*!
     * \brief Return the thickness of the active layer.
     */
    Scalar layerThickness() const
    {
        return layerThickness_;
    }

    /*!
     * \brief Return the mass fraction of the specified grain class in the active layer.
     *
     * \param grainClassIdx index of the grain class starting with 0
     */
    Scalar massFraction(int grainClassIdx) const
    {
        return massFractions_[grainClassIdx];
    }

    /*!
     * \brief Return the sediment mass in the active layer.
     */
    Scalar sedimentMass() const
    {
        Scalar sedimentMass = 0.0;
        for (int i = 0; i < numGrainClasses; i++)
        {
            sedimentMass += this->sedimentMass(i);
        }
        return sedimentMass;
    }

    /*!
     * \brief Return the sediment mass in the active layer for the specified grain class.
     *
     * \param grainClassIdx index of the grain class starting with 0
     */
    Scalar sedimentMass(int grainClassIdx) const
    {
        return priVars_[grainClassIdx];
    }

    /*!
     * \brief Return the tangent of the secondary currents angle.
     *
     */
    Scalar tanAngleSecondaryCurrents() const
    {
        return tanAngleSecondaryCurrents_;
    }

    /*!
     * \brief Return water velocity component inside the sub-control volume.
     *
     * \param directionIndex index of the direction starting at x = 0
     */
    Scalar velocity(int directionIndex) const
    {
        return velocity_[directionIndex];
    }

    /*!
     * \brief Return water velocity vector inside the sub-control volume.
     */
    Dune::FieldVector<Scalar, 2> velocity() const
    {
        return velocity_;
    }

    /*!
     * \brief Return water detph h inside the sub-control volume.
     *
     */
    Scalar waterDepth() const
    {
        return waterDepth_;
    }

private:
    // data members
    std::vector<Dune::FieldVector<Scalar, 2>> bedloadDischarge_;
    Scalar bottomActiveLayer_;
    Dune::FieldVector<Scalar, 2> bottomShearStress_;
    Scalar frictionParameter_;
    Scalar harmonicAverageGrainDensity_;
    Scalar layerThickness_;
    std::vector<Scalar> massFractions_;
    PrimaryVariables priVars_;
    Scalar tanAngleSecondaryCurrents_;
    Dune::FieldVector<Scalar, 2> velocity_;
    Scalar waterDepth_;
};

} // end namespace Dumux

#endif
