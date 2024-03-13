// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMConstraint
 */

#ifndef DUMUX_CONSTRAINT_MODEL_HH
#define DUMUX_CONSTRAINT_MODEL_HH

// The property tag is simply an empty struct with the name `PNMConstraintModel`
namespace Dumux::Properties::TTag {
//! The pnm constraint model tag
struct PNMConstraintModel {};
} // end namespace Dumux::Properties::TTag

namespace Dumux {

/*!
 * \ingroup PNMConstraintModel
 * \brief Adds I/O fields specific to the constraint model
 */
class PNMConstraintIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& volVars){ return volVars.theta(); },
                              "Theta");
    }

    template <class ModelTraits = void, class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        return "Theta";
    }
};

/*!
 * \ingroup PNMConstraintModel
 * \brief Volume averaged quantities required by the model
 */
template <class Traits>
class PNMConstraintVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
    }

    /*!
     * \brief Returns the image intensity of the control volume.
     */
    Scalar theta() const
    { return priVars_[Indices::thetaIdx]; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

/*!
 * \ingroup PNMConstraintModel
 * \brief Element-wise calculation of the residual and its derivatives
 */
template<class TypeTag>
class PNMConstraintModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // For a simple algebraic eq. we don't have a storage
        NumEqVector storage(0.0);
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // For a simple algebraic eq. we don't have a storage
        // ToDo also call flux constraint
        NumEqVector flux(0.0);
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        // add constraint equation
        source += problem.volumeConstraint(element, fvGeometry, elemVolVars, scv);

        return source;
    }

};
} // end namespace Dumux

// The model properties
namespace Dumux::Properties {

// The type of the local residual is the class defined above
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMConstraintModel>
{ using type = PNMConstraintModelLocalResidual<TypeTag>; };

// The default scalar type is double
// we compute with double precision floating point numbers
template<class TypeTag>
struct Scalar<TypeTag, TTag::PNMConstraintModel>
{ using type = double; };

// The model traits specify some information about our equation system.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMConstraintModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int thetaIdx = 0;
            static constexpr int constraintEqIdx = 0;
        };

        static constexpr int numEq() { return 1; }
    };
};

// The primary variable vector has entries of type `Scalar` and is
// as large as the number of equations
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PNMConstraintModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMConstraintModel> { using type = PNMConstraintIOFields; };

// The volume variables that are used in our model
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMConstraintModel>
{
    struct Traits
    {
        using ModelTraits
            = GetPropType<TypeTag, Properties::ModelTraits>;
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };

    using type = PNMConstraintVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
