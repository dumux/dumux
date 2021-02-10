/*!
 * \file
 * \ingroup FreeflowTurbulenceModel
 *
 * \brief
 *
 * In addition to the momentum and mass/mole balance equations that form the RANS equations,
 * <B> turbulence models </B> are solved to provide closure approximating the unresolved reynolds stress term
 *
 * During Reynolds averaging where the Navier stokes equations become the RANS eqauations,
 * an additional term arises due to the non averaging invariant convection term
 * (\f$ \nabla \cdot {\boldsymbol{v}} {\boldsymbol{v}}^T\f$)
 *
 * This additional term (\f$ \nabla \cdot \boldsymbol{v}^\prime \boldsymbol{v}^\prime \f$),
 * is a tensor, with 9 entries in 3D systems.
 *
 * Using the linear eddy viscosity assumption presented by Bousinesseq, this tensor is assumed to be
 * isotropic and homogenous, and can be approximated as a viscous stress. This adds an additional term
 * to the viscous term, which is as follows:
 *
 * \f[
 *     \tau = - (\nu_{mol} + \nu_{turb}) \nabla \boldsymbol{v} + \boldsymbol{v}^T
 * \f]
 *
 * In order to define the turbulent eddy viscosity \f$ \nu_{turb} \f$ various turbulence models can be used.
 * They are found within the associated headers.
 */


#ifndef DUMUX_RANS_TURBULENCE_MODEL_HH
#define DUMUX_RANS_TURBULENCE_MODEL_HH

#include "indices.hh"
#include "iofields.hh"

namespace Dumux::Properties {

namespace TTag {
//! The type tag for the base Reynolds-Averaged Navier-Stokes model
struct RANS { using InheritsFrom = std::tuple<ModelProperties>; };
}

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal RANS Models
///////////////////////////////////////////////////////////////////////////
/*!
 * \ingroup RANSModel
 * \brief Traits for the Reynolds-averaged Navier-Stokes model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct RANSModelTraits : NavierStokesModelTraits<dimension>
{
    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }
};

//! The model traits of the isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RANS>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = RANSModelTraits<dim>;
};

//! The model traits of the isothermal model
template<class TypeTag>
struct GridVariables<TypeTag, TTag::RANS>
{
private:
    using GridGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>;
    using GridVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GridFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using NavierStokesGridVariables = FVGridVariables<GridGeometry, GridVolumeVariables, GridFluxVariablesCache>;
public:
    using type = RANSGridVariables<NavierStokesGridVariables>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::RANS> { using type = RANSIOFields; };

}
#endif
