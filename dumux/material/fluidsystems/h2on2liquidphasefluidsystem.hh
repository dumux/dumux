#ifndef DUMUX_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HH_OLD
#define DUMUX_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HH_OLD

#warning this header is deprecated, use dumux/material/fluidsystems/h2on2.hh instead

#include <dumux/material/fluidsystems/h2on2.hh>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 *
 * \brief A one-phase (water-phase) fluid system with water and nitrogen as components.
 *
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace Dumux::FluidSystems::.
 * An adapter class using Dumux::FluidSystem<TypeTag> is also provided
 * at the end of this file.
 */
template <class Scalar, bool useComplexRelations = true>
DUNE_DEPRECATED_MSG("Class Dumux::FluidSystems::H2ON2LiquidPhase is deprecated. Use Dumux::FluidSystems::H2ON2 instead.")
class H2ON2LiquidPhase
    : public H2ON2<Scalar, useComplexRelations>
{};

} // end namespace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
/*!
 * \brief A one-phase fluid system with water and nitrogen as components.
 *
 * This is an adapter to use Dumux::H2ON2LiquidPhaseFluidSystem<TypeTag>, as is
 * done with most other classes in Dumux.
 */
template<class TypeTag>
DUNE_DEPRECATED_MSG("Class Dumux::H2ON2LiquidPhaseFluidSystem is deprecated. Use Dumux::H2ON2FluidSystem instead.")
class H2ON2LiquidPhaseFluidSystem
: public H2ON2FluidSystem<TypeTag>
{};
#endif

} // end namespace

#endif
