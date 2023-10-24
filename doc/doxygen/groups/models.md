# Doxygen subgroups of the model group

<!-- PorousmediumflowModels -->

@defgroup PorousmediumflowModels Porous medium flow (Darcy)
@brief Single and multi-phase models for flow and transport in porous materials
@ingroup Models

<!-- PorousmediumflowModels subgroups -->

@defgroup OnePModel 1p
@brief Single-phase (immiscible) Darcy flow
@copydoc dumux/porousmediumflow/1p/model.hh
@ingroup PorousmediumflowModels

@defgroup OnePNCModel 1pnc
@brief Single-phase, multi-component Darcy flow
@copydoc dumux/porousmediumflow/1pnc/model.hh
@ingroup PorousmediumflowModels

@defgroup OnePNCMinModel 1pncmin
@brief Single-phase, multi-component Darcy flow with mineralization
@copydoc dumux/porousmediumflow/1pncmin/model.hh
@ingroup PorousmediumflowModels

@defgroup TwoPModel 2p
@brief Two-phase (immiscible) Darcy flow
@copydoc dumux/porousmediumflow/2p/model.hh
@ingroup PorousmediumflowModels

@defgroup TwoPOneCModel 2p1c
@brief Two-phase, one-component Darcy flow
@copydoc dumux/porousmediumflow/2p1c/model.hh
@ingroup PorousmediumflowModels

@defgroup TwoPTwoCModel 2p2c
@brief Two-phase, two-component Darcy flow
@copydoc dumux/porousmediumflow/2p2c/model.hh
@ingroup PorousmediumflowModels

@defgroup TwoPNCModel 2pnc
@brief two-phase, multi-component Darcy flow
@copydoc dumux/porousmediumflow/2pnc/model.hh
@ingroup PorousmediumflowModels

@defgroup TwoPNCMinModel 2pncmin
@brief Two-phase, multi-component Darcy flow with mineralization
@copydoc dumux/porousmediumflow/2pncmin/model.hh
@ingroup PorousmediumflowModels

@defgroup ThreePModel 3p
@brief Three-phase (immiscible) Darcy flow
@copydoc dumux/porousmediumflow/3p/model.hh
@ingroup PorousmediumflowModels

@defgroup ThreePThreeCModel 3p3c
@brief Three-phase, three-component Darcy flow
@copydoc dumux/porousmediumflow/3p3c/model.hh
@ingroup PorousmediumflowModels

@defgroup ThreePWaterOilModel 3pwateroil
@brief Three-phase, two-component Darcy flow with water (liquid & gas) and oil
@copydoc dumux/porousmediumflow/3pwateroil/model.hh
@ingroup PorousmediumflowModels

@defgroup BoxDFMModel boxdfm
@brief Vertex-centered, continuous-pressure, conforming lower-dimensional discrete-fracture model
@copydoc dumux/porousmediumflow/boxdfm/model.hh
@ingroup PorousmediumflowModels

@defgroup CO2Model CO2
@brief Two-phase, two-component Darcy flow specialized for supercritical CO<sub>2</sub> storage
@copydoc dumux/porousmediumflow/co2/model.hh
@ingroup PorousmediumflowModels

@defgroup MineralizationModel mineralization
@brief Model adding components that can precipitate as a solid phase to a standard Darcy flow model
@copydoc dumux/porousmediumflow/mineralization/model.hh
@ingroup PorousmediumflowModels

@defgroup MPNCModel mpnc
@brief Generalized multi-phase, multi-component Darcy flow
@copydoc dumux/porousmediumflow/mpnc/model.hh
@ingroup PorousmediumflowModels

@defgroup NonEquilibriumModel nonequilibrium
@brief Model that adds nonequilibrium equations to another porous medium flow model (only used in MPNCModel currently)
@copydoc dumux/porousmediumflow/nonequilibrium/model.hh
@ingroup PorousmediumflowModels

@defgroup ThermalNonEquilibriumModel thermal-nonequilibrium
@brief Model that adapts the energy localresidual to thermal nonequilibrium
@copydoc dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh
@ingroup NonEquilibriumModel

@defgroup NIModel non-isothermal
@brief Model that adds an energy equation (thermal equilibrium) to another porous medium flow model
@copydoc dumux/porousmediumflow/nonisothermal/model.hh
@ingroup PorousmediumflowModels

@defgroup RichardsModel Richards
@brief Richards flow
@copydoc dumux/porousmediumflow/richards/model.hh
@ingroup PorousmediumflowModels

@defgroup ExtendedRichardsModel extended Richards' equation
@brief extended Richards' equation
@copydoc dumux/porousmediumflow/richardsextended/model.hh
@ingroup PorousmediumflowModels

@defgroup RichardsNCModel Richards multi-component
@brief Richards multi-component flow
@copydoc dumux/porousmediumflow/richardsnc/model.hh
@ingroup PorousmediumflowModels

@defgroup SolidEnergyModel solid-energy
@brief Energy equation for the solid (general heat equation)
@copydoc dumux/porousmediumflow/solidenergy/model.hh
@ingroup PorousmediumflowModels

@defgroup TracerModel tracer
@brief Multi-component advection-diffusion-reaction model with given velocity field
@copydoc dumux/porousmediumflow/tracer/model.hh
@ingroup PorousmediumflowModels

<!-- FreeflowModels -->

@defgroup FreeflowModels Free flow (Navier-Stokes)
@brief Single-phase models based on the Navier-Stokes equation
@ingroup Models

<!-- FreeflowModels subgroups -->

@defgroup NavierStokesModel Navier-Stokes
@brief Single-phase Navier-Stokes flow
@copydoc dumux/freeflow/navierstokes/model.hh
@ingroup FreeflowModels

@defgroup RANSModel Reynolds-Averaged Navier-Stokes
@brief Single-phase Reynolds-Averaged Navier-Stokes flow
@copydoc dumux/freeflow/rans/model.hh
@ingroup FreeflowModels

@defgroup ZeroEqModel 0-Eq. Models
@brief Zero-equation or algebraic turbulence models
@copydoc dumux/freeflow/rans/zeroeq/model.hh
@ingroup RANSModel

@defgroup OneEqModel 1-Eq. Models
@brief One-equation turbulence model by Spalart-Allmaras
@copydoc dumux/freeflow/rans/oneeq/model.hh
@ingroup RANSModel

@defgroup TwoEqModel 2-Eq. Models
@brief Two-equation turbulence models
@ingroup RANSModel

@defgroup KEpsilonModel K-epsilon model
@brief K-epsilon model
@copydoc dumux/freeflow/rans/twoeq/kepsilon/model.hh
@ingroup TwoEqModel

@defgroup KOmegaModel K-omega model
@brief K-omega model
@copydoc dumux/freeflow/rans/twoeq/komega/model.hh
@ingroup TwoEqModel

@defgroup LowReKEpsilonModel Low-Re k-epsilon model
@brief Low-Re k-epsilon model
@copydoc dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh
@ingroup TwoEqModel

@defgroup SSTModel SST model
@brief SST model
@copydoc dumux/freeflow/rans/twoeq/sst/model.hh
@ingroup TwoEqModel

@defgroup FreeflowNCModel Compositional
@brief Single-phase multi-component free-flow flow models
@copydoc dumux/freeflow/compositional/navierstokesncmodel.hh
@ingroup FreeflowModels

@defgroup FreeflowNIModel Nonisothermal
@brief An energy equation adaptor for isothermal free-flow models
@copydoc dumux/freeflow/nonisothermal/model.hh
@ingroup FreeflowModels

<!-- ShallowWaterModels -->

@defgroup ShallowWaterModels Shallow water flow
@brief Two-dimensional shallow water flow (depth-averaged)
@copydoc dumux/freeflow/shallowwater/model.hh
@ingroup Models

<!-- GeomechanicsModels -->

@defgroup GeomechanicsModels Geomechanics
@brief Models taking into account solid deformation
@ingroup Models

<!-- GeomechanicsModels subgroups -->

@defgroup Elastic Solid mechanics w/o fluid pressure
@brief Models linear elastic deformation of a solid. Disregards fluid pressure.
@copydoc dumux/geomechanics/elastic/model.hh
@ingroup GeomechanicsModels

@defgroup PoroElastic Solid mechanics with fluid pressure
@brief Models linear elastic deformation of a solid. Takes fluid pressure into account.
@copydoc dumux/geomechanics/poroelastic/model.hh
@ingroup GeomechanicsModels

@defgroup Hyperelastic Hyperelastic solid mechanics w/o fluid pressure
@brief Models nonlinear deformation of an elastic solid.
@copydoc dumux/geomechanics/hyperelastic/model.hh
@ingroup GeomechanicsModels

<!-- PoreNetworkModels -->

@defgroup PoreNetworkModels Pore network
@brief Single and multi-phase models for flow and transport in pore networks
@ingroup Models

<!-- PoreNetworkModels subgroups -->

@defgroup PNMOnePModel 1p
@brief Single-phase (immiscible) flow
@copydoc dumux/porenetwork/1p/model.hh
@ingroup PoreNetworkModels

@defgroup PNMOnePNCModel 1pnc
@brief Single-phase, multi-component flow
@copydoc dumux/porenetwork/1pnc/model.hh
@ingroup PoreNetworkModels

@defgroup PNMTwoPModel 2p
@brief Two-phase (immiscible) flow
@copydoc dumux/porenetwork/2p/model.hh
@ingroup PoreNetworkModels

@defgroup PNMSolidEnergyModel solidenergy
@brief Energy equation for the solid (heat equation)
@copydoc dumux/porenetwork/solidenergy/model.hh
@ingroup PoreNetworkModels
