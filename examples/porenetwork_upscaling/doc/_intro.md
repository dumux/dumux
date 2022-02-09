# Determining the upscaled properties of a pore network

__In this example, you will learn how to__

* simulate creeping/non-creeping flow on a pore network by applying a pressure gradient in a given direction
* perform an upscaling in order to determine the flow properties of the porous medium such as:
the Darcy (interinsic) single-phase permeability $`\mathbf{K}`$ [m$`^2`$] using the creeping flow simulation, the Forchheimer permeability $`\mathbf{K}`$ [m$`^2`$] and the Forchheimer coefficient $`\mathbf{\beta}`$ [m$`^{-1}`$] using the non-creeping flow simulation.


__Result__.
As a result of the creeping flow simulation of this example, you will get the Darcy (interinsic) single-phase permeabilities for each spatial direction $`K_{xx}`$, $`K_{yy}`$, $`K_{zz}`$ [m$`^2`$] as direct output on your terminal as following:

```
X-direction:
-- Darcy permeability = 3.326e-12 m^2

```

The pressure distribution throughout the pore network as well as pore-network characteristics will also be written to a vtp output file that can be viewed with ParaView.
Figure 1 shows the pressure distribution within the pore network for the case of flow in x-direction.

<figure>
    <center>
        <img src="img/result.png" alt="Pore-network pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Pressure distribution within the pore network for flow in x-direction. </figcaption>
    </center>
</figure>

The non-creeping flow simulation additionally gives you Forchheimer permeability and coefficient for each spatial direction as in the example below for the x-direction:

```
X-direction:
-- Darcy permeability = 3.326e-12 m^2
-- Forchheimer permeability = 3.270e-12 m^2
-- Forchheimer coefficient = 8.666e+04 m^-1
```

Furthermore, the ratio of apparent permeability to Darcy permeability is plotted versus the Forchheimer number for three spatial dimension as figure 2 shows.

<figure>
    <center>
        <img src="img/permeability_ratio_versus_forchheimer_number.png" alt="Permeability ratio vs. Forchheimer number" width="60%"/>
        <figcaption> <b> Fig.2 </b> - Variation of apparent to darcy permeability ratio versus Forchheimer number. </figcaption>
    </center>
</figure>

After building the executable, use the keyword `Problem.AssumeCreepingFlow` in params.input file to select the flow type to be simulated (i.e. set it true for creeping flow and false for non-creeping flow simulations). Then, run the simulation with `./example_pnm1p_upscaling`.

__Table of contents__. This description is structured as follows:

[[_TOC_]]


## Problem setup

We consider a single-phase problem within a randomly generated pore network of 20x20x20 pores cubed from which some of the pore throats were [deleted randomly](https://doi.org/10.1029/2010WR010180).
The inscribed pore body radii follow a truncated log-normal distribution.
To calculate the upscaled properties, $`15`$ pressure differences in the range of $`1`$ to $`10^{10}`$ Pa are applied sequentially in every direction, while all lateral sides are closed.
The resulting mass flow rates are then used  to determine upscaled properties as described [later](upscalinghelper.md).

## Mathematical and numerical model

In this example we are using the single-phase pore-network model of DuMu<sup>x</sup>. We require mass conservation at each pore body $`i`$:

```math
 \sum_j Q_{ij} = 0,
```
where $`Q_{ij}`$ is the discrete volume flow rate in a throat connecting pore bodies $`i`$ and $`j`$. In case of creeping flow, the pore network model considers a Hagen-Poiseuille-type law to relate the volume flow from on pore body to another to discrete pressure drops $`\Delta p = p_i - p_j`$ between the pore bodies usig throat conductance, $`g_{ij}`$.

```math
 Q_{ij} = g_{ij} (p_i - p_j),
```
or 

```math
 (p_i - p_j) = Q_{ij} / g_{ij} .
```

In the simulation of non-creeping flow, to capture inertial effects in fluid flow, an extension of the Hagen-Poiseuille-type law which includes expansion and contraction of flow moving from a throat to a pore and vice versa is used. 

```math
  (p_i - p_j) = Q_{ij} / g_{ij}  + (C_{exp} + C_{cont})Q^2,
```

where $`C_{exp}`$ and $`C_{cont}`$ are expansion and contraction coefficients. 

# Implementation & Post processing

In the following, we take a closer look at the source files for this example.
We will discuss the different parts of the code in detail subsequently.

```
└── pnmpermeabilityupscaling/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameters
    ├── properties.hh           -> compile time settings for the simulation
    ├── problem.hh              -> boundary & initial conditions for the simulation
    ├── spatialparams.hh        -> parameter distributions for the simulation
    └── upscalinghelper.hh      -> helper for the upscaling calculations and writing the results
```
