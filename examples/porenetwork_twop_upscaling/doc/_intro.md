# Upscaling two-phase flow properties using pore network

__In this example, you will learn how to__

* simulate two-phase flow on a pore network by ???applying a pressure gradient in a given direction???
* perform an upscaling in order to determine the flow properties of the porous medium such as:
capillary pressure -- saturation and relative permeability -- saturation curve.


__Result__.
As a result of the two-phase flow simulation of this example, you will get capillary pressure -- saturation and relative permeability -- saturation curve

???The pressure distribution throughout the pore network as well as pore-network characteristics will also be written to a vtp output file that can be viewed with ParaView.???
Figure 1 shows the saturation distribution within the pore network for the case of flow in x-direction.

<figure>
    <center>
        <img src="img/result.png" alt="Pore-network pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - ???Pressure distribution within the pore network for flow in x-direction.??? </figcaption>
    </center>
</figure>


* TODO: add a figure for capillary pressure - saturation curve
* TODO: add a figure for relative permeability - saturation curve


After building the executable, run the simulation with `./example_pnm2p_upscaling`.

__Table of contents__. This description is structured as follows:

[[_TOC_]]

*TODO: review the documentation

## Problem setup

We consider a single-phase problem within a randomly generated pore network of 20x20x20 pores cubed from which some of the pore throats were [deleted randomly](https://doi.org/10.1029/2010WR010180).
The inscribed pore body radii follow a truncated log-normal distribution.
To calculate the upscaled properties, $`10`$ pressure differences in the range of $`10`$ to $`10^{10}`$ Pa/m are applied sequentially in directions specified by the user by setting the parameter `Problem.Directions` in the `params.input` file, while all lateral sides are closed.
The resulting mass flow rates are then used to determine upscaled properties as described [later](upscalinghelper.md).

## Mathematical and numerical model

In this example, we are using the single-phase pore-network model of DuMu<sup>x</sup>. We require mass conservation at each pore body $`i`$:

```math
 \sum_j Q_{ij} = 0,
```
where $`Q_{ij}`$ is the discrete volume flow rate in a throat connecting pore bodies $`i`$ and $`j`$. In case of creeping flow, the pore network model considers a Hagen-Poiseuille-type law to relate the volume flow from on pore body to another to discrete pressure drops $`\Delta p = p_i - p_j`$ between the pore bodies using throat conductance, $`g_{ij}`$.

```math
 Q_{ij} = g_{ij} (p_i - p_j),
```
or 

```math
 (p_i - p_j) = Q_{ij} / g_{ij} .
```

In the simulation of non-creeping flow, to capture inertial effects in fluid flow, an extension of the Hagen-Poiseuille-type law which includes the expansion and contraction of the pore space when moving from a throat to a pore and vice versa. 

```math
  (p_i - p_j) = Q_{ij} / g_{ij}  + (C_{exp} + C_{cont})Q^2,
```

where $`C_{exp}`$ and $`C_{cont}`$ are the expansion and the contraction coefficient, respectively. 

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
