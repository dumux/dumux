> :footprints: This documented example was originally contributed by: [Martin Schneider](https://www.iws.uni-stuttgart.de/en/institute/team/Schneider/), Dennis Gläser, and [Timo Koch](https://timokoch.github.io/).

# Rotation-symmetric pressure distribution

__In this example, you will learn how to__

* solve a rotation-symmetric problem one-dimensionally
* perform a convergence test against an analytical solution
* apply the `Rotational Extrusion` filters in [ParaView](https://www.paraview.org/) for a two-dimensional visualization of the one-dimensional results


__Result__. With the `Rotational Extrusion` and the `Warp By Scalar` filters in [ParaView](https://www.paraview.org/),
the pressure distribution of this example looks as shown in the following picture:

<figure>
    <center>
        <img src="img/result.png" alt="Rotation-symmetric pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Rotation-symmetric pressure distribution on a disc (warped to 3D). </figcaption>
    </center>
</figure>


__Table of contents__. This description is structured as follows:

[[_TOC_]]


## Problem setup

We consider a single-phase problem that leads to a rotation-symmetric pressure distribution.
The following figure illustrates the setup:

<figure>
    <center>
        <img src="img/setup.svg" alt="Rotation-symmetric setup" width="60%"/>
        <figcaption> <b> Fig.2 </b> - Setup for the rotation-symmetric problem. The pressure boundary conditions are shown by the colored lines and the simulation domain is depicted in grey.</figcaption>
    </center>
</figure>

This could, for example, represent a cross section of an injection/extraction well in a homogeneous
and isotropic porous medium, where the well with radius $`r_1`$ is cut out and the
injection/extraction pressure $`p_1`$ is prescribed as a Dirichlet boundary condition. At the outer
radius $`r_2`$, we set the pressure $`p_2`$. In the polar coordinates $`r`$ and $`\varphi`$, the
solution to this problem is independent of the angular coordinate $`\varphi`$ and can be reduced to
a one-dimensional problem in the radial coordinate $`r`$. Therefore, in this example, we want to
solve the problem on a one-dimensional computational domain as illustrated by the orange line in
the above figure.

## Mathematical model

In this example we are using the single-phase model of DuMu<sup>x</sup>, which considers Darcy's law to relate
the Darcy velocity $`\textbf u`$ to gradients of the pressure $`p`$. In the case of rotational
symmetry, the mass balance equation for the fluid phase can be transformed using polar coordinates:

```math
-\frac{1}{r} \frac{\partial}{\partial r} \left( r  \frac{\varrho k}{\mu} \frac{\partial p}{\partial r} \right) = 0,
```

where we identify the Darcy velocity in radial direction $`u_r = -\frac{k}{\mu} \frac{\partial p}{\partial r}`$,
and where $`k`$ is the permeability of the porous medium, $`\mu`$ is the dynamic viscosity of the
fluid, and $`\varrho`$ is the fluid density.

## Discretization

We employ a finite-volume scheme to spatially discretize the mass balance equation shown above.
Let us consider a discretization of the one-dimensional domain into control volumes
$`K_i = \left[ r_i, r_{i+1} \right]`$. The discrete equation describing mass conservation inside a control volume
$`K_i`$ is obtained by integration and reads:

```math
    - 2 \pi r_{i+1} \left( \varrho u_r \right)_{r_{i+1}}
    + 2 \pi r_i \left( \varrho u_r \right)_{r_i}
    = 0.
```

For this type of equation, the implementation of the finite-volume schemes in DuMu<sup>x</sup> is based on
the general form:

```math
\sum_{\sigma \in \mathcal{S}_K} | \sigma | \left( \varrho \textbf u \cdot \textbf n \right)_\sigma = 0,
```

where $`\sigma`$ are the faces of the control volume and where the notation
$`( \cdot )_\sigma`$ was used to denote quantities evaluated for a face $`\sigma`$.
The area of a face is denoted with $`| \sigma |`$. Thus, comparing the two equations
we identify $`| \sigma | = 2 \pi r_\sigma`$ for the case of rotational symmetry
on a disc. Here, $`r_\sigma`$ refers to the radius at which the face is situated
in the one-dimensional discretization.

In DuMu<sup>x</sup>, this is rotational extrusion is approximated by using modified control volume
volumes and control volume face areas for the mid-point integration rule.

# Implementation & Post processing
