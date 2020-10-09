# Part 3: Upscaling helper

The upscaling helper evaluates the pore-network simulation results for each direction $`i`$ and calculates the upscaled intrinsic permeability in this direction using:

```math
 K_i = v_{\mathrm{Darcy},i} / \nabla p_i ~ \mu.
```
$`\nabla p_i`$ is a given pressure gradient in $`i`$-direction and $`\mu`$ the fluid dynamic viscosity.

We evaluate the Darcy velocity as

```math
     v_{\mathrm{Darcy},i} = \frac{q_{\mathrm{mass,tot},i} / \varrho}{A_{\mathrm{tot},i}}
```

where $`q_{\mathrm{mass,tot},i}`$ is the total mass flow leaving the network over the REV's boundary with area
$`A_{\mathrm{tot},i}`$ in $`i`$-direction. $`\varrho `$ is the fluid mass density.

The code documentation is structured as follows:

[[_TOC_]]
