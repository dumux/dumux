# Part 3: Upscaling helper

The upscaling helper evaluates for each direction the pore-network simulation results and calculates the upscaled instrinsic permeability in this direction using:

```math
 K = v_\mathrm{Darcy} / \nabla p ~ \mu.
```
$`\nabla p`$ is a given pressure gradient and $`\mu`$ the fluid dynamic viscosity.

We evaluate the Darcy velocity as

```math
     v_\mathrm{Darcy} = \frac{q_\mathrm{mass,tot} / \varrho}{A_\mathrm{tot}}
```

where $`q_\mathrm{mass,tot}`$ is the total mass flow leaving the network over the REV's boundary with area
$`A_\mathrm{tot}`$. $`\varrho `$ is the fluid mass density.


The code documentation is structured as follows:

[[_TOC_]]
