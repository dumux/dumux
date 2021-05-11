The upscaling helper evaluates the pore-network simulation results for each direction $`i`$ and calculates the upscaled properties in this direction. Firstly, it evaluates the the Apparent velocity as:

```math
     v_{\mathrm{Apparent},i} = \frac{q_{\mathrm{mass,tot},i} / \varrho}{A_{\mathrm{tot},i}}
``` 

where $`q_{\mathrm{mass,tot},i}`$ is the total mass flow leaving the network over the REV's boundary with area
$`A_{\mathrm{tot},i}`$ in $`i`$-direction. $`\varrho `$ is the fluid mass density. Then, we calculate upscaled permeability as:

```math
 K_i = v_{\mathrm{Apparent},i} / \nabla p_i ~ \mu.
```
$`\nabla p_i`$ is a given pressure gradient in $`i`$-direction and $`\mu`$ the fluid dynamic viscosity. In creeping flow simulation, calculated permeability, $`K_i`$, is Darcy (intrinsic) permeability, $`K_D`$ of the system. 
To simulate non-creeping flow, we use Forchheimer's equation to upscale the properties.
```math
 \nabla p_i = \frac{\mu}{K_f} v_{\mathrm{Apparent},i} + \varrho \beta v_{\mathrm{Apparent},i}^2, 
```
where $`K_f`$ is Forchheimer permeability and $`\beta`$ is Forchheimer coefficient. To calculate upscaled properties, we rearrange Forchehimer's equation and find the linear regression line of $`\nabla p_i v_{\mathrm{Apparent},i}/\mu `$ versus $`\varrho v_{\mathrm{Apparent},i}/\mu `$. Using the intercept and the slope of the regreesion line, we can respectively calculate Forchheimr permeability and coefficient. We compute Darcy (intrinsic) permeability as the maximum permeability of the sample of data of the system which happens when pressure gradient is small enough such that inertial effect are negligible. Considering a slight difference between Darcy (intrinsic) permeability and Forchheimer permeability, in many applications they can be used interchangeabely. Here, we distinguish between them and calculate them separately.

The code documentation is structured as follows: