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
where $`K_f`$ is Forchheimer permeability and $`\beta`$ is Forchheimer coefficient. Although some researchers for the sake of simplicity assumes that $`K_f = K_i`$, they are not exactly the same properties. As the velocity increases, the flow regime in a porous medium shifts from Darcy to Forchheimer regime. This change in the flow regime causes that the pressure drop ,which in Darcy flow just includes the viscous dissipation, becomes a combination of both viscous dissipation (the first term in the Forchheimer equation) and the inertial term (the second term in Forchheimer equation). Considering the shift in flow regime, a porous medium having a Forchheimer flow regime shows a different viscous dissipation than when the same porous medium experiences a Darcy fllow. In other words, moving from Darcy to Forchheimr regime establishes a new velocity field in the porous medium which needs a new viscous dissipation and also an inertial term. Furthermore, the first and second terms in Forchheimer equation have strong influence on each other. For more detail, We refer to the study conducted by [Dukhan and Minjeur (2010)](https://link.springer.com/article/10.1007/s10934-010-9393-1). To calculate upscaled properties, we rearrange Forchehimer's equation and find the linear regression line of $`\nabla p_i v_{\mathrm{Apparent},i}/\mu `$ versus $`\varrho v_{\mathrm{Apparent},i}/\mu `$. Using the intercept and the slope of the regreesion line, we can respectively calculate Forchheimr permeability and coefficient. It should be noted that the calculation of the Forchheimer permeability is highly affected by the pressure range applied to the porous medium as well as the number of sample points that are used in the regression process. We compute Darcy (intrinsic) permeability as the maximum permeability of the sample of data of the system which happens when pressure gradient is small enough such that inertial effects are negligible. To ensure such a small pressure gradient, it is recommended to use more than 10 pressure smaple points which can be set in the input file. As mentioned before, considering a slight difference between Darcy (intrinsic) permeability and Forchheimer permeability, in many applications they can be used interchangeabely. Here, however, we distinguish between them, calculate and report them separately.

The code documentation is structured as follows: