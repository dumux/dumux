The upscaling helper evaluates the pore-network simulation results for each direction $`i`$ and calculates the upscaled properties in this direction. Firstly, it evaluates the apparent velocity as:

```math
     v_{\mathrm{Apparent},i} = \frac{q_{\mathrm{mass,tot},i} / \varrho}{A_{\mathrm{tot},i}}
``` 

where $`q_{\mathrm{mass,tot},i}`$ is the total mass flow leaving the network over the REV's boundary with area
$`A_{\mathrm{tot},i}`$ in $`i`$-direction. $`\varrho `$ is the fluid mass density. Then, we calculate upscaled permeability as:

```math
 K_D = v_{\mathrm{Apparent},i}  ~ \mu / \nabla p_i.
```
$`\nabla p_i`$ is a given pressure gradient in $`i`$-direction and $`\mu`$ the fluid dynamic viscosity. In the creeping flow simulation, the calculated permeability, $`K_D`$, is the Darcy (intrinsic) permeability of the system.

To simulate non-creeping flow, we use Forchheimer's equation to upscale the properties. 
```math
 \nabla p_i = \frac{\mu}{K_f} v_{\mathrm{Apparent},i} + \varrho \beta v_{\mathrm{Apparent},i}^2, 
```
where $`K_f`$ is the Forchheimer permeability and $`\beta`$ is the Forchheimer coefficient.

Although it is sometimes assumed for the sake of simplicity that $`K_f = K_D`$, they are not exactly the same properties. As the velocity increases, the flow regime in a porous medium shifts from the Darcy to the Forchheimer regime. This change in the flow regime means that in addition to viscous dissipation (the first term in the Forchheimer equation), we also need to consider inertia effects (the second term in the Forchheimer equation). Moreover, the shift in the flow regime, means a different viscous dissipation than what the same porous medium experiences a Darcy flow. In other words, moving from the Darcy to the Forchheimer regime establishes a new velocity field in the porous medium which causes the difference in viscous dissipation in addition to inertia effects. The first and second terms in the Forchheimer equation are correlated. For more detail, we refer to the study conducted by [Dukhan and Minjeur (2010)](http://dx.doi.org/10.1007/s10934-010-9393-1).

In this example the end of the creeping flow (Darcy regime) is defined as the moment when the pressure drop calculated by the Dracy (intrinsic) permeability $`\Delta p_i = v_{\mathrm{Apparent},i} \mu l/ K_D`$ becomes less than 99% of the total pressure drop [Muljadi et al.](https://doi.org/10.1016/j.advwatres.2015.05.019).
To calculate upscaled properties, we rearrange Forchehimer's equation as:

```math
 \frac{\nabla p_i v_{\mathrm{Apparent},i}}{\mu} = \frac{1}{K_f} + \frac{\varrho v_{\mathrm{Apparent},i}}{\mu} \beta .
```
Finding the linear regression line of $`\nabla p_i v_{\mathrm{Apparent},i}/\mu `$ versus $`\varrho v_{\mathrm{Apparent},i}/\mu `$ and using the intercept and the slope of the regression line, we can respectively calculate the Forchheimer permeability and coefficient. It should be noted that the calculation of the Forchheimer permeability can be affected by the pressure range applied to the porous medium as well as the number of sample points that are used in the regression process. We compute the Darcy (intrinsic) permeability as the maximum permeability of the sample of data of the system which happens when the pressure gradient is small enough such that inertial effects are negligible. To ensure such a small pressure gradient, we set the first pressure gradient to be applied as $'10 Pa/m'$. However, this value can be adapted using the keyword `Problem.MinimumPressureGradient` in params.input. it is recommended to use more than 10 pressure sample points which can be set in the input file. As mentioned before, considering a slight difference between Darcy (intrinsic) permeability and Forchheimer permeability, in many applications they can be used interchangeabely. Here, however, we distinguish between them, calculate and report them separately.

The code documentation is structured as follows:
