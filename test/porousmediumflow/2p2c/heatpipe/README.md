# Heatpipe Effect {#benchmark-heatpipe}

## One-dimensional non-isothermal two-phase two-component flow

**Problem Description**

The Heatpipe Effect is described by a nonisothermal water-gas system in a porous medium, in which the heat transfer processes convection, conduction, and diffusion as well as capillary forces play an essential role. Udell and Fitch (1985) \cite Udell:1985 provide a semi-analytical solution for this system, which is practical for comparing with numerical results (e.g., see  Emmert, 1997 \cite emmertpromo).

A one-dimensional horizontal porous column is considered. A constant heat flux is applied at the right boundary. Due to the heat flux, the system is heated until boiling temperature is reached and steam is produced at the right-hand boundary. This causes a pressure gradient in the gas phase and the steam flows away from the heat source. After reaching cooler regions of the column, the steam condenses and sets free its latent heat of vaporization. After a while, a non-uniform saturation profile is obtained with a gradient from the cooler to the hot end of the heatpipe.

According to the capillary pressure–saturation relationship a gradient of the capillary pressure into the same direction is produced. Hence, the pressure gradients of the phases have opposite directions and a circulation flow is created. After a stationary system state has been reached, three regions can be distinguished, each of them associated with a dominant heat-transport process.

![Schematic description](heatpipe-schematic-description.png){html: width=80%}


**Semi-analytical Reference Solution**

Udell and Fitch (1985) \cite Udell:1985 derive four coupled first-order differential equations for pressure, saturation, temperature and gas-phase mole fraction. These equations are solved by numerical integration by means of a fourth-order Runge–Kutta method. The numerical simulation of the heatpipe system was carried out with the BOX discretization method. Note that the choice of BOX or CVFE makes no difference in the present one-dimensional case. For the simulations, the upwind-technique for weighting the mobilities and mole fractions at the interface between two boxes has been switched off since upwinding cannot be considered in the analytical solution.


**Setup**

A one-dimensional horizontal porous column is considered:

A constant heat flux of $q = 100\ \mathrm{W}$ is imposed at the right boundary of the horizontal column. Zero-flux (Neumann) boundary conditions are prescribed for all mass components. The initial water saturation in the entire domain is $S_w = 0.5$.

At the left boundary, Dirichlet boundary conditions are applied for the gas-phase pressure $p_g = 101330\ \mathrm{Pa}$, the effective water saturation $S_{we} = 1.0$ and the air mole fraction in the gas phase $x_g^a = 0.71$.

The corresponding equilibrium temperature is determined from $S_{we}$, $p_g$, and $x_g^a$, yielding $T = 68.6^\circ\mathrm{C}$.

The temperature is computed according to

```math
T =
\frac{
T_0 \left(
1 + \dfrac{p_c - x_g^a p_g}{h_v^w \rho_w}
\right)
}{
1 - T_0 \dfrac{R}{h_v^w M^w}
\ln\left(
\dfrac{p_g (1 - x_g^a)}{p_0}
\right)
}.
```

The following model parameters were used for the simulation run:

| Parameter                                    | Symbol    | Value                   | Unit  |
|----------------------------------------------|-----------|-------------------------|-------|
| Permeability                                 | $K$       | $1.0\text{e-}12$        | m²    |
| Porosity                                     | $\phi$    | $0.4$                   | -     |
| Residual wetting-phase saturation            | $S_{wr}$  | $0.15$                   | -     |
| Heat conductivity of the fully saturated porous medium | $\lambda_{pm}^{S_w=1}$ | $1.13$ | W/(m*K) |
| Heat conductivity of the dry porous medium | $\lambda_{pm}^{S_w=0}$ | $0.582$     | W/(m*K) |
| Soil grain density                           | $\varrho_s$   | $2600$              | kg/m³  |
| Specific heat capacity of the soil grains    | $c_s$     | $700$                   | J/(kg*K) |
| Density of water                             | $\varrho_w$  | $958.4$              | kg/m³ |
| Dynamic viscosity of water                   | $\mu_w$  | $2.938\text{e-}4$        | Pa*s |
| Dynamic viscosity of air                     | $\mu_g^a$  | $2.08\text{e-}5$       | Pa*s |
| Dynamic viscosity of steam                   | $\mu_g^w$  | $1.2\text{e-}5$        | Pa*s |

A function according to  Fatt and Klikoff (1959) \cite Fatt:1959 is chosen for the relative permeability-saturation relationship:
```math
k_{rg} = (1 - S_e)^3 \quad \text{for steam (gas phase)} \nonumber 
```
```math
k_{rw} = S_e^3 \quad \text{for water} \, 
```
with the effective water-phase saturation
```math
S_e = \frac{S_w - S_{wr}}{1-S_{wr}}\; .
```
For the capillary pressure-saturation relationship, the following function of Leverett (1941) \cite lev1 is used:
```math
p_c = p_0 \cdot \gamma \cdot 1.417(1-S_e) - 2.120(1-S_e)^2 + 1.263(1-S_e)^3 \, .
```

The surface tension $\gamma$ at $T = 100.5$$^\circ$C is 0.05878 Nm$^{-1}$ and $p_0 = \sqrt{\phi/K}$ applies for the scaling pressure. A constant value of 0.5 is assigned to the tortuosity $\tau$ and the binary diffusion constant $D_g^{aw}$ takes the value $2.6 \cdot 10^{-6}$ m$^2$/s.

The dimension of the model domain in x--direction is chosen at 2.4 m. However, this is not important for the length of the heatpipe after the stationary state has been reached as long as the domain is sufficicently large for the heatpipe to be built. We used a discretization length of $\Delta x = 0.04$ m.

The initial conditions were chosen to $p_g = 101330$ Pa, $S_w = 0.5$ und $T = 70^\circ$C.


**Result**

To run the test and produce plots below, execute:
```bash
python3 compile_run_plot.py
```
The script expects PyVista and Matplotlib to be available for post-processing.

The results of the numerical simulation with the semi-analytical solution are compared below. In all cases, the curves match very well. A length of $\approx$ 2.0m is obtained consistently for the heatpipe. 

Due to the complex interaction of different physical processes and the excellent agreement between the simulation and the semi-analytical solution, we state that the verification of the numerical model for an air--water system was 
successful.

With the given initial conditions, the model could reproduce the heating at the right-hand boundary from 70$^\circ$C up to boiling temperature. Also, the gradual extension of the heatpipe region until the stationary system state was reached was modeled correctly. The disappearance of the water phase associated with a change of the phase state and a substitution of the primary variables was carried out well. Fig. 4 shows clearly that the pressure gradients of both phases have different signs.


**References**

M. Emmert. Numerische Simulation von isothermen/ nichtisothermen Mehrphasen-prozessen unter Berücksichtigung der Veränderung der Fluideigenschaften. PhD thesis, Institut für Wasserbau, Universität Stuttgart, 1997.

I. Fatt and W.A. Klikoff. Effect of fractional wettability on multiphase flow through porous media. AIME Transactions, 216:246, 1959.

M.C. Leverett. Capillary Behavior in Porous Solids, volume 142. AIME Transactions, 1941.

K.S. Udell and J.S. Fitch, editors. Heat and mass transfer in capillary porous media considering evaporation, condensation and non-condensible gas effects, Denver, CO, August 1985. Paper presented at 23rd ASME/AIChE National Heat Transfer Conference.