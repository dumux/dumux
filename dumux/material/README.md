# Derivation of Fugacities and Fugacity Coefficients from Henry's Constant for Air-Water System

## 1. Basic Assumptions and Definitions

We start with the assumption of equilibrium between liquid and vapor phases for a given component:

$f_i^L = f_i^V$

This can be written using the approach:

$\phi_i x_i p^L = p^V y_i = p_i^V$

where:
- $\phi_i$ is the fugacity coefficient
- $x_i$ is the mole fraction in the liquid phase
- $y_i$ is the mole fraction in the vapor phase
- $p$ is the total pressure
- $p_i^V$ is the partial pressure of component $i$ in the vapor phase

## 2. Henry's Constant

The Henry's coefficient is often defined as:

$K_H = \frac{y_i p^V}{x_i}=K_Dp^V$

## 3. Relating Fugacity Coefficient to Henry's Constant

From the definitions above, we can relate the fugacity coefficient to Henry's constant:

$\phi_i x_i p^L = p^V y_i$

Dividing both sides by $x_i$:

$\phi_i  p^L = \frac{y_i}{x_i}p^V$

Substituting the definition of Henry's constant:

$\phi_i p^L= K_H$
$\phi_i=\frac{K_H}{p^L}$

This provides an ansatz for fugacity irrespective of the equilibrium. Many data sources for vapor-liquid equilibrium (VLE) use these definitions.
For different definitions of the Henry coefficient one needs to be carefull. For all water-something the IAPWS calculation for $K_D$ times vapor pressure is used.
## 4. Application to Air-Water System

For an air-water system, we consider the main components of air: N₂, O₂, Ar, and CO₂. The rest are ignored for simplicity. The mole fractions of these components in the lumped component air are:

- O₂: 0.209
- Ar: 0.00934
- CO₂: 0.0004
- N₂: 1 - (0.209 + 0.00934 + 0.0004) = 0.78126

We treat N₂ as the balance to ensure the concentrations sum to 1:


## 5. Fugacity Equilibria Equations

Writing down the fugacity equilibria for all components:

$\phi_{N_2} x_{N_2} p^L = p^V y_{N_2}$

$\phi_{O_2} x_{O_2} p^L = p^V y_{O_2}$

$\phi_{Ar} x_{Ar} p^L = p^V y_{Ar}$

$\phi_{CO_2} x_{CO_2} p^L = p^V y_{CO_2}$

## 6. Derivation of $\phi_{air-water}$ and $K_{H,air}$

Let's sum up all equations and introduce $x_{air}$ as the sum of all concentrations in the liquid phase and $y_{air}$ as the sum in the gas phase:

Furthermore, we split up

$y_i=y_{i,air} y_{air}$

yielding in:

$\sum_i x_i=x_{air}=\frac{p^V y_{air}}{p^L}(\sum_i \frac{y_{i,air}}{\phi_i})$

Rearranging and renaming:

$(\sum_i \frac{y_{i,air}}{\phi_i})^{-1}=\phi_{air}$

Yields:

$\phi_{air} x_{air} p^L=p^V y_{air}$

But note, that the composition of air in water is different in comparison to gas.
Similary, one can derive for $K_{h,air}$:
From:

$(\sum_i \frac{y_{i,air}}{\phi_i})^{-1}=\phi_{air}$

with:
$\phi_i=K_{H,i}/p^L$

$K_{H,air}=(\sum_i \frac{y_{i,air}}{K_{H,i}})^{-1}$