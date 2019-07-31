#Shallow water flow with bottom friction
This example shows how the shallow water flow model of Dumux can be
applied to simulate steady subcritical flow including bottom friction (bed shear stress).

## Shallow water model
The shallow water equations (SWEs) are given as:

\f[
\frac{\partial \mathbf{U}}{\partial t} +
\frac{\partial \mathbf{F}}{\partial x} + \\
\frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
\f]

with  <B>U</B>, <B>F</B>, <B>G</B> are defined as

\f[
\mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
\mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 \\ huv \end{bmatrix},
\mathbf{G} = \begin{bmatrix} hv \\ huv \\ hv^2  + \frac{1}{2} gh^2 \end{bmatrix}
\f]

Z is the bedSurface, h the water depth, u the velocity in
x-direction and v the velocity in y-direction, g is the constant of gravity.

The source terms for the bed friction <B>S_b</B> and bed slope
<B>S_f</B> are given as
\f[
\mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
               \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
\mathbf{S_f} = \begin{bmatrix} 0 \\ -ghS_{fx} \\ -ghS_{fy}\end{bmatrix}.
\f]

For this example, a cell-centered finte volume method (cctpfa) is applied to solve the SWEs
in combination with a fully-implicit time discretization. For large time
 * step sizes (CFL > 1) this can lead to a strong semearing of sharp fronts.
 * This can be seen in the movement of short traveling waves (e.g. dam break
 * waves). Nevertheless the fully implicit time discretization showes oftenFor large time
 * step sizes (CFL > 1) this can lead to a strong semearing of sharp fronts.

 * no drawbacks in cases where no short waves are considered. Thus the model
 * can be a good choice for simulating flow in rivers and channels, where the
 * fully-implicit discretization allows large time steps and reduces the
 * overall computation time drastically.

##Problem set-up
