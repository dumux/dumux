#Shallow water flow with bottom friction
This example shows how the shallow water flow model can be
applied to simulate steady subcritical flow including
bottom friction (bed shear stress).


## Shallow water model
The shallow water equations (SWEs) are given as:

$$
\frac{\partial \mathbf{U}}{\partial t} +
\frac{\partial \mathbf{F}}{\partial x} +
\frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
$$

where $U$, $F$ and $G$  defined as

$$
\mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
\mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 \\ huv \end{bmatrix},
\mathbf{G} = \begin{bmatrix} hv \\ huv \\ hv^2  + \frac{1}{2} gh^2 \end{bmatrix}
$$

Z is the bedSurface, h the water depth, u the velocity in
x-direction and v the velocity in y-direction, g is the constant of gravity.

The source terms for the bed friction $$S_b$$ and bed slope
$$S_f$$ are given as
$$
\mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
               \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
\mathbf{S_f} = \begin{bmatrix} 0 \\ -ghS_{fx} \\ -ghS_{fy}\end{bmatrix}.
$$

For this example, a cell-centered finite volume method (cctpfa) is applied to solve the SWEs
in combination with a fully-implicit time discretization. For cases where no sharp fronts or
traveling waves occur it is possible to apply time steps larger than CFL number = 1 to reduce
the computation time. Even a steady state solution is considered an implicit time stepping method
is applied.

##Problem set-up
The model domain is given by a rough channel with a slope of 0.001. Bottom friction is considered by applying
the friction law of Manning (Manning n = 0.025). At the lateral sides no friction is considered and  a
no-flow no slip boundary condition is applied. This is the default boundary condition for the shallow water model.

The domain is 1000 meters long and 10 meters wide. At the left border a discharge boundary condition
is applied as inflow boundary condition with q = -1.0 ($m^2 s^{-1}. At the right border a water fixed depth boundary condition
is applied for the outflow. Normal flow is assumed, therefor the water depth at the right border is calculated after
the of Gaukler-Manning-Strickler equation:

 $$ v_m = 1/n * R_{hy}^{2/3} * I_s^{1/2}$$

Where the mean velocity $v_m$ is given as

$$ v_m = \frac{q}{h}$$,

$n$ is the friction value after Manning. $R_{hy} the hydraulic radius, which is assumed to be equal to
the water depth. $I_s$ is the bed slope and $q$ the unity inflow discharge

The water depth h can be calculated as
$$h = \left(\frac{n*q}{\sqrt{I_s}} \right)^{3/5}$$

The formula of Gaukler Manning and Strickler is also used to calculate the analytic solution. All parameters
for the simulation are given in the file *params.input*.









