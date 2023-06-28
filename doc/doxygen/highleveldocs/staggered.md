# Staggered Navier-Stokes Discretization

## Conservation of Mass

The PDE describing the conservation of Mass can be written in the following way:
\f{equation}{
 \frac{\partial \rho}{\partial t}
- \nabla \cdot (\rho\boldsymbol{u})
- q
= 0
\f}

In DuMuX the used grid to dicretize the equation is the following:

<div align="center">
  <img src="staggered_mass.svg"  width="50%">
</div>

Describing the equation is almost straight-forward. One just needs to think about upwinding in the flux Term:

\f{equation}{
 \frac{\partial \rho}{\partial t}
= \frac{\rho^{t+1}_C-\rho^{t}_C}{\Delta t}
\f}

For the advective flux
\f{equation}{
(\rho\boldsymbol{u})
\f}
One needs can describe for example the flux on a north face:
\f[\rho v_N \f]
For the east facing slope:
\f[\rho v_E \f]
Again one needs to think about upwinding the density $ \rho $.

## Conservation of linear Momentum in the y- Direction

\f{equation}{
    \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}})) - \nabla p + \varrho \textbf{g} - \textbf{f}
\f}

The grid to discetize this equation is now staggered. As seen in the overview.

<div align="center">
  <img src="staggered_y_dir.svg"  width="50%">
</div>

Describing the storage term yields:

\f{equation}{
    \frac{\partial (\varrho \textbf{v})}{\partial t}
    = \frac{(\rho \cdot u_C)^{t+1} -(\rho \cdot u_C)^{t}}{\Delta t}
\f}
The advective flux
\f{equation}{
(\varrho \textbf{v} \textbf{v}^{\text{T}})
\f}
For a north facing slope to discretize the y-component of the momomentum, one needs to discretize:
\f[\rho v \cdot v \f]
This can be described:
\f[\rho (\frac{v_N + v_C}{2})^2 \f]

The pressure term in the y-Direction is the following:

\f{equation}{
   - \nabla p = - \frac{p_N - p_S}{\Delta y}
\f}

To evaluate the stress tensor one needs to evaluate the stress tensor at the faces.

\f{equation}{
  (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}})) = \mu \begin{bmatrix}
   2\frac{\partial u}{\partial x} &
   \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x} \\
   \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x} &
   2\frac{\partial v}{\partial y}
   \end{bmatrix} \\
\f}
Evaluating this for the north face only leaves one term:
\f[ 2\frac{\partial v}{\partial y} \f]
which can be described in the form:
 \f[ \frac{v_N - v_C}{L} \f]

Evaluating this for the east face also leaves one with one term:
\f[ \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x} \f]
which can be described by:
\f[\frac{u_{NE}-u_{SE}}{\Delta y} + \frac{v_E - v_C}{L}\f]
