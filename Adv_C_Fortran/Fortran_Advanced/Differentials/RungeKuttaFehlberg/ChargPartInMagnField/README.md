# Charged Particle Motion in Earth's Magnetic Field (Fortran, RKF45)

This program simulates the trajectory of a charged particle (proton) moving in the Earth's dipole magnetic field by numerically integrating the equations of motion using the Runge-Kutta-Fehlberg 4(5) (RKF45) adaptive method.

---

## 1. **Physical Model**

The motion of a charged particle in a magnetic field is governed by the Lorentz force:

$$
m \frac{d\vec{v}}{dt} = q (\vec{v} \times \vec{B})
$$

where:
- $ m $: mass of the particle (kg)
- $ q $: charge of the particle (Coulombs)
- $ \vec{v} $: velocity vector (m/s)
- $ \vec{B} $: magnetic field vector (Tesla)

### **Earth's Magnetic Field**  
The Earth's field is approximated as a dipole aligned with the z-axis:

$$
\vec{B} = \frac{\mu_0}{4\pi} \frac{1}{r^3} [3(\vec{m} \cdot \hat{r})\hat{r} - \vec{m}]
$$

where:
- $ \vec{m} $: Earth's magnetic moment (A·m²), taken as $ (0,0,7.94 \times 10^{22}) $
- $ \mu_0/4\pi = 10^{-7} $ (SI units)
- $ r $: distance from Earth's center

---

## 2. **Differential Equations Being Solved**

The system is advanced as six first-order ODEs for position and velocity:

- **State vector:**  
  $ y = [x, y, z, v_x, v_y, v_z] $

- **ODEs:**  
  $$
  \begin{align*}
  \frac{dx}{dt} &= v_x \\
  \frac{dy}{dt} &= v_y \\
  \frac{dz}{dt} &= v_z \\
  \frac{dv_x}{dt} &= \frac{q}{m} (v_y B_z - v_z B_y) \\
  \frac{dv_y}{dt} &= \frac{q}{m} (v_z B_x - v_x B_z) \\
  \frac{dv_z}{dt} &= \frac{q}{m} (v_x B_y - v_y B_x) \\
  \end{align*}
  $$

  The right-hand side is implemented in the subroutine `func`.

---

## 3. **Program Structure and Variables**

### **Main Program: `charged_particle_motion`**

- **Constants:**
  - `neqn`: Number of equations (6)
  - `nt`: Number of time steps (e.g., 26000)
  - `dt`: Time step (s)
  - `relerr`, `abserr`: Relative and absolute error tolerances for RKF45

- **State/Time Variables:**
  - `y(6)`: State vector [x, y, z, vx, vy, vz]
  - `yp(6)`: Derivative vector (for RKF45 interface)
  - `t`: Current time (s)
  - `tout`: Target time for integration (s)
  - `flag`: RKF45 control flag (see below)

- **Physical Parameters:**
  - `q`: Particle charge (Coulombs), proton
  - `m`: Particle mass (kg), proton
  - `qm`: Charge-to-mass ratio (C/kg)
  - `R_maa`: Earth's mean radius (m)
  - `B0`: Not used in the current code, included for compatibility

- **Integration Loop:**
  - For each step, advance from `t` to `tout = t + dt` using `rkf45`.
  - Results are written to the file `proton.dat`.
  - Prints status to screen every 1000 steps.

---

### **Subroutines**

#### **`func(t, y, yp)`**
- Computes derivatives of the state vector.
- **Input:** 
  - `t`: Current time (unused)
  - `y`: Current state [x, y, z, vx, vy, vz]
- **Output:** 
  - `yp`: Derivatives [vx, vy, vz, dvx/dt, dvy/dt, dvz/dt]
- **Physics:** Implements the Lorentz force equation as above.

#### **`calc_magn_field(x, y, z, B)`**
- Computes Earth's dipole magnetic field at position (x, y, z).
- **Input:** 
  - `x, y, z`: Position in meters (SI)
- **Output:** 
  - `B(3)`: Magnetic field components [Bx, By, Bz] (Tesla)
- Uses Earth's dipole moment and the standard formula for a magnetic dipole field.

#### **`rkf45`**
- Advances the solution from `t` to `tout` for the system of ODEs.
- Adaptive step size control based on `relerr` and `abserr`.
- You must provide this routine (already in your project, not shown here).

---

## 4. **Input/Output**

- **Initial conditions** are set directly in the code for:
    - `y(1) = 3.0D7` : Start at 30,000 km from Earth's center on x-axis
    - `y(2) = 0.0D0`, `y(3) = 0.0D0`
    - `y(4) = 0.0D0`, `y(5) = 1.0D7`, `y(6) = 2.0D7` : Initial velocity (vx, vy, vz)

- **Output file:**  
  - `proton.dat`: Columns are time, position (x,y,z), velocity (vx,vy,vz); units are SI (meters, m/s).

- **Console output:**  
  - Prints step number and current state every 1000 steps.

---

## 5. **RKF45 Usage**

- **Flag control:**  
  - At the start, `flag = 1` (integrate to `tout` exactly).  
  - After the first step, `flag = 2` (let RKF45 manage step size, not guaranteed to hit `tout` exactly).
  - If you want output at *precise* intervals (`dt`), always reset `flag = 1` before each call.

- **Error handling:**  
  - If `flag` returns negative or 6/8, the integration stops.

---

## 6. **Units**

- **Distance:** meters (m)
- **Velocity:** meters/second (m/s)
- **Time:** seconds (s)
- **Magnetic field:** Tesla (T)

---

## 7. **File List**

- `charged_particle_motion.f90` : Main program and ODE definitions
- `rkf45.f90` : Runge-Kutta-Fehlberg ODE integrator (must be provided!)

---

## 8. **References**

- [Runge-Kutta-Fehlberg method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
- [Dipole magnetic field theory](https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field)
- Fortran documentation and syntax

---

## 9. **Customization**

- To simulate a different particle, change `q`, `m`, and initial conditions.
- To start at a different point, change `y(1:3)`.
- To simulate longer or shorter, adjust `nt` and `dt`.

---

**Short summary:**  
This code simulates the 3D motion of a proton in Earth's dipole field using high-accuracy adaptive ODE integration. The code is modular and easy to adapt for other particles or field models.
