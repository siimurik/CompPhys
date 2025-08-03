# 📘 README: Orbital Dynamics Jacobian for Stiffness Detection

This project analyzes the stiffness of an orbital dynamics system by constructing the Jacobian matrix of a system of ordinary differential equations (ODEs) and calculating eigenvalues using LAPACK with OpenBLAS.

## System Overview

We consider a simplified 2D gravitational model:

### Variables
- `y(1)` → x position
- `y(2)` → y position
- `y(3)` → x velocity (vx)
- `y(4)` → y velocity (vy)

### Equations of Motion
```
ydot(1) = y(3)                         ! dx/dt = vx
ydot(2) = y(4)                         ! dy/dt = vy
ydot(3) = -GM * y(1) / r³              ! dvx/dt = -GM*x/r³
ydot(4) = -GM * y(2) / r³              ! dvy/dt = -GM*y/r³
```

Where:
- `GM` is the gravitational constant multiplied by mass
- `r² = x² + y²`
- `r³ = (r²)^(3/2)`
- `r⁵ = (r²)^(5/2)`

---

## Jacobian Matrix

We define the Jacobian `J(i,j)` as the partial derivative of `ydot(i)` with respect to `y(j)`. The matrix is 4×4:

$$
J = 
\begin{bmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
\frac{\partial \dot{v}_x}{\partial x} & \frac{\partial \dot{v}_x}{\partial y} & 0 & 0 \\
\frac{\partial \dot{v}_y}{\partial x} & \frac{\partial \dot{v}_y}{\partial y} & 0 & 0
\end{bmatrix}
$$

### Partial Derivatives (Simplified)
```fortran
r2 = y(1)**2 + y(2)**2
r3 = r2**1.5
r5 = r2**2.5

J(3,1) = -GM * (1.0d0 / r3 - 3.0d0 * y(1)**2 / r5)
J(3,2) =  3.0d0 * GM * y(1) * y(2) / r5
J(4,1) = J(3,2)
J(4,2) = -GM * (1.0d0 / r3 - 3.0d0 * y(2)**2 / r5)
```

All other entries are 0 except those related to velocity as described.

---

## Usage

You can plug the Jacobian into your stiffness evaluation code to analyze system behavior using LAPACK’s `dgeev` routine.

## Notes

- Ensure `GM` is properly defined for your simulation.
- The Jacobian is derived analytically for performance and stability.
- Use this inside a `compute_jacobian(y, J)` subroutine for modularity.

---

Great question—it gets right to the heart of how differential systems are analyzed mathematically. The structure of the Jacobian matrix isn’t arbitrary; it’s rooted in how we describe the sensitivity of a system of equations to its variables.

---

## What Is the Jacobian?

In a system of **first-order differential equations** like:

$$
\dot{\mathbf{y}} = \mathbf{f}(\mathbf{y}, t)
$$

where:
- $\mathbf{y} = [y_1, y_2, ..., y_n]^T$ is the state vector,
- $\mathbf{f}(\mathbf{y}, t)$ defines the dynamics,

the **Jacobian matrix** is the matrix of **first-order partial derivatives**:

$$
J_{ij} = \frac{\partial f_i}{\partial y_j}
$$

This means each element $J_{ij}$ measures how the $i^\text{th}$ equation responds to a change in the $j^\text{th}$ variable—basically, the **local sensitivity** of your system.

---

## Why That Specific Structure for Your ODEs?

Your system involves gravitational dynamics in 2D space:

- $\dot{x} = v_x$
- $\dot{y} = v_y$
- $\dot{v}_x = -GM \cdot x / r^3$
- $\dot{v}_y = -GM \cdot y / r^3$

Since velocity directly affects position, the upper part of the Jacobian has ones in columns 3 and 4. Meanwhile, acceleration ($\dot{v}_x, \dot{v}_y$) depends **nonlinearly** on position ($x, y$), so those rows contain **derivatives of inverse cubic expressions**, making their structure more complex. That’s why the lower-left block is filled and the lower-right block is zero—because velocity does not directly influence acceleration in Newtonian gravity.

---

## Intuition Behind Layout

It's helpful to think of the Jacobian as a block matrix:

$$
J =
\begin{bmatrix}
0 & I \\
A(\mathbf{y}) & 0
\end{bmatrix}
$$

- Top left: zeros (position equations don’t depend on position)
- Top right: identity (position depends directly on velocity)
- Bottom left: nonlinear derivatives (acceleration depends on position)
- Bottom right: zeros (acceleration doesn’t depend on velocity)

---

Nice setup. You're dealing with a classic low-pass RLC filter driven by a source voltage $ V_b(t) $, and you're interested in the currents through the inductor ($ i_L $), resistor ($ i_R $), and capacitor ($ i_C $). Let’s derive the system into something you could simulate or analyze numerically.

---

## 🔁 System Equations

You're given:

1. $ V = i_R R $
2. $ i_C = C \dfrac{dV}{dt} $
3. $ L \dfrac{di_L}{dt} = V_b - V $
4. $ i_L = i_R + i_C $

Let’s rearrange this into a state-space system. Let’s choose **$ V $** and **$ i_L $** as the primary state variables:

---

## 📐 Define State Variables

Let:
- $ x_1 = V $
- $ x_2 = i_L $

Then:

- From (2): $ i_C = C \dfrac{dx_1}{dt} $
- From (1): $ i_R = \dfrac{x_1}{R} $
- From (4): $ x_2 = \dfrac{x_1}{R} + C \dfrac{dx_1}{dt} $

We now solve for $ \dfrac{dx_1}{dt} $:

### 🔄 Step 1: Express $ \dfrac{dx_1}{dt} $
$$
x_2 = \dfrac{x_1}{R} + C \dfrac{dx_1}{dt} \Rightarrow \dfrac{dx_1}{dt} = \dfrac{x_2 - x_1 / R}{C}
$$

### 🔄 Step 2: From (3) express $ \dfrac{dx_2}{dt} $
$$
L \dfrac{dx_2}{dt} = V_b - x_1 \Rightarrow \dfrac{dx_2}{dt} = \dfrac{V_b - x_1}{L}
$$

---

## 🧮 Final System (Two ODEs)

$$
\boxed{
\begin{aligned}
\dfrac{dx_1}{dt} &= \dfrac{x_2 - x_1 / R}{C} \\
\dfrac{dx_2}{dt} &= \dfrac{V_b(t) - x_1}{L}
\end{aligned}
}
$$

Once you solve this system numerically (e.g. with Runge-Kutta), you can recover:
- $ i_R = \dfrac{x_1}{R} $
- $ i_C = C \dfrac{dx_1}{dt} $
- $ i_L = x_2 $

---

🧮 Great question! To check for stiffness in your RLC system, you’d construct the **Jacobian matrix** of the ODE system you’re solving. Let’s break this down.

---

## 🔧 **The System of Equations**

You’re tracking:
- $ x_1 = V $ → capacitor voltage
- $ x_2 = i_L $ → inductor current

From your RK4 function, your system is:

$$
\begin{align*}
\frac{dx_1}{dt} &= \frac{x_2 - x_1 / R}{C} \\
\frac{dx_2}{dt} &= \frac{V_b(t) - x_1}{L}
\end{align*}
$$

You’re solving this system:
$$
\frac{d\mathbf{x}}{dt} = \mathbf{f}(t, \mathbf{x})
$$
where $ \mathbf{x} = \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} $

---

## 🧪 **Jacobian Matrix Setup**

The **Jacobian** $ J = \frac{\partial \mathbf{f}}{\partial \mathbf{x}} $ is:

$$
J =
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} \\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2}
\end{bmatrix}
=
\begin{bmatrix}
-\frac{1}{RC} & \frac{1}{C} \\
-\frac{1}{L} & 0
\end{bmatrix}
$$

This matrix tells you how sensitive your system is to perturbations in the state variables. It’s especially useful for:
- Estimating **stiffness** (large eigenvalue spread means stiffness)
- Setting up **implicit methods** or **adaptive solvers**

---

## 🧠 **Stiffness Check Tip**

To check stiffness numerically:
- Compute eigenvalues of $ J $ at sample points.
- If there's a **large ratio between smallest and largest eigenvalues**, the system is stiff.
- You can also simulate with small vs. large timesteps and compare results for clues.

---


I'll create a Fortran program that solves the RLC low-pass filter system using the DOPRI5 solver, based on your reference code. First, let me derive the system of ODEs from your circuit equations.

From your equations:
1. $V = i_R R$
2. $i_C = C \frac{dV}{dt}$
3. $L \frac{di_L}{dt} = V_b - V$
4. $i_L = i_R + i_C$

We can derive:
- From (1): $i_R = \frac{V}{R}$
- From (4) and above: $i_C = i_L - \frac{V}{R}$
- From (2): $C \frac{dV}{dt} = i_L - \frac{V}{R}$, so $\frac{dV}{dt} = \frac{1}{C}(i_L - \frac{V}{R})$
- From (3): $\frac{di_L}{dt} = \frac{1}{L}(V_b - V)$

This gives us the system:
- $\frac{dV}{dt} = \frac{1}{C}(i_L - \frac{V}{R})$
- $\frac{di_L}{dt} = \frac{1}{L}(V_b - V)$This Fortran program solves your RLC low-pass filter system using the DOPRI5 solver. Here's what the code does:

## Key Features:

1. **Circuit Analysis**: The program automatically calculates important circuit characteristics:
   - Natural frequency: ω_n = 1/√(LC) ≈ 31.62 rad/s
   - Damping ratio: ζ = 0.5R√(C/L) ≈ 1.58 (overdamped)
   - Time constant: τ = L/R = 0.01 s

2. **System of ODEs**: Converts your circuit equations into:
   - dV/dt = (1/C)(i_L - V/R)
   - di_L/dt = (1/L)(V_b - V)

3. **Initial Conditions**: Starts with V = 0V and i_L = 0A (circuit at rest)

4. **Integration**: Uses DOPRI5 to solve from t = 0 to t = 5τ = 0.05s (enough time to see settling)

5. **Output**: Provides detailed time-series data showing voltage and current evolution

## Expected Behavior:

With your parameters (Vb=24V, L=1H, R=100Ω, C=1mF), the system is **overdamped** (ζ > 1), so you'll see:
- No oscillations
- Smooth exponential approach to steady state
- Final voltage approaches 24V
- Final current approaches 0A

The circuit acts as a low-pass filter, and at DC (t→∞), the capacitor blocks current flow, so the final voltage equals the input voltage.

To compile and run:
```bash
gfortran rlc_dopri5.f90 dopri5.f -o rlc_dopri5
./rlc_dopri5
```

The program will output a table showing time, voltage, and current evolution, plus detailed statistics about the integration and circuit behavior.