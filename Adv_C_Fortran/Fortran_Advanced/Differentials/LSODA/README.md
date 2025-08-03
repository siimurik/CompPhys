# DLSODE - Livermore Solver for Ordinary Differential Equations

## Overview

DLSODE is a robust, double-precision FORTRAN subroutine for solving systems of first-order ordinary differential equations (ODEs) with initial conditions. Developed at Lawrence Livermore National Laboratory by Alan C. Hindmarsh, DLSODE can handle both stiff and nonstiff systems efficiently.

**Key Features:**
- Solves initial value problems of the form: `dy/dt = f(t,y)` where `y` is a vector
- Handles both stiff and nonstiff systems
- Multiple integration methods (Adams for nonstiff, BDF for stiff)
- Automatic step size and order control
- User-selectable error tolerances
- Optional Jacobian matrix specification

## Mathematical Problem

DLSODE solves systems of the form:

```
dy/dt = f(t,y)
```

or in component form:

```
dy(i)/dt = f(i,t,y(1),y(2),...,y(NEQ))  for i = 1,...,NEQ
```

with initial conditions `y(t₀) = y₀`.

## Subroutine Call

```fortran
CALL DLSODE(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, 
           ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
```

## Parameters

### Input Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `F` | EXTERNAL | User-supplied subroutine defining the ODE system |
| `NEQ` | INTEGER | Number of first-order ODEs |
| `TOUT` | DOUBLE PRECISION | Next output point desired |
| `ITOL` | INTEGER | Tolerance type flag (1 or 2) |
| `RTOL` | DOUBLE PRECISION | Relative tolerance parameter |
| `ATOL` | DOUBLE PRECISION | Absolute tolerance parameter (scalar or array) |
| `ITASK` | INTEGER | Task flag (typically 1 for normal integration) |
| `IOPT` | INTEGER | Optional input flag (0=no, 1=yes) |
| `LRW` | INTEGER | Length of RWORK array |
| `LIW` | INTEGER | Length of IWORK array |
| `JAC` | EXTERNAL | User-supplied Jacobian routine (if needed) |
| `MF` | INTEGER | Method flag specifying solution technique |

### Input/Output Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `Y` | DOUBLE PRECISION array | Solution vector (input: initial values, output: solution at TOUT) |
| `T` | DOUBLE PRECISION | Independent variable (input: initial time, output: final time) |
| `ISTATE` | INTEGER | Integration state flag |

### Work Arrays

| Parameter | Type | Description |
|-----------|------|-------------|
| `RWORK` | DOUBLE PRECISION array | Real work array |
| `IWORK` | INTEGER array | Integer work array |

## Method Flags (MF)

The method flag `MF = 10*METH + MITER` specifies the solution approach:

### METH (Integration Method)
- **1**: Adams method (for nonstiff problems)
- **2**: BDF (Backward Differentiation Formula) method (for stiff problems)

### MITER (Iteration Method)
- **0**: Functional iteration (no Jacobian)
- **1**: Chord iteration with user-supplied full Jacobian
- **2**: Chord iteration with internally generated full Jacobian
- **3**: Chord iteration with internally generated diagonal Jacobian
- **4**: Chord iteration with user-supplied banded Jacobian
- **5**: Chord iteration with internally generated banded Jacobian

### Common MF Values
- **MF = 10**: Nonstiff problems, no Jacobian
- **MF = 21**: Stiff problems, user-supplied full Jacobian
- **MF = 22**: Stiff problems, internally generated full Jacobian
- **MF = 24**: Stiff problems, user-supplied banded Jacobian
- **MF = 25**: Stiff problems, internally generated banded Jacobian

## Required Work Array Sizes

### RWORK Array Length (LRW)
- **MF = 10**: `20 + 16*NEQ`
- **MF = 21 or 22**: `22 + 9*NEQ + NEQ²`
- **MF = 24 or 25**: `22 + 10*NEQ + (2*ML + MU)*NEQ`

### IWORK Array Length (LIW)
- **MF = 10**: `20`
- **MF = 21, 22, 24, or 25**: `20 + NEQ`

## User-Supplied Subroutines

### Function Subroutine F
```fortran
SUBROUTINE F(NEQ, T, Y, YDOT)
INTEGER NEQ
DOUBLE PRECISION T, Y(*), YDOT(*)
```
- **Input**: `NEQ`, `T`, `Y`
- **Output**: `YDOT(i) = f(i,T,Y(1),...,Y(NEQ))`

### Jacobian Subroutine JAC (if MF = 21 or 24)
```fortran
SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
INTEGER NEQ, ML, MU, NROWPD
DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
```
- For full Jacobian (MF = 21): Load `PD(i,j)` with `∂f(i)/∂y(j)`
- For banded Jacobian (MF = 24): Load `PD(i-j+MU+1,j)` with `∂f(i)/∂y(j)`

## Error Control

Error is controlled by the weighted RMS norm:
```
||e||_RMS = sqrt(Σ(e(i)/EWT(i))²/NEQ) ≤ 1
```

Where the error weights are:
```
EWT(i) = RTOL*|Y(i)| + ATOL     (if ITOL = 1)
EWT(i) = RTOL*|Y(i)| + ATOL(i)  (if ITOL = 2)
```

## ISTATE Return Values

### Successful Returns
- **1**: No integration performed (TOUT = T)
- **2**: Successful integration

### Error Returns
- **-1**: Excessive work done (increase MXSTEP)
- **-2**: Too much accuracy requested (relax tolerances)
- **-3**: Illegal input detected
- **-4**: Repeated error test failures
- **-5**: Repeated convergence failures (check Jacobian)
- **-6**: Error weight became zero

## Basic Usage Example

```fortran
PROGRAM EXAMPLE
    EXTERNAL FEX, JEX
    INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
    INTEGER IWORK(23)
    DOUBLE PRECISION Y(3), T, TOUT, RTOL, ATOL(3), RWORK(58)
    
    ! Problem setup
    NEQ = 3
    Y(1) = 1.0D0
    Y(2) = 0.0D0  
    Y(3) = 0.0D0
    T = 0.0D0
    TOUT = 0.4D0
    
    ! Tolerance setup
    ITOL = 2
    RTOL = 1.0D-4
    ATOL(1) = 1.0D-6
    ATOL(2) = 1.0D-10
    ATOL(3) = 1.0D-6
    
    ! Method setup
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    LRW = 58
    LIW = 23
    MF = 21  ! Stiff method with user-supplied Jacobian
    
    ! Integration call
    CALL DLSODE(FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
   *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
    
    ! Check results
    IF (ISTATE .LT. 0) THEN
        WRITE(*,*) 'Error: ISTATE =', ISTATE
    ELSE
        WRITE(*,*) 'Solution at T =', T
        WRITE(*,*) 'Y =', Y(1), Y(2), Y(3)
    ENDIF
END

SUBROUTINE FEX(NEQ, T, Y, YDOT)
    INTEGER NEQ
    DOUBLE PRECISION T, Y(3), YDOT(3)
    YDOT(1) = -0.04D0*Y(1) + 1.0D4*Y(2)*Y(3)
    YDOT(3) = 3.0D7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
    RETURN
END

SUBROUTINE JEX(NEQ, T, Y, ML, MU, PD, NRPD)
    INTEGER NEQ, ML, MU, NRPD
    DOUBLE PRECISION T, Y(3), PD(NRPD,3)
    PD(1,1) = -0.04D0
    PD(1,2) = 1.0D4*Y(3)
    PD(1,3) = 1.0D4*Y(2)
    PD(2,1) = 0.04D0
    PD(2,3) = -PD(1,3)
    PD(3,2) = 6.0D7*Y(2)
    PD(2,2) = -PD(1,2) - PD(3,2)
    RETURN
END
```

## Optional Inputs and Outputs

### Optional Inputs (when IOPT = 1)
- `RWORK(5)`: Initial step size (H0)
- `RWORK(6)`: Maximum step size (HMAX)
- `RWORK(7)`: Minimum step size (HMIN)
- `IWORK(5)`: Maximum order (MAXORD)
- `IWORK(6)`: Maximum steps (MXSTEP, default 500)
- `IWORK(7)`: Maximum "H = 0" warnings (MXHNIL, default 10)

### Optional Outputs
- `RWORK(11)`: Last step size used (HU)
- `RWORK(12)`: Next step size (HCUR)
- `RWORK(13)`: Current T value reached (TCUR)
- `IWORK(11)`: Number of steps taken (NST)
- `IWORK(12)`: Number of function evaluations (NFE)
- `IWORK(13)`: Number of Jacobian evaluations (NJE)
- `IWORK(14)`: Last order used (NQU)

## Performance Tips

1. **Choose appropriate MF**: Use MF=10 for nonstiff problems, MF=22 for stiff problems when Jacobian computation is difficult
2. **Provide Jacobian when possible**: User-supplied Jacobians (MF=21,24) are more efficient than numerical approximations
3. **Set reasonable tolerances**: Too tight tolerances waste computation; too loose tolerances give inaccurate results
4. **Use banded Jacobian for sparse systems**: MF=24,25 for systems with banded Jacobian structure

## Common Issues and Solutions

1. **ISTATE = -1**: Increase `IWORK(6)` (MXSTEP) if legitimate
2. **ISTATE = -2**: Relax tolerances (increase RTOL, ATOL)
3. **ISTATE = -4,-5**: Check problem setup, Jacobian correctness, or try different MF
4. **Poor performance**: Ensure correct stiff/nonstiff classification and appropriate MF choice

## References

- Hindmarsh, A.C., "ODEPACK, A Systematized Collection of ODE Solvers," in Scientific Computing, R.S. Stepleman et al. (Eds.), North-Holland, Amsterdam, 1983, pp. 55-64.
- Based on GEAR and GEARB packages
- Part of the ODEPACK collection

## Notes

- DLSODE is not re-entrant
- Work arrays should not be altered between calls for the same problem
- For maximum efficiency on Cray systems, use CFT77 compiler
- The solver automatically handles step size and order selection
