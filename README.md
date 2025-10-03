# PDELab - Partial Differential Equation Solver

A C++ implementation of numerical methods for solving partial differential equations (PDEs) including heat equations, wave equations, and Laplace equations.

## Features

### Supported PDEs
- **Heat Equation**: Diffusion processes with multiple numerical schemes
- **Wave Equation**: Hyperbolic PDE for wave propagation
- **Laplace Equation**: Elliptic boundary value problems

### Numerical Methods
- **Explicit Schemes**: Forward Time Central Space (FTCS)
- **Implicit Schemes**: Backward Euler with Thomas algorithm
- **Crank-Nicolson**: Second-order accurate in time
- **Jacobi Iteration**: For elliptic problems

## Quick Start

### Prerequisites
- C++ compiler with C++11 support or later
- Standard library support

### Compilation
```bash
g++ -o pde pde.cpp -std=c++11
```

### Running
```bash
./pde
```

The program will generate CSV output files for visualization:
- `heat_explicit.csv` - Heat equation with explicit scheme
- `heat_implicit.csv` - Heat equation with implicit scheme
- `heat_cn.csv` - Heat equation with Crank-Nicolson
- `wave_explicit.csv` - Wave equation solution
- `laplace.csv` - Laplace equation solution

## Mathematical Background

### Heat Equation
Solves: ∂u/∂t = α ∂²u/∂x²

**Schemes implemented:**
- **Explicit (FTCS)**: Stability condition dt ≤ 0.5 dx²/α
- **Implicit**: Unconditionally stable
- **Crank-Nicolson**: Second-order accurate, unconditionally stable

### Wave Equation
Solves: ∂²u/∂t² = c² ∂²u/∂x²

**Scheme implemented:**
- **Explicit**: Second-order central differences with CFL condition c·dt/dx ≤ 1

### Laplace Equation
Solves: ∂²u/∂x² = f(x)

**Scheme implemented:**
- **Jacobi Iteration**: Iterative solver for elliptic boundary value problems

## Code Structure

### Key Components

#### `PDEParams` Structure
```cpp
struct PDEParams {
    double L = 1.0;      // Domain length
    int Nx = 101;        // Number of spatial points
    double T = 0.1;      // Final time
    double alpha = 0.01; // Diffusion coefficient
    double c = 1.0;      // Wave speed
    double dt = 0.0;     // Time step (auto if 0)
    int Nt = 0;          // Number of time steps
};
```

#### Main Functions
- `heat_explicit()` - Explicit heat equation solver
- `heat_implicit()` - Implicit heat equation solver  
- `heat_crank_nicolson()` - Crank-Nicolson heat equation solver
- `wave_explicit()` - Explicit wave equation solver
- `laplace_jacobi()` - Jacobi solver for Laplace equation
- `thomas_solve()` - Tridiagonal matrix solver

## Customization

### Boundary Conditions
Functions accept lambda expressions for boundary conditions:
```cpp
auto zeroBC = [](double t){ return 0.0; };
auto sinBC = [](double t){ return sin(t); };
```

### Initial Conditions
```cpp
auto ic = [](double x){ return sin(M_PI * x); };
auto gaussianIC = [](double x){ return exp(-10*(x-0.5)*(x-0.5)); };
```

### Source Terms (for Poisson equation)
```cpp
auto source = [](double x){ return -2*M_PI*M_PI*sin(M_PI*x); };
```

## Example Usage

```cpp
PDEParams params;
params.L = 1.0;     // Unit domain
params.Nx = 101;    // 101 grid points
params.T = 0.1;     // Solve until t=0.1
params.alpha = 0.01; // Diffusion coefficient

auto ic = [](double x){ return sin(M_PI * x); };
auto bc = [](double t){ return 0.0; };

heat_explicit(params, ic, bc, bc, "output.csv");
```

## Output Format

All solvers output CSV files with columns:
- `x`: Spatial coordinate
- `u`: Solution value at final time

These can be easily imported into plotting software like Python/matplotlib, MATLAB, or Excel.

## Stability Considerations

- **Explicit schemes**: Require stability conditions to prevent numerical instability
- **Implicit schemes**: Unconditionally stable but require solving linear systems
- **CFL conditions**: Automatically enforced when `dt = 0` (auto time step)

## Contributing

Feel free to contribute by:
- Adding new PDE types
- Implementing additional numerical schemes
- Improving documentation
- Adding visualization tools

## License

This project is open source. Feel free to use and modify for educational and research purposes.

## References

- Numerical Methods for Partial Differential Equations
- Finite Difference Methods for PDEs
- Computational Fluid Dynamics literature

---

*Built for educational purposes and numerical experimentation with PDEs.*
