# N-Body Gravitational Simulation

A C++17 implementation of an N-body gravitational simulation using the Velocity Verlet integration scheme. This project demonstrates numerical methods for solving classical mechanics problems, specifically orbital dynamics of the Sun-Earth-Moon system.

## Features

- **Velocity Verlet Integration**: Symplectic integrator that conserves energy over long timescales
- **Modular Design**: Clean separation between physics, I/O, and data structures
- **Realistic Parameters**: Uses NASA JPL data for celestial body masses and orbital parameters
- **Visualization**: Python script generates animated orbital trajectories

## Quick Start

```bash
# Build the simulation
make

# Run simulation (outputs to simulation_output.csv)
make run

# Generate visualization (requires Python + matplotlib)
make plot

# Or do everything at once
make simulate
```

## Requirements

**C++ Compilation:**
- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)

**Visualization (optional):**
```bash
pip install matplotlib pandas pillow
```

## Physics Background

### Newton's Law of Universal Gravitation

Every particle in the universe attracts every other particle with a force proportional to the product of their masses and inversely proportional to the square of the distance between them:

$$\vec{F}_{ij} = G \frac{m_i m_j}{|\vec{r}_{ij}|^2} \hat{r}_{ij}$$

where:
- $G = 6.67430 \times 10^{-11} \, \text{m}^3 \text{kg}^{-1} \text{s}^{-2}$ is the gravitational constant
- $\vec{r}_{ij}$ is the displacement vector from body $i$ to body $j$

### The N-Body Problem

For $N$ interacting bodies, each body experiences gravitational acceleration from all others:

$$\vec{a}_i = \sum_{j \neq i} G \frac{m_j}{|\vec{r}_{ij}|^3} \vec{r}_{ij}$$

This system of coupled differential equations has no general closed-form solution for $N > 2$, necessitating numerical integration.

### Velocity Verlet Algorithm

The Velocity Verlet method is a second-order symplectic integrator. Unlike simple Euler integration, it:

1. **Conserves energy** (crucial for orbital mechanics)
2. **Is time-reversible** (a fundamental property of Hamiltonian systems)
3. **Has O(dt²) local error** vs O(dt) for Euler

The algorithm:

```
1. x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²
2. Compute a(t+dt) from new positions
3. v(t+dt) = v(t) + ½·[a(t) + a(t+dt)]·dt
```

### Why Verlet Over Runge-Kutta?

While RK4 has higher local accuracy (O(dt⁴)), Verlet's symplectic nature makes it superior for long-term orbital integration:

| Method | Local Error | Energy Drift | Time-Reversible |
|--------|-------------|--------------|-----------------|
| Euler | O(dt) | Unbounded | No |
| RK4 | O(dt⁴) | Linear | No |
| Verlet | O(dt²) | Bounded | Yes |

For orbital simulations spanning years, Verlet maintains bounded energy error while RK4 accumulates secular drift.

## Project Structure

```
.
├── main.cpp              # Simulation implementation
├── Makefile              # Build configuration
├── plot.py               # Visualization script
├── README.md             # This file
├── simulation_output.csv # Generated trajectory data
└── orbit_animation.gif   # Generated visualization
```

## Code Architecture

### Core Data Structures

```cpp
struct Vec3 {
    double x, y, z;
    // Operator overloading for physics-like syntax
};

struct Body {
    std::string name;
    double mass;
    Vec3 position, velocity, acceleration;
};
```

### Key Functions

| Function | Purpose |
|----------|---------|
| `compute_acceleration()` | Calculates gravitational force summation |
| `velocity_verlet_step()` | Advances system by one timestep |
| `create_sun_earth_moon_system()` | Initializes realistic orbital parameters |

## Simulation Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Time step | 3600 s (1 hour) | Resolves lunar orbit (~27 days) |
| Duration | 1 year | Shows complete Earth orbit |
| Output interval | 24 hours | Balances resolution vs file size |

## Extending the Simulation

### Adding Bodies

```cpp
// In create_sun_earth_moon_system():
bodies.emplace_back("Mars", 6.417e23,
    Vec3(2.279e11, 0, 0),  // Semi-major axis
    Vec3(0, 24077, 0));    // Orbital velocity
```

### Higher-Order Methods

For improved accuracy, consider implementing:
- **Forest-Ruth** (4th order symplectic)
- **Yoshida** (4th/6th order symplectic)
- **Gauss-Legendre** (implicit, A-stable)

## Performance

The simulation has O(N²) complexity per timestep. For large N, consider:
- Barnes-Hut tree algorithm: O(N log N)
- Fast Multipole Method: O(N)

## Output Format

`simulation_output.csv` contains:

```csv
time,Sun_x,Sun_y,Sun_z,Earth_x,Earth_y,Earth_z,Moon_x,Moon_y,Moon_z
0.000000e+00,0.000000e+00,0.000000e+00,...
```

All positions in SI units (meters).

## References

1. Verlet, L. (1967). "Computer Experiments on Classical Fluids". *Physical Review*, 159(1), 98-103.
2. Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration*. Springer.
3. NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/

## License

MIT License - See LICENSE file for details.

## Author

Makar Ulesov
