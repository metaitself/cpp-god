# Schwarzschild Black Hole Raytracer

A minimal, single-file C++ relativistic raytracer that simulates null geodesics (photon trajectories) and accretion disk imaging around a static Schwarzschild black hole using explicit Fourth-Order Runge-Kutta (RK4) numerical integration.

The project contains zero external dependencies and renders real-time ray-traced frames directly to standard output as an ASCII projection.

---

## Theoretical Background

The spacetime geometry around a non-rotating, spherically symmetric mass $M$ is described by the **Schwarzschild Metric**:

$$ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta d\phi^2)$$

Calculations operate in geometrized units where $G = c = M = 1$.

### Geodesic Equations & Effective Acceleration

Photons follow null geodesics ($ds^2 = 0$). In a 3D Cartesian formulation, the spatial trajectory can be integrated using an effective relativistic acceleration term derived from the conserved angular momentum $\vec{L} = \vec{r} \times \vec{v}$:

$$\vec{a} = - \frac{3 M L^2}{r^5} \vec{r} = - \frac{1.5 \, r_s \, L^2}{r^5} \vec{r}$$

where:
- $r_s = 2M$ is the Schwarzschild radius (event horizon).
- $L^2 = |\vec{r} \times \vec{v}|^2$ is the specific angular momentum magnitude squared.
- $\vec{r}$ and $\vec{v}$ represent the spatial position and coordinate velocity vectors of the photon.

---

## Implementation Details

- **Numerical Integration**: Photon paths are integrated backwards from the observer camera using explicit **RK4** with an adaptive-safe fixed time step (`DT = 0.05`).
- **Event Horizon Capture**: Rays entering $r \le 1.02 \, r_s$ are flagged as captured by the singularity and rendered as shadow regions.
- **Accretion Disk Intersection**: A geometrically thin, optically thick disk is positioned at the equatorial plane ($z = 0$), bounded between the Innermost Stable Circular Orbit ($r_{\text{ISCO}} = 6M$) and outer radius ($r_{\text{max}} = 12M$).
- **Intensity & Redshift Profile**: Disk radiation intensity follows a power-law falloff profile:
  $$I(r) \approx \left(\frac{r_{\text{ISCO}}}{r}\right)^{1.5}$$
- **Camera Model**: Pinhole camera matrix with customizable field of view (FOV), target tracking, and normalized screen-space ray distribution.

---

## Build and Run

A standard C++17 or C++20 compliant compiler is required.

### GCC / Clang (Linux & macOS)

```bash
g++ -O3 -std=c++20 main.cpp -o blackhole
./blackhole
