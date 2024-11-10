# Functions Documentation

## 1. Lane-Emden System Functions

### `lane_emden_sys(xi, S, n)`
Defines the Lane-Emden system of differential equations.

- **Parameters:**
  - `xi (float)`: The dimensionless radial coordinate.
  - `S (list)`: A list containing `y1 = θ(ξ)` and `y2 = dθ/dξ`.
  - `n (float)`: The polytropic index.
- **Returns:** A list `[dy1/dxi, dy2/dxi]` representing derivatives.

### `lane_emden_solver(n)`
Solves the Lane-Emden system for a given polytropic index `n`.

- **Parameters:**
  - `n (float)`: The polytropic index.
- **Returns:** A tuple `(x, y1, y2)` where:
  - `x`: Radial coordinate array.
  - `y1`: Values of `θ(ξ)`.
  - `y2`: Values of `dθ/dξ`.

## 2. Chemical Composition Function

### `chemical_composition(xi, X_s, r_c, d_c, t, t_tot)`
Calculates hydrogen abundance at a given radius.

- **Parameters:**
  - `xi (array)`: Dimensionless radial coordinates.
  - `X_s (float)`: Initial hydrogen abundance.
  - `r_c (float)`: Radius-related parameter.
  - `d_c (float)`: Diffusion-related parameter.
  - `t (float)`: Current time.
  - `t_tot (float)`: Total time.
- **Returns:** Array `X` representing hydrogen abundance at each radial point.

## 3. Density, Pressure, and Temperature Calculation

### `rho_P_T_central(n, R, M, X, Z)`
Calculates the central density, pressure, and temperature in a polytropic star model.

- **Parameters:**
  - `n (float)`: Polytropic index.
  - `R (float)`: Radius of the star.
  - `M (float)`: Mass of the star.
  - `X (array)`: Hydrogen abundance.
  - `Z (float)`: Metallicity.
- **Returns:** Tuple `(rho_c, P_c, T_c)` where:
  - `rho_c`: Central density.
  - `P_c`: Central pressure.
  - `T_c`: Central temperature.

## 4. Density, Pressure, and Temperature Calculation

### `rho_P_T(n, R, M, X, Z)`
Calculates the density, pressure, and temperature at each radius in a polytropic star.

- **Parameters:**
  - `n (float)`: Polytropic index.
  - `R (float)`: Radius of the star.
  - `M (float)`: Mass of the star.
  - `X (array)`: Hydrogen abundance.
  - `Z (float)`: Metallicity.
- **Returns:** Tuple `(rho, P, T)` for density, pressure, and temperature arrays.

## 5. Emissivity Calculation

### `emissivity(rho, T, X, Z)`
Calculates total emissivity from proton-proton and CNO cycle reactions.

- **Parameters:**
  - `rho (array)`: Density array.
  - `T (array)`: Temperature array.
  - `X (array)`: Hydrogen abundance array.
  - `Z (array)`: Metallicity array.
- **Returns:** Tuple `(e_total, e_pp, e_cno)` for total, proton-proton, and CNO emissivity.

## 6. Luminosity Calculation

### `luminosity(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol)`
Calculates star luminosity using the Lane-Emden solution, chemical composition, and emissivity.

- **Parameters:**
  - `n (float)`: Polytropic index.
  - `X_s (float)`: Initial hydrogen abundance.
  - `r_c (float)`: Radius-related parameter.
  - `d_c (float)`: Diffusion-related parameter.
  - `R (float)`: Radius of the star.
  - `time (float)`: Current time.
  - `total_time (float)`: Total time.
  - `M (float)`: Mass of the star.
  - `Z (float)`: Metallicity.
  - `L_sol (float)`: Solar luminosity constant.
- **Returns:** Float representing the calculated luminosity.

## 7. Luminosity Comparison Function

### `y(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol, L)`
Calculates a logarithmic comparison between the calculated luminosity and a given luminosity.

- **Parameters:** Same as `luminosity`, with an additional:
  - `L (float)`: Given luminosity for comparison.
- **Returns:** Float for `log10(lumin / (L / L_sol))`.

## 8. Polytropic Index Finder

### `index(f, a, b, M, L, Y, Z_X, R, L_sol, r_c, d_c, time, total_time)`
Finds the polytropic index using the bisection method.

- **Parameters:**
  - `f (function)`: Function for the root-finding problem.
  - `a (float)`, `b (float)`: Bounds for the bisection method.
  - Other parameters: Stellar properties and conditions.
- **Returns:** Polytropic index.

## 9. Monte Carlo Simulation

### `monte_carlo(M, dM, R, dR, L, dL, Z_X, dZ_X, Y, dY, r_n, d_c, tau, tau_t, N, L_sol, a, b)`
Performs Monte Carlo simulations to calculate polytropic index distributions.

- **Parameters:**
  - `M`, `dM`, `R`, `dR`, `L`, `dL`: Mean and uncertainty for stellar properties.
  - `Z_X`, `dZ_X`, `Y`, `dY`: Metallicity to hydrogen ratio and helium mass fraction with uncertainties.
  - `r_n`, `d_c`, `tau`, `tau_t`: Star parameters and timescales.
  - `N (int)`: Number of iterations.
  - `L_sol (float)`: Solar luminosity constant.
  - `a`, `b`: Bounds for the bisection method.
- **Returns:** Array of calculated polytropic indices.

## 10. Visualization Functions

### `lane_emden_plots()`
Generates plots of Lane-Emden solutions for various polytropic indices.

- **Parameters:** None
- **Returns:** None (Displays plots)

### `monte_carlo_plots(mean, std, N, nn)`
Generates visualizations of Monte Carlo simulation results and compares with a Gaussian distribution.

- **Parameters:**
  - `mean (float)`: Mean of the Monte Carlo results.
  - `std (float)`: Standard deviation.
  - `N (int)`: Number of points for Gaussian visualization.
  - `nn (np.ndarray)`: Simulated polytropic index values.
- **Returns:** None (Displays plots)

### `mass_lum_plots(mass_e, mass_t, mass_sun, lum_sun, lum_e, lum_t)`
Generates plots comparing the luminosity as a function of mass for the Sun, Epsilon Eridani, and Theta Persei A.

- **Parameters:**
  - Arrays for mass distributions and corresponding luminosities of the three stars.
- **Returns:** None (Displays plots)


### `Examples`

```

#Sol
M_sol = 1.988475e33
L_sol = 3.828e33
R_sol = 6.957e10
Z_sol = 0.02857
Y_sol = 0.28


# Epsilon Eridani
M_e = 0.82 * M_sol
dM_e = 0.02 * M_sol
R_e = 0.738 * R_sol
dR_e= 0.0003 * R_sol
L_e = 0.32 * L_sol
dL_e = 0.01 * L_sol
metal_e = -0.08
dmetal_e = 0.01
Z_e = Z_sol*10**metal_e
Y_e = 0.2423
dY_e = 0.0054
total_time_e = 2*4.6*M_e*(1/L_e)


# Theta Persei
M_t = 1.138 * M_sol
dM_t = 0.010 * M_sol
R_t = 1.319 * R_sol
dR_t = 0.011 * R_sol
L_t = 2.235 * L_sol
dL_t = 0.040 * L_sol
metal_t = -0.03
dmetal_t = 0.09
Z_t = Z_e = Z_sol*10**metal_t
Y_t = 0.2423
dY_t = 0.0054
total_time_t = 2*4.6*M_e*(1/L_t)


# Monte Carlo plots (Theta Persei A)
nn_t = monte_carlo(M_t, dM_t, R_t, dR_t, L_t, dL_t,  Z_t, dmetal_t, Y_t, dY_t, 0.2, 0.03, total_time_t/2, total_time_t, 3000, L_sol, a=1, b=3.9)
std_t = np.std(nn_t)
mean_t = np.mean(nn_t)
monte_carlo_plots(mean_t, std_t, N=3000, nn_t)


# Indexes
poly_index_t = index(y, 2.5, 4.5, M_t, L_t, Y_t, Z_t, R_t, L_sol, 0.2, 0.03, total_time_t/3, total_time_t)
poly_index_e = index(y, 1.2, 4.5, M_e, L_e, Y_e, Z_e, R_e, L_sol, 0.2, 0.03, total_time_e/2, total_time_e)
poly_index_sun = index(y, 1, 4, M_sol, L_sol, Y_sol, Z_sol, R_sol,  L_sol, 0.2, 0.03, 1, 2)

# Theta
Xs_t = (1-Y_t)/(1+Z_t)
Z = 1-Y_t-Xs_t
x_t, y, z = lane_emden_solver(poly_index_t)
X_t = chemical_composition(x_t, Xs_t, 0.2, 0.03, total_time_t/2, total_time_t)

# Epsilon
Xs_e = (1-Y_e)/(1+Z_e)
Z = 1-Y_e-Xs_e
x_e, y, z = lane_emden_solver(poly_index_e)
X_e = chemical_composition(x_e, Xs_e, 0.2, 0.03, total_time_e/2, total_time_e)

# Sun
Xs_sun = (1-Y_sol)/(1+Z_sol)
Z = 1-Y_sol-Xs_sun
x_sol, y, z = lane_emden_solver(poly_index_sun)
X_sun = chemical_composition(x_sol, Xs_sun, 0.2, 0.03, 1, 2)

# Thermodynamic
rho_sun, P_sun, T_sun = rho_P_T(poly_index_sun, R_sol, M_sol, X_sun, Z_sol)
rho_e, P_e, T_e = rho_P_T(poly_index_e, R_e, M_e, X_e, Z_e)
rho_t, P_t, T_t = rho_P_T(poly_index_t, R_t, M_t, X_t, Z_t)


rho_sun_c, P_sun_c, T_sun_c = rho_P_T_central(poly_index_sun, R_sol, M_sol, X_sun, Z_sol)
rho_e_c, P_e_c, T_e_c = rho_P_T_central(poly_index_e, R_e, M_e, X_e, Z_e)
rho_t_c, P_t_c, T_t_c = rho_P_T_central(poly_index_t, R_t, M_t, X_t, Z_t)

# Luminosity
lum_sun = luminosity(poly_index_sun, Xs_sun, 0.2, 0.03, R_sol, 1, 2,  M_sol, Z_sol, L_sol )
lum_e = luminosity(poly_index_e, Xs_e, 0.2, 0.03, R_e, total_time_e/3, total_time_e, M_e, Z_e, L_sol )
lum_t = luminosity(poly_index_t, Xs_t, 0.2, 0.03, R_t, total_time_t/2, total_time_t, M_t, Z_t, L_sol )

# Mass
mass_e = mass(poly_index_e, M_e)
mass_t = mass(poly_index_t, M_t)
mass_sun = mass(poly_index_sun, M_sol)


# Emissivity
e_t, epp_t, ecno_t = emissivity(rho_t, T_t, X_t, Z_t)
e_e, epp_e, ecno_e = emissivity(rho_e, T_e, X_e, Z_e)
e_sun, epp_sun, ecno_sun = emissivity(rho_sun, T_sun, X_sun, Z_sol)


# Plots
thermo_plots(x_sol, x_e, x_t, rho_sun, rho_e, rho_t, P_sun, P_e, P_t, T_sun, T_e, T_t)
emissivity_plots(x_sol, x_e, x_t, e_t, epp_t, ecno_t, e_e, epp_e, ecno_e, e_sun, epp_sun, ecno_sun)
mass_lum_plots(mass_e, mass_t, mass_sun, lum_sun, lum_e, lum_t)


```