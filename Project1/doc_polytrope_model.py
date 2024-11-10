# Imports
import numpy as np
from scipy.integrate import odeint, solve_ivp, cumulative_simpson, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import bisect
import matplotlib.pyplot as plt
from scipy.stats import norm

# Global vars
G = 6.67408e-8
R_g = 8.314511e7

# Lane-Emden System
def lane_emden_sys(xi, S, n):
    """
    Defines the Lane-Emden system of differential equations.
    
    Parameters:
    xi : float
        The dimensionless radial coordinate.
    S : list
        A list containing two elements: y1 and y2, where y1 = θ(ξ) and y2 = dθ/dξ.
    n : float
        The polytropic index.
    
    Returns:
    list
        A list of two elements representing the derivatives [dy1/dxi, dy2/dxi].
    """
    y1, y2 = S
    dy1_dxi = y2
    dy2_dxi = - np.abs(y1)**n - 2*y2/xi 
    return [dy1_dxi, dy2_dxi]

def lane_emden_solver(n):
    """
    Solves the Lane-Emden system for a given polytropic index n.
    
    Parameters:
    n : float
        The polytropic index.
    
    Returns:
    tuple
        A tuple (x, y1, y2) where x is the radial coordinate, y1 is θ(ξ), and y2 is dθ/dξ.
    """
    # xi Range
    xi_0 = 1e-7  # Prevent numerical singularity
    xi_range = np.linspace(xi_0, 40, 100000)

    # Initial Conditions
    y1_0 = 1 - xi_0**2 / 2  # θ(ξ₀) ≈ 1 - ξ₀²/2
    y2_0 = -xi_0           # dθ/dξ(ξ₀) ≈ -ξ₀
    lane_em_0 = [y1_0, y2_0]

    # Solve Lane-Emden system
    sol = odeint(lane_emden_sys, lane_em_0, xi_range, args=(n,), tfirst=True)
    
    x = np.append([0], xi_range)
    y1 = np.append([1], sol[:, 0])
    y2 = np.append([0], sol[:, 1])

    zero_index = 1 + int(np.argwhere(y1 < 0)[0])

    x = x[:zero_index]
    y1 = y1[:zero_index]
    y2 = y2[:zero_index]

    # Interpolation to find the point where y1 = 0
    interp_func_y1 = interp1d(y1[-2:], x[-2:], kind='linear')
    root = interp_func_y1(0)  # The value of x where y1 crosses zero

    # Update x and y1 with the point where y1 = 0
    x[-1] = root
    y1[-1] = 0

    # Interpolation for y2 at the new point
    interp_func_y2 = interp1d(x[-2:], y2[-2:], kind='linear')
    y2[-1] = interp_func_y2(root)  # Get the interpolated value of y2 at the new point

    return x, y1, y2

# Chemical Composition (Hydrogen and Helium)
def chemical_composition(xi, X_s, r_c, d_c, t, t_tot):
    """
    Calculates the chemical composition (Hydrogen content) at a given radius.

    Parameters:
    xi : array
        The dimensionless radial coordinates.
    X_s : float
        The initial hydrogen abundance.
    r_c : float
        A radius-related parameter.
    d_c : float
        A diffusion-related parameter.
    t : float
        The current time.
    t_tot : float
        The total time.

    Returns:
    X : array
        The hydrogen abundance at each radial point.
    """
    X_c = X_s - X_s * t / t_tot
    X = X_c + (X_s - X_c) / (1 + np.exp((r_c - xi) / xi[-1] / d_c))
    return X

def rho_P_T_central(n, R, M, X, Z):
    """
    Calculates the central density, pressure, and temperature at each radius in a polytropic star model.

    Parameters:
    n : float
        The polytropic index.
    R : float
        The radius of the star.
    M : float
        The mass of the star.
    X : array
        The hydrogen abundance at each radial point.
    Z : float
        The metallicity of the star.

    Returns:
    tuple
        A tuple (rho_c, P_c, T_c) where rho_c is the central density, P_c is the central pressure, and T_c is the central temperature.
    """
    x, y1, y2 = lane_emden_solver(n)
    mu = 4/(5*X-Z+3)
    rho_c = -(M*x[-1])/(4*np.pi*y2[-1]*R**3)
    T_c = -(mu*G*M)/(R_g*R*(n+1)*x[-1]*y2[-1])
    P_c = (G*M**2)/(R**4 * 4*np.pi*(n+1)*y2[-1]**2)
    return rho_c, P_c, T_c

# Thermodynamic functions (density, pressure, temperature)
def rho_P_T(n, R, M, X, Z):
    """
    Calculates the density, pressure, and temperature at each radius in a polytropic star model.

    Parameters:
    n : float
        The polytropic index.
    R : float
        The radius of the star.
    M : float
        The mass of the star.
    X : array
        The hydrogen abundance at each radial point.
    Z : float
        The metallicity of the star.

    Returns:
    tuple
        A tuple (rho, P, T) where rho is the density, P is the pressure, and T is the temperature.
    """
    x, y1, y2 = lane_emden_solver(n)
    mu = 4 / (5 * X - Z + 3)
    rho_c = -(M * x[-1]) / (4 * np.pi * y2[-1] * R**3)
    T_c = -(mu * G * M) / (R_g * R * (n + 1) * x[-1] * y2[-1])
    P_c = (G * M**2) / (R**4 * 4 * np.pi * (n + 1) * y2[-1]**2)
    rho = rho_c * y1**n
    T = T_c * y1
    P = P_c * y1**(n + 1)
    return rho, P, T

def emissivity(rho, T, X, Z):
    """
    Calculates the total emissivity from proton-proton and carbon-nitrogen-oxygen reactions.

    Parameters:
    rho : array
        The density at each radial point.
    T : array
        The temperature at each radial point.
    X : array
        The hydrogen abundance at each radial point.
    Z : array
        The metallicity at each radial point.

    Returns:
    tuple
        A tuple (e_total, e_pp, e_cno) where e_total is the total emissivity, e_pp is the proton-proton emissivity, and e_cno is the carbon-nitrogen-oxygen emissivity.
    """
    T6 = T / 1e6

    T_og_length = len(T6)  # Original length of T6 array
    T_new_length = 1 + int(np.argwhere(T6 < 1e-3)[0])  # Avoid numerically unstable values
    size_dif = T_og_length - T_new_length

    # Update array sizes for emissivity calculations
    T6 = T6[:T_new_length]
    X = X[:T_new_length]
    rho = rho[:T_new_length]

    # Formula variables
    alpha = 1.2e17 * ((1 - X - Z) / (4 * X))**2 * np.exp(-100 * T6**(-1/3))
    eps_0 = 2.38e6 * X**2 * rho * T6**(-2/3) * (1 + 0.0123 * T6**(1/3) + 0.0109 * T6**(2/3) + 0.00095 * T6) * np.exp(-33.80 * T6**(-1/3) + 0.27 * rho**(1/2) * T6**(-3/2))
    phi_alpha = 1 - alpha + np.sqrt(alpha * (alpha + 2))
    F1 = (np.sqrt(alpha + 2) - np.sqrt(alpha)) / (np.sqrt(alpha + 2) + 3 * np.sqrt(alpha))
    F2 = (1 - F1) / (1 + 8.94e15 * (X / (4 - 3 * X)) * T6**(-1/6) * np.exp(-102.65 * T6**(-1/3)))
    F3 = 1 - F1 - F2

    # Emissivity from proton-proton (e_pp) and carbon-nitrogen-oxygen (e_cno) reactions
    e_pp = (eps_0 / 0.980) * phi_alpha * (0.980 * F1 + 0.960 * F2 + 0.721 * F3)
    e_cno = 8.67e27 * Z * X * rho * T6**(-2/3) * (1 + 0.0027 * T6**(1/3) - 0.00778 * T6**(2/3) - 0.000149 * T6) * np.exp(-152.28 * T6**(-1/3))

    # Fill the arrays back to their original sizes with zeros
    zeros = np.zeros(size_dif)
    e_pp = np.concatenate((e_pp, zeros))
    e_cno = np.concatenate((e_cno, zeros))

    # Total emissivity
    e_total = e_pp + e_cno

    return e_total, e_pp, e_cno


def luminosity(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol):
    """
    Calculates the luminosity of a star based on the Lane-Emden equation, chemical composition,
    thermodynamic functions, and emissivity.

    Parameters:
    n : float
        The polytropic index.
    X_s : float
        The initial hydrogen abundance.
    r_c : float
        A radius-related parameter.
    d_c : float
        A diffusion-related parameter.
    R : float
        The radius of the star.
    time : float
        The current time.
    total_time : float
        The total time for the process.
    M : float
        The mass of the star.
    Z : float
        The metallicity of the star.
    L_sol : float
        The solar luminosity constant.

    Returns:
    float
        The luminosity of the star.
    """
    # xi, theta(xi), theta'(xi)
    xi, y1, y2 = lane_emden_solver(n)

    # Chemical composition X, Y
    X = chemical_composition(xi, X_s, r_c, d_c, time, total_time)

    # Thermodynamic functions
    rho, P, T = rho_P_T(n, R, M, X, Z)

    # Emissivity
    ems_tot, e_pp, e_cno = emissivity(rho, T, X, Z)

    # Constant of the luminosity formula
    const = -1/(xi[-1]**2 * y2[-1])

    # Integrand of the luminosity formula
    integrand = xi**2 * y1**n * M * ems_tot / L_sol

    # Cumulative integral
    integral = cumulative_simpson(integrand, x=xi, initial=0)
    return const * integral


def y(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol, L):
    """
    Calculates a logarithmic quantity that compares the calculated luminosity with a given luminosity.

    Parameters:
    n : float
        The polytropic index.
    X_s : float
        The initial hydrogen abundance.
    r_c : float
        A radius-related parameter.
    d_c : float
        A diffusion-related parameter.
    R : float
        The radius of the star.
    time : float
        The current time.
    total_time : float
        The total time for the process.
    M : float
        The mass of the star.
    Z : float
        The metallicity of the star.
    L_sol : float
        The solar luminosity constant.
    L : float
        The luminosity to compare with.

    Returns:
    float
        The logarithmic comparison of the luminosity: log10(lumin / L).
    """
    lumin = luminosity(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol)[-1]
    y_value = np.log10(lumin / (L / L_sol))
    return y_value


def index(f, a, b, M, L, Y, Z_X, R, L_sol, r_c, d_c, time, total_time):
    """
    Finds the value of the polytropic index by solving an equation using the bisection method.

    Parameters:
    f : function
        The function to find the root of (polytropic index equation).
    a : float
        The lower bound for the bisection method.
    b : float
        The upper bound for the bisection method.
    M : float
        The mass of the star.
    L : float
        The luminosity of the star.
    Y : float
        The helium mass fraction.
    Z_X : float
        The ratio of metallicity to hydrogen.
    R : float
        The radius of the star.
    L_sol : float
        The solar luminosity constant.
    r_c : float
        A radius-related parameter.
    d_c : float
        A diffusion-related parameter.
    time : float
        The current time.
    total_time : float
        The total time for the process.

    Returns:
    float
        The value of the polytropic index found using the bisection method.
    """
    X_s = (1 - Y) / (1 + Z_X)
    Z = 1 - Y - X_s
    index_value = bisect(f, a, b, args=(X_s, r_c, d_c, R, time, total_time, M, Z, L_sol, L))
    return index_value


def monte_carlo(M, dM, R, dR, L, dL, Z_X, dZ_X, Y, dY, r_n, d_c, tau, tau_t, N, L_sol, a, b):
    """
    Performs a Monte Carlo simulation to calculate the distribution of polytropic indices based on 
    input uncertainties in various stellar parameters.

    Parameters:
    M : float
        The mass of the star.
    dM : float
        The uncertainty in the mass of the star.
    R : float
        The radius of the star.
    dR : float
        The uncertainty in the radius of the star.
    L : float
        The luminosity of the star.
    dL : float
        The uncertainty in the luminosity of the star.
    Z_X : float
        The ratio of metallicity to hydrogen.
    dZ_X : float
        The uncertainty in Z_X.
    Y : float
        The helium mass fraction.
    dY : float
        The uncertainty in Y.
    r_n : float
        A radius-related parameter.
    d_c : float
        A diffusion-related parameter.
    tau : float
        The current time.
    tau_t : float
        The total time.
    N : int
        The number of Monte Carlo iterations.
    L_sol : float
        The solar luminosity constant.
    a : float
        The lower bound for the bisection method.
    b : float
        The upper bound for the bisection method.

    Returns:
    np.ndarray
        An array of polytropic indices calculated from the Monte Carlo simulation.
    """
    nn = np.zeros(N)
    Mm = M + dM * np.random.normal(0, 1, N)
    Rr = R + dR * np.random.normal(0, 1, N)
    Ll = L + dL * np.random.normal(0, 1, N)
    Z_Xx = Z_X + dZ_X * np.random.normal(0, 1, N)
    Yy = Y + dY * np.random.normal(0, 1, N)
    for i in range(N):
        print(f"Iteration {i + 1}/{N}")
        try:
            # Attempt to calculate the polytropic index with the current parameters
            nn[i] = index(y, a, b, Mm[i], Ll[i], Yy[i], Z_Xx[i], Rr[i], r_n, d_c, tau, tau_t, L_sol)
        except ValueError as e:
            # If there's a ValueError, skip this iteration and print a message
            print(f"Warning: Skipping iteration {i + 1} due to error: {e}")
            # Optionally, set this entry to NaN or some other marker
            nn[i] = np.nan
    return nn

def lane_emden_plots():
    """
    Generates visualizations of the Lane-Emden solutions for different polytropic indices (n) 
    and their corresponding profiles.

    This function creates two plots:
    1. A plot showing the Lane-Emden solutions (theta) as a function of xi for different values of n.
    2. A plot of the scaled xi (xi/xi_s) showing the behavior of theta for different values of n.

    Parameters:
    None

    Returns:
    None
    """
    # List of n values
    ns = [0, 1, 1.5, 3, 4, 4.5]

    # Lane-Emden solution for different values of n
    for n in ns:
        xi_range, y1, y2 = lane_emden_solver(n)
        plt.plot(xi_range, y1, label=f"n = {n}")
        plt.ylim(0, 1)
        plt.xlim(0, 15)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$\theta(\xi)$')
        plt.title('Lane-Emden Solutions for different values of n')
        plt.legend()
        plt.grid(True)
    plt.show()

    # xi_s for different values of n
    for n in ns:
        xi_range, y1, y2 = lane_emden_solver(n)
        xs = xi_range / xi_range[-1]
        plt.plot(xs, y1, label=f"n = {n}")
        plt.ylim(0, 1)
        plt.xlim(0, 1)
        plt.xlabel(r'$\xi$/$\xi_S$')
        plt.ylabel(r'$\theta$')
        plt.title(r'$\xi_S$ for different values of n')
        plt.legend()
        plt.grid(True)
    plt.show()


def monte_carlo_plots(mean, std, N, nn):
    """
    Generates visualizations of the results from a Monte Carlo simulation, showing the distribution 
    of polytropic indices (n) and comparing it with a Gaussian distribution.

    This function creates a plot with:
    1. A histogram of the Monte Carlo simulation results.
    2. A Gaussian distribution curve based on the provided mean and standard deviation.
    3. A vertical line indicating the mean value.

    Parameters:
    mean : float
        The mean value of the Monte Carlo simulation (expected value of the polytropic index).
    std : float
        The standard deviation of the Monte Carlo simulation.
    N : int
        The number of points used to generate the Gaussian distribution for visualization.
    nn : np.ndarray
        The array of simulated polytropic index values from the Monte Carlo simulation.

    Returns:
    None
    """
    # Gaussian
    # Generate x values for the curve
    x_gauss = np.linspace(mean - 5 * std, mean + 5 * std, N)
    # Calculate the PDF of the normal distribution
    y_gauss = norm.pdf(x_gauss, mean, std)

    # Plots
    # Mean vertical line
    plt.axvline(x=mean, color='red', linewidth=2, label='Mean')

    # Gaussian
    plt.plot(x_gauss, y_gauss, color='orange', linewidth=2, label='Gaussian')

    # Histogram
    plt.hist(nn, 40, color='blue', edgecolor='black', density=True)

    plt.xlabel('n')
    plt.ylabel('N(n)')
    plt.title('Monte Carlo Simulation Results for Index')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend()
    plt.show()

def mass_lum_plots(mass_e, mass_t, mass_sun, lum_sun, lum_e, lum_t):

    """
    Plot the normalized luminosity as a function of normalized mass for three stars: the Sun, Epsilon Eridani, 
    and Theta Persei A.
    
    Parameters:
    mass_e (array-like): Array representing the mass distribution for Epsilon Eridani.
    mass_t (array-like): Array representing the mass distribution for Theta Persei A.
    mass_sun (array-like): Array representing the mass distribution for the Sun.
    lum_sun (array-like): Array representing the luminosity distribution for the Sun.
    lum_e (array-like): Array representing the luminosity distribution for Epsilon Eridani.
    lum_t (array-like): Array representing the luminosity distribution for Theta Persei A.
    
    Returns:
    None: The function displays a plot of normalized luminosity as a function of normalized mass.
    """
    
    # Normalize Mass (using your specified normalization formula)
    mass_sun_normalized = (mass_sun - np.min(mass_sun)) / (np.max(mass_sun) - np.min(mass_sun))
    mass_e_normalized = (mass_e - np.min(mass_e)) / (np.max(mass_e) - np.min(mass_e))
    mass_t_normalized = (mass_t - np.min(mass_t)) / (np.max(mass_t) - np.min(mass_t))

    # Normalize Luminosity (using your specified normalization formula)
    lum_sun_normalized = (lum_sun - np.min(lum_sun)) / (np.max(lum_sun) - np.min(lum_sun))
    lum_e_normalized = (lum_e - np.min(lum_e)) / (np.max(lum_e) - np.min(lum_e))
    lum_t_normalized = (lum_t - np.min(lum_t)) / (np.max(lum_t) - np.min(lum_t))
    
    # Plot
    plt.title(r'Luminosity as a function of Mass')
    plt.plot(mass_sun_normalized, lum_sun_normalized, color='red', linestyle='--', label='Sun')
    plt.plot(mass_e_normalized, lum_e_normalized, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(mass_t_normalized, lum_t_normalized, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.ylabel(r'Lr/L')
    plt.xlabel(r'm/M')
    plt.legend()
    plt.show()

def emissivity_plots(xi_sol, xi_e, xi_t, e_t, epp_t, ecno_t, e_e, epp_e, ecno_e, e_sun, epp_sun, ecno_sun):

    """
    Plot the emissivity profiles as functions of the normalized radial coordinate (xi) for three stars: 
    the Sun, Epsilon Eridani, and Theta Persei A. Includes total emissivity and separate plots for 
    CNO cycle and PP chain contributions.
    
    Parameters:
    xi_sol (array-like): Radial coordinate for the Sun.
    xi_e (array-like): Radial coordinate for Epsilon Eridani.
    xi_t (array-like): Radial coordinate for Theta Persei A.
    e_t (array-like): Total emissivity for Theta Persei A.
    epp_t (array-like): Emissivity by PP chains for Theta Persei A.
    ecno_t (array-like): Emissivity by CNO cycle for Theta Persei A.
    e_e (array-like): Total emissivity for Epsilon Eridani.
    epp_e (array-like): Emissivity by PP chains for Epsilon Eridani.
    ecno_e (array-like): Emissivity by CNO cycle for Epsilon Eridani.
    e_sun (array-like): Total emissivity for the Sun.
    epp_sun (array-like): Emissivity by PP chains for the Sun.
    ecno_sun (array-like): Emissivity by CNO cycle for the Sun.
    
    Returns:
    None: The function displays three plots showing total emissivity, CNO cycle emissivity, and PP chain emissivity 
    as functions of normalized xi.
    """

    xi_sol_normalized = (xi_sol - np.min(xi_sol)) / (np.max(xi_sol) - np.min(xi_sol))
    xi_e_normalized = (xi_e - np.min(xi_e)) / (np.max(xi_e) - np.min(xi_e))
    xi_t_normalized = (xi_t - np.min(xi_t)) / (np.max(xi_t) - np.min(xi_t))

    plt.title(r'Emissivity as a function of $\xi$')
    plt.plot(xi_sol_normalized, e_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, e_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, e_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\epsilon (erg/g/s$')
    plt.legend()
    plt.show()

    plt.title(r'Emissivity by CNO cycle as a function of $\xi$')
    plt.plot(xi_sol_normalized, ecno_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, ecno_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, ecno_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\epsilon_{cno} (erg/g/s$')
    plt.legend()
    plt.show()

    plt.title(r'Emissivty by PP chains as a function of $\xi$')
    plt.plot(xi_sol_normalized, epp_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, epp_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, epp_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\epsilon_{pp}  (erg/g/s$')
    plt.legend()
    plt.show()

def thermo_plots(xi_sol, xi_e, xi_t,  rho_sun, rho_e, rho_t, P_sun, P_e, P_t, T_sun, T_e, T_t):

    """
    Plot the thermodynamic profiles (density, temperature, and pressure) as functions of the normalized 
    radial coordinate (xi) for three stars: the Sun, Epsilon Eridani, and Theta Persei A.
    
    Parameters:
    xi_sol (array-like): Radial coordinate for the Sun.
    xi_e (array-like): Radial coordinate for Epsilon Eridani.
    xi_t (array-like): Radial coordinate for Theta Persei A.
    rho_sun (array-like): Density profile for the Sun.
    rho_e (array-like): Density profile for Epsilon Eridani.
    rho_t (array-like): Density profile for Theta Persei A.
    P_sun (array-like): Pressure profile for the Sun.
    P_e (array-like): Pressure profile for Epsilon Eridani.
    P_t (array-like): Pressure profile for Theta Persei A.
    T_sun (array-like): Temperature profile for the Sun.
    T_e (array-like): Temperature profile for Epsilon Eridani.
    T_t (array-like): Temperature profile for Theta Persei A.
    
    Returns:
    None: The function displays three plots showing density, temperature, and pressure as functions of normalized xi.
    """
    
    # Normalize xi values to range between 0 and 1
    xi_sol_normalized = (xi_sol - np.min(xi_sol)) / (np.max(xi_sol) - np.min(xi_sol))
    xi_e_normalized = (xi_e - np.min(xi_e)) / (np.max(xi_e) - np.min(xi_e))
    xi_t_normalized = (xi_t - np.min(xi_t)) / (np.max(xi_t) - np.min(xi_t))

    plt.title(r'Density as a function of $\xi$')
    plt.plot(xi_sol_normalized, rho_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, rho_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, rho_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\rho (g/cm^3)$')
    plt.legend()
    plt.show()

    plt.title(r'Temperature as a function of $\xi$')
    plt.plot(xi_sol_normalized, T_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, T_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, T_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$T (K)$')
    plt.legend()
    plt.show()

    plt.title(r'Pressure as a function of $\xi$')
    plt.plot(xi_sol_normalized, P_sun, color='red', linestyle='--', label='Sun')
    plt.plot(xi_e_normalized, P_e, color = 'orange', linestyle='--',label='Epsilon Eridani')
    plt.plot(xi_t_normalized, P_t, color = 'blue', linestyle='--',label=r'$\theta$ Persei A')
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$P (Ba)$')
    plt.legend()
    plt.show()

