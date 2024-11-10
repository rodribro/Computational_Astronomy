
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

# Lane-Emden
def lane_emden_sys(xi, S, n):
    y1, y2 = S
    dy1_dxi = y2
    dy2_dxi = - np.abs(y1)**n - 2*y2/xi 

    return [dy1_dxi, dy2_dxi]



def lane_emden_solver(n):

    # xi Range
    xi_0 = 1e-7 # prevent numerical singularity
    xi_range = np.linspace(xi_0, 40, 100000)


    # Initial Conditions
    y1_0 = 1-xi_0**2/2           # θ(ξ₀) ≈ 1 - ξ₀²/2
    y2_0 = -xi_0                        # dθ/dξ(ξ₀) ≈ -ξ₀
    lane_em_0 = [y1_0, y2_0]

    # Solve Lane-Emden system
    sol = odeint(lane_emden_sys, lane_em_0, xi_range, args=(n,), tfirst=True)
    
    x = np.append([0], xi_range)
    y1 = np.append([1], sol[:, 0])
    y2 = np.append([0], sol[:, 1])

    zero_index = 1 + int(np.argwhere(y1<0)[0])

    x = x[:zero_index]
    y1 = y1[:zero_index]
    y2 = y2[:zero_index]

    # Interpolação linear para encontrar o ponto onde y1 = 0
    interp_func_y1 = interp1d(y1[-2:], x[-2:], kind='linear')
    root = interp_func_y1(0)  # O valor de x onde y1 cruza zero

    # Atualizar x e y1 com o ponto onde y1 = 0
    x[-1] = root
    y1[-1] = 0

    # Para y2, fazer uma interpolação semelhante para obter o valor correspondente de y2
    interp_func_y2 = interp1d(x[-2:], y2[-2:], kind='linear')
    y2[-1] = interp_func_y2(root)  # Obter o valor interpolado de y2 no novo ponto x

    return x, y1, y2

# Chemical compostion (H and He)
def chemical_composition(xi, X_s, r_c, d_c, t, t_tot):

    # Hidrogen
    X_c = X_s - X_s * t/t_tot
    X = X_c + (X_s - X_c)/(1+np.exp((r_c-xi)/xi[-1]/(d_c)))
    # X = X_c + (X_s - X_c)/(1+np.exp((r_n - xi/xi[-1])/(d_c)))
    

    return X

# Thermodynamic functions (general and center)
'''
def density(M, R, n):
    xi, y1, y2 = lane_emden_solver(n)
    central_density = -(M*xi[-1])/(4*np.pi*y2[-1]*R**3)
    density = central_density * y1**n

    return density, central_density

def pressure(M, R, n):
    xi, y1, y2 = lane_emden_solver(n)
    central_pressure = (G*M**2)/(R**4 * 4*np.pi*(n+1)*y2[-1]**2)
    pressure = central_pressure * y1**(n+1)

    return pressure, central_pressure

def temperature(M, R, X, Z, n): # R_g = gas constant, G  (define both as global variables?)
    xi, y1, y2 = lane_emden_solver(n)
    mu = 4/(5*X-Z+3) # needs to come from chemical composition array 
    central_temperature = mu*G*M/(R_g*R*(n+1)*xi[-1]*y2[-1])
    temperature = central_temperature *  y1

    return temperature, central_temperature

    GENERATING OVERFLOW PROBLEMS FOR LUMINOSITY
'''

def rho_P_T_central(n, R, M, X, Z):
    x, y1, y2 = lane_emden_solver(n)
    mu = 4/(5*X-Z+3)
    rho_c = -(M*x[-1])/(4*np.pi*y2[-1]*R**3)
    T_c = -(mu*G*M)/(R_g*R*(n+1)*x[-1]*y2[-1])
    P_c = (G*M**2)/(R**4 * 4*np.pi*(n+1)*y2[-1]**2)
    return rho_c, P_c, T_c

def rho_P_T(n, R, M, X, Z):
    x, y1, y2 = lane_emden_solver(n)
    mu = 4/(5*X-Z+3)
    rho_c = -(M*x[-1])/(4*np.pi*y2[-1]*R**3)
    T_c = -(mu*G*M)/(R_g*R*(n+1)*x[-1]*y2[-1])
    P_c = (G*M**2)/(R**4 * 4*np.pi*(n+1)*y2[-1]**2)
    rho = rho_c* y1**n
    T = T_c * y1
    P = P_c * y1**(n+1)
    return rho, P, T

def emissivity(rho, T, X, Z):
    T6 = T/1e6

    T_og_length = len(T6) # original length of T6 array
    T_new_length = 1 + int(np.argwhere(T6<1e-3)[0]) # the remaining values of T6 might cause numerically unstable values of emissivity 
    size_dif = T_og_length - T_new_length

    # Update array sizes for emissivity calculations
    T6 = T6[:T_new_length]
    X = X[:T_new_length]
    rho = rho[:T_new_length]

    # Formula variables
    alpha = 1.2e17*((1-X-Z)/(4*X))**2 * np.exp(-100*T6**(-1/3))
    eps_0 = 2.38e6* X** 2 * rho* T6** (-2/3) * (1 + 0.0123*T6** (1/3) + 0.0109*T6** (2/3) + 0.00095*T6) * np.exp(-33.80*T6** (-1/3) + 0.27*rho** (1/2)* T6**(-3/2))
    phi_alpha = 1-alpha + np.sqrt(alpha*(alpha+2))
    F1 = (np.sqrt(alpha+2) - np.sqrt(alpha))/(np.sqrt(alpha+2) + 3*np.sqrt(alpha))
    F2 = (1-F1)/(1+8.94e15 * (X/(4-3*X))* T6**(-1/6) * np.exp(-102.65*T6**(-1/3)))
    F3 = 1-F1-F2

    # Emissivity: proton-proton (e_pp), carbon-nitrogen-oxigen (e_cno)
    e_pp = (eps_0/0.980) * phi_alpha * (0.980*F1 + 0.960*F2 + 0.721*F3)
    e_cno = 8.67e27*Z*X*rho*T6** (-2/3) * (1 + 0.0027*T6** (1/3) - 0.00778*T6** (2/3) - 0.000149*T6) * np.exp(-152.28*T6**(-1/3))

    # Fill the arrays back to their original sizes with zeros (further calculations using emissivity won't be affected)
    zeros = np.zeros(size_dif)
    e_pp = np.concatenate((e_pp, zeros))
    e_cno = np.concatenate((e_cno, zeros))

    # Total emissivity
    e_total = e_pp + e_cno

    return e_total, e_pp, e_cno

def luminosity(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol ):

    # xi, theta(xi), theta'(xi)
    xi, y1, y2 = lane_emden_solver(n)

    # chemical composition X, Y
    X = chemical_composition(xi, X_s, r_c, d_c, time, total_time)

    # Thermodynamic functions
    rho, P, T = rho_P_T(n, R, M, X, Z)

    # Emissivity
    ems_tot, e_pp, e_cno = emissivity(rho, T, X, Z)

    # Constant of the luminosity formula
    const = -1/(xi[-1]**2 * y2[-1])

    # Integrand of the luminosity formula
    integrand = xi**2 * y1**n * M * ems_tot/L_sol

    # Cumulative integral
    integral = cumulative_simpson(integrand, x=xi, initial=0)
    return const * integral


def y(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol, L):
    lumin = luminosity(n, X_s, r_c, d_c, R, time, total_time, M, Z, L_sol)[-1]
    y = np.log10(lumin/(L/L_sol))
    return y

def index(f, a, b, M, L, Y, Z_X, R, L_sol, r_c, d_c, time, total_time):
    X_s = (1-Y)/(1+Z_X)
    Z = 1-Y-X_s 
    index = bisect(f, a, b, args=(X_s, r_c, d_c, R, time, total_time, M, Z, L_sol, L))
    return index

def monte_carlo(M, dM, R, dR, L, dL, Z_X, dZ_X, Y, dY, r_n, d_c, tau, tau_t, N, L_sol,a,b):
    nn = np.zeros(N)
    Mm = M + dM*np.random.normal(0, 1, N)
    Rr = R + dR*np.random.normal(0, 1, N)
    Ll = L + dL*np.random.normal(0, 1, N)
    Z_Xx = Z_X + dZ_X*np.random.normal(0, 1, N)
    Yy = Y + dY*np.random.normal(0, 1, N)
    for i in range(N):
        print(f"Iteration {i + 1}/{N}")
        try:
            # Attempt to calculate the polytropic index with the current parameters
            nn[i] = index(y, a, b, Mm[i], Ll[i], Yy[i], Z_Xx[i], Rr[i], r_n, d_c, tau, tau_t, L_sol)
        except ValueError as e:
            # If there's a ValueError, skip this iteration and print a message
            print(f"Warning: Skipping iteration {i + 1} due to error: {e}")
            #nn[i] = np.nan  # Optionally, set this entry to NaN or some other marker
    return nn


# Lane-Emden visualizations
def lane_emden_plots():
    # List of n values
    ns = [0,1, 1.5, 3, 4, 4.5]

    # Lane-Emden solution for different values of n
    for n in ns:
        xi_range, y1, y2 = lane_emden_solver(n)
        plt.plot(xi_range, y1, label=f"n = {n}")
        plt.ylim(0,1)
        plt.xlim(0,15)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$\theta(\xi)$')
        plt.title('Lane-Emden Solutions for different values of n')
        plt.legend()
        plt.grid(True)
    plt.show()

    # xi_s for different values of n
    for n in ns:
        xi_range, y1, y2 = lane_emden_solver(n)
        xs = xi_range/xi_range[-1]
        plt.plot(xs, y1, label=f"n = {n}")
        plt.ylim(0,1)
        plt.xlim(0,1)
        plt.xlabel(r'$\xi$/$\xi_S$')
        plt.ylabel(r'$\theta$')
        plt.title(r'$\xi_S$ for different values of n')
        plt.legend()
        plt.grid(True)
    plt.show()

# Monte Carlo visualization
def monte_carlo_plots(mean, std, N, nn): 
    # Gaussian
    # Generate x values for the curve
    x_gauss = np.linspace(mean - 5*std, mean + 5*std, N)
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