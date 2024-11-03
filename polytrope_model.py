
# Imports
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt

# Lane-Emden
def lane_emden_sys(xi, S, n):
    y1, y2 = S
    dy1_dxi = y2
    dy2_dxi = -2*y2/xi - np.abs(y1)**n

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
    
    dec1 = (y1[-2] - y1[-1])/(x[-2]-x[-1])
    root = -(y1[-2] - dec1*x[-2])/dec1
    dec2 = (y2[-1] - y2[-2])/(x[-1] - x[-2])
    ponto_sup = y2[-2]-dec2*x[-2]+dec2*x[-1]
    x[-1] = root
    y1[-1] = 0
    y2[-1] = ponto_sup

    return x, y1, y2

# Chemical compostion (H and He)
def chemical_composition(xi, X_s, r_c, d_c, R, t, t_tot, Z_s):

    # Hidrogen
    X_c = X_s - X_s * (t/t_tot)
    X = X_c + (X_s - X_c)/(1+np.exp((r_c-xi)/(R*d_c)))
    
    # Helium
    Y = 1 - Z_s - X

    return X, Y

# Thermodynamic functions (general and center)

def density(M, R, n):
    xi, y1, y2 = lane_emden_solver(n)
    a_n = -xi[-1]/(3*y2[-1])
    central_density = M/(4/3*np.pi*R**3) * a_n
    density = central_density * y1**n

    return density, central_density

def pressure(G, M, R, n):
    xi, y1, y2 = lane_emden_solver(n)
    c_n = 1 / (4 * np.pi * (n+1) * y2[-1]**2)
    central_pressure = (G*M**2)/R**4 * c_n
    pressure = central_pressure * y1**(n+1)

    return pressure, central_pressure

def temperature(G, M, R, X, Z, R_g, n): # R_g = gas constant, G  (define both as global variables?)
    xi, y1, y2 = lane_emden_solver(n)
    mu = 4/(5*X-Z+3)
    b_n = 1 / (n+1) * xi[-1] * (-y2[-1])
    central_temperature = mu*G*M/(R_g*R) * b_n
    temperature = central_temperature *  y1

    return temperature, central_temperature

def emissivity(rho, T, X, Z):
    T6 = T/1e6

    T_og_length = len(T6) # original length of T6 array
    T_new_length = 1 + int(np.argwhere(T6<1e-3)[0]) # the remaining values of T6 might cause numerically unstable values of emissivity 
    size_dif = T_og_length - T_new_length

    # Update array sizes for emissivity calculations
    T6 = T6[:T_new_length]
    X = X[:T_new_length]
    rho = rho[:T_new_length]
    e_0 = 2.38e106*X**2*rho

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