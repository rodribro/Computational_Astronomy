
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