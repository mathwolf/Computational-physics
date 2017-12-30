# Python 2.7.2

# Steven Wray
# Physics 4620/HW2
# Lotka-Volterra model of predator/prey relationship

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as fig
from mpl_toolkits.mplot3d import Axes3D

STEADY_STATE_FLAG = 0

def rk4(x, t, tau, derivsRK, param):
    # 4th order Runge-Kutta integrator
    # Input
    #   x           current value of dependent variable
    #   t           independent variable
    #   tau         stepsize
    #   derivsRK    righthand side of ODE
    #               returns dx/dt
    #               format derivsRK (x, param)
    #   param       extra parameters passed to derivsRK
    # Output
    #   xout        new value of x after step of size tau

    half_tau = 0.5 * tau
    f1 = derivsRK(x, t, param)
    xtemp = x + half_tau * f1
    f2 = derivsRK(xtemp, t, param)
    xtemp = x + half_tau * f2
    f3 = derivsRK(xtemp, t, param)
    xtemp = x + tau * f3
    f4 = derivsRK(xtemp, t, param)
    xout = x + (tau/6.) * (f1 + f4 + 2.*(f2+f3))
    return xout


def rka(x, t, tau, err, derivsRK, param):
    # Adaptive Runge-Kutta routine
    # Inputs
    #   x       current value of dependent variable
    #   t       independent variable (time)
    #   tau     step size
    #   err     desired local truncation error
    #   derivsRK    righthand side of ODE
    #               returns dx/dt
    #   param   extra parameters passed to derivsRK
    # Outputs
    #   x_small  new value of dependent variable
    #   t       new value of independent variable
    #   tau     suggested step size for next call to rka

    # Save initial values
    t_save = t
    x_save = x
    # Safety factors
    safe1 = 0.9
    safe2 = 4.

    # Loop over max attempts to satisfy error bound
    max_try = 100
    for i in range(max_try):
        # Take two small time steps
        half_tau = 0.5 * tau
        x_temp = rk4(x_save, t_save, half_tau, derivsRK, param)
        t = t_save + half_tau
        x_small = rk4(x_temp, t, half_tau, derivsRK, param)

        # Take single big time step
        t = t_save + tau
        x_big = rk4(x_save, t_save, tau, derivsRK, param)

        # Compute the estimated truncation error
        error_ratio = 0
        eps = 10 ** (-16)    
        for j in range(x_big.size):
            scale = err * (abs(x_small[j]) + abs(x_big[j])) / 2.
            x_diff = x_small[j] - x_big[j]
            ratio = abs(x_diff) / (scale + eps)
            if ratio > error_ratio:
                error_ratio = ratio

        # Estimate new tau value
        tau_old = tau
        if error_ratio > 0.:
            tau = safe1 * tau_old * error_ratio ** (-0.20)
        else:
            tau = 0
        tau = max(tau, tau_old/safe2)
        tau = min(tau, safe2*tau_old)

        # If error is acceptable, return computed values
        if error_ratio < 1:
            return x_small, t, tau

    # Issue error message if bound is never satisfied
    print("Error: Adaptive RK routine failed.")
    sys.exit()


def lotkaVolterraRK(s, t, param):
    # Returns right hand side of Lotka-Volterra ODEs
    # Inputs
    #   s           state vector [x, y]
    #   t           time
    #   param       parameter vector [a, b, c, d, e]
    # Output
    #   deriv       state vector [x', y']

    # Unravel input vectors for clarity
    x = s[0]
    y = s[1]
    a = param[0]
    b = param[1]
    c = param[2]
    d = param[3]
    e = param[4]

    # Compute derivatives
    x_prime = (a - b*x - c*y) * x
    y_prime = (-d + e*x) * y

    # Return state vector of derivatives
    deriv = np.array([x_prime, y_prime], dtype = float)
    return(deriv)


def main():
    # Compute trajectories for the Lotka-Volterra equations using
    # an adaptive Runge-Kutta routine

    # Set initial state of system
    s = input("Enter the initial x value (prey population): ")
    x = float(s)
    s = input("Enter the initial y value (predator population): ")
    y = float(s)
    state = [x, y]

    # Set the parameters 
    A = 10.0
    B = 10.0**(-5)
    C = 0.1
    D = 10.0
    E = 0.1
    param = [A, B, C, D, E]

    tau = 1     # Initial guess for timestep
    err = 0.001 # Error tolerance for adaptive routine

    # Set the number of timesteps
    time = 0
    s = input("Enter desired number of steps: ")
    nstep = int(s)

    # Record initial values for plot
    x_plot_list = [state[0]]
    y_plot_list = [state[1]]
    time_plot_list = [time]
    tau_plot_list = [tau]

    # Loop
    for i in range(nstep):
        # Update state using adaptive Runge-Kutta
        state, time, tau = rka(state, time, tau, err, lotkaVolterraRK, param)

        if i % 10 == 0:
            s = 'Finished ' + str(i) + ' steps out of ' + str(nstep) + '.'
            print(s)

        # Record for plotting
        x_plot_list.append(state[0])
        y_plot_list.append(state[1])
        time_plot_list.append(time)
        tau_plot_list.append(tau)
        
    # Print max and min steps returned by rka
    s = 'Adaptive time step: Max = ' + str(max(tau_plot_list)) + \
        ', Min = ' + str(min(tau_plot_list)) + '\n'
    print(s)

    # Plot x(t) versus y(t)
    plt.plot(x_plot_list, y_plot_list, '-')
    plt.xlabel('Prey population')
    plt.ylabel('Predator population')
    plt.title('Lotka-Volterra model phase plot')
    plt.show()
    

        
        
main()
