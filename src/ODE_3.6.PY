# Program to compute the orbit of a comet

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as fig

# Set physical parameters
GM = 4*(math.pi**2)         # Newton's gravitational constant times
                            # mass of the sun (au^3/yr^2)
MASS_COMET = 1.             # Mass of the comet

ADAPT_ERR = 0.001           # Error constant for adaptive RK
MAX_ITER = 100000           # Maximum iterations when looking for aphelion

# Constants representing different numerical methods
EULER = 1
EULER_CROMER = 2
RUNGE_KUTTA = 3
ADAPTIVE_RK = 4


def norm(v):
    # Calculate the norm of a vector
    # Inputs
    #   v   a vector
    #Output
    #   n    vector norm as float
    n = math.sqrt(np.dot(v,v))
    return(n)


def gravrk(s, t, GM):
    # Return right hand side of Kepler ODE
    # Called by RK routine
    # Inputs
    #   s   state vector [r[0] r[1] v[0] v[1]]]
    #   t   time (not used in this model)
    #   GM  parameter GM (gravitational constant * solar mass)
    # Output
    #   deriv   state vector [dr[0]/dt dr[1]/dt dv[0]/dt dv[1]/dt]

    # Create position and velocity vectors
    r = np.array([s[0], s[1]], dtype=float)
    v = np.array([s[2], s[3]], dtype=float)
    # Compute acceleration
    a = (-1 * GM / (norm(r)**3)) * r
    # Return state vector of derivatives
    deriv = np.array([v[0], v[1], a[0], a[1]], dtype=float)
    return deriv


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
        error_ratio = 0.
        eps = 10 ** (-16)    
        for j in range(x_big.size):
            scale = err * (abs(x_small[j]) + abs(x_big[j])) / 2.
            x_diff = x_small[j] - x_big[j]
            ratio = abs(x_diff) / (scale + eps)
            if ratio > error_ratio:
                error_ratio = ratio

        # Estimate new tau value
        tau_old = tau
        tau = safe1 * tau_old * error_ratio ** (-0.20)
        tau = max(tau, tau_old/safe2)
        tau = min(tau, safe2*tau_old)

        # If error is acceptable, return computed values
        if error_ratio < 1:
            return x_small, t, tau

    # Issue error message if bound is never satisfied
    print("Error: Adaptive RK routine failed.")
    sys.exit()
    


def orbit(initial_distance, initial_velocity, n_iter, tau, \
          numerical_method):
    # Calculate the orbit of a comet around the sun
    # Produce a list of polar coordinates describing the
    # comet's position.  Also produce a list of kinetic and
    # potential energy for the points of the orbit.
    
    # Inputs
    #   initial_distance    starting position of comet (AU)
    #   initial_velocity    starting velocity (AU/yr)
    #   n_iter              total data points to calculate -
    #                       if n_iter = -1, compute orbit to
    #                       next aphelion
    #   tau                 timestep (yr)
    #   numerical_method    one of:
    #                       EULER
    #                       EULER_CROMER
    #                       RUNGE_KUTTA
    #                       ADAPTIVE_RK
    # Outputs - arrays of values. One value for each timestep.
    #   radius_list         polar coordinate R
    #   theta_list          polar coordinate theta
    #   kinetic_list        kinetic energy
    #   potential_list      gravitational potential energy
    
    time = 0

   # Create vectors for distance, velocity, acceleration
    r = np.array([initial_distance, 0.0], dtype=float)
    v = np.array([0.0, initial_velocity], dtype=float)
    a = np.array([0.0, 0.0], dtype=float)
    state = np.array([r[0], r[1], v[0], v[1]], dtype=float)
    
    # Lists of polar coordinates used in plotting
    radius_list = [norm(r)]
    theta_list = [math.atan2(r[1], r[0])]

    # Lists of time/energy values for plotting
    time_list = [time]
    kinetic_list = [0.5 * MASS_COMET * norm(v)**2]
    potential_list = [-1 * GM * MASS_COMET / norm(r)]

    i = 1   # Total data points calculated

    while True:
        # Calculate the values for the next step using selected
        # numerical method
        if numerical_method == EULER:
            a = (-1 * GM / (norm(r)**3)) * r
            r = r + tau*v
            v = v + tau*a
            time = time + tau
        elif numerical_method == EULER_CROMER:
            a = (-1 * GM / (norm(r)**3)) * r
            v = v + tau*a
            r = r + tau*v
            time = time + tau
        elif numerical_method == RUNGE_KUTTA:
            state = rk4(state, time, tau, gravrk, GM)
            r = [state[0], state[1]]
            v = [state[2], state[3]]
            time = time + tau
        else:   # ADAPTIVE_RK
            state, time, tau = rka(state, time, tau, ADAPT_ERR, gravrk, GM)
            r = [state[0], state[1]]
            v = [state[2], state[3]]
        
        # Write new values to list
        radius_list.append(norm(r))
        theta_list.append(math.atan2(r[1], r[0]))
        time_list.append(time)
        kinetic_list.append(0.5 * MASS_COMET * norm(v)**2)
        potential_list.append(-1 * GM * MASS_COMET/ norm(r))

        i = i + 1

        if n_iter == -1:
            # If we have reached aphelion, break out of loop
            if i >= 3 and \
               radius_list[-1] - radius_list[-2] < 0 and \
               radius_list[-2] - radius_list[-3] > 0:
                break
            if i > MAX_ITER:
                # If we have not returned to apehlion at this point,
                # the orbit is probably not closed
                radius_list[0] = -1     # Set flag for hyperbolic orbit
                break
        else:
            if i > n_iter:
                break

    return(radius_list, theta_list, time_list, kinetic_list, potential_list)


def plot_orbit():
    # Plot the orbit and energy for an arbitrary comet.
    
    # Get initial position and velocity from user
    s = input("Enter initial radial distance (AU): ")
    initial_distance = float(s)
    s = input("Enter initial tangential velocity (AU/yr): ")
    initial_velocity = float(s)

    # Get total interations, timestep, and numerical method from user
    s = input("Enter number of steps: ")
    n_iter = int(s)
    s = input("Enter time step (yr): ")
    tau = float(s)
    print("Numerical methods:")
    print("\t1. Euler")
    print("\t2. Euler-Cromer")
    print("\t3. Runge-Kutta")
    print("\t4. Adaptive RK")
    s = input("Enter numerical method (1-4): ")
    numerical_method = int(s)

    r_list, theta_list, time_list, kinetic_list, potential_list = \
            orbit(initial_distance, initial_velocity, n_iter, tau, \
                  numerical_method)

    # Create list for total energy plot
    total_energy_list = []
    for i in range(len(kinetic_list)):
        total_energy_list.append(kinetic_list[i] + potential_list[i])
    
    # Plot orbit
    plt.figure(1)
    plt.subplot(211, polar='True')
    plt.plot(theta_list, r_list, 'ro')
    plt.title('Orbital Path')

    plt.subplot(212)
    line1 = plt.plot(time_list, kinetic_list, 'b--')
    line2 = plt.plot(time_list, potential_list, 'g--')
    line3 = plt.plot(time_list, total_energy_list, 'k')
    plt.ylabel('Energy  [M AU^2/yr^2]')
    plt.xlabel('Time [yr]')
    plt.figlegend((line1, line2, line3), ('kinetic', 'potential', \
                  'total'), 'lower right')
    plt.title('Energy')

    plt.show()


def find_timestep():
    # Get input from user
    print("\nUsing a comet with an aphelion of 35 AU.")
    print("Based on aphelion of Halley's comet.")
    s = input("Enter velocity at aphelion (AU/yr): ")
    initial_velocity = float(s)
    print("Numerical methods:")
    print("\t1. Euler-Cromer")
    print("\t2. Runge-Kutta")
    s = input("Enter numerical method (1-4): ")
    if s == 1:
        numerical_method = EULER_CROMER
    else:
        numerical_method = RUNGE_KUTTA
        
    for i in range (40,0,-1):
        tau = i
        print(tau)
        r_list, theta_list, time_list, kinetic_list, potential_list = \
                orbit(35, initial_velocity, -1, tau, numerical_method)
        if r_list[0] == -1:     # Calculated orbit is not closed
            continue
        initial_energy = kinetic_list[0] + potential_list[0]
        final_energy = kinetic_list[-1] + potential_list[-1]
        if abs(initial_energy - final_energy)/initial_energy < 0.01:
            s = '\nTimestep of ' + str(tau) + ' years required for 1% error.\n'
            print(s)
            return

    
    for i in range(1,6):
        for j in range (9,0,-1):
            tau = j * 10**(-1 * i)
            print(tau)
            r_list, theta_list, time_list, kinetic_list, potential_list = \
                    orbit(35, initial_velocity, -1, tau, numerical_method)
            if r_list[0] == -1:     # Calculated orbit is not closed
                continue
            initial_energy = kinetic_list[0] + potential_list[0]
            final_energy = kinetic_list[-1] + potential_list[-1]
            if abs(initial_energy - final_energy)/initial_energy < 0.01:
                s = '\nTimestep of ' + str(tau) + ' years required for 1% error.\n'
                print(s)
                return

    print("\nTimestep smaller than 10^-5 required or orbit not closed.\n")
    return
            
    


def main():
    while True:
        # Ask user for desired function
        print("Functions:")
        print("\t1. Plot the orbit of an arbitrary comet.")
        print("\t2. Find the largest timestep that conserves total energy ")
        print("\t\tto 1% per orbit. (Exercise 3.6, 3.11)")
              
        print("\t9. Exit")
        s = input("Enter a number (1-9): ")
        n = int(s)
        if n == 1:
            plot_orbit()
        elif n==2:
            find_timestep()
        else:       # Exit program
            break
 

main()
