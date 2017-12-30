# pendulum.py
# Python 2.7.2

# Steven Wray
# Physics 4620
# Project 1: Simple pendulum

# 2.18: Use a numerical method to create a phase space plot
# of the pendulum's motion.

import math
import numpy as np
import matplotlib.pyplot as plt

EULER = 1
VERLET = 2
EULER_CROMER = 3
MIDPOINT = 4
LEAPFROG = 5

def main():

    # Get parameters from user
    print("Numerical methods")
    print("\t1. Euler")
    print("\t2. Verlet")
    print("\t3. Euler-Cromer")
    print("\t4. Midpoint")
    print("\t5. Leapfrog")
    s = input("Enter a numerical method: ")
    numerical_method = int(s)
    s = input("Initial angle [degrees]: ")
    theta0 = float(s)
    s = input("Enter time step [seconds]: ")
    tau = float(s)
    
    # Set physical constants and variables
    G_OVER_L = 1
    time = 0            # Initial time
    reversals = 0       # Count reversals of pendulum
    max_steps = 10000    # Maximum timesteps

    # Set initial position of pendulum
    theta = theta0 * math.pi / 180  # radians
    omega = 0                       # angular velocity
    
    # Take a backward step to start Verlet
    # and leapfrog methods
    accel = - G_OVER_L * math.sin(theta)
    if numerical_method == VERLET:
        theta_old = theta - omega * tau + 0.5 * tau**2 * accel
    if numerical_method == LEAPFROG:
        omega = omega - tau * accel

    # Lists to store data for plotting
    theta_list = [theta0]
    omega_list = [omega]
    time = time + tau

    i = 0
    while i < max_steps:
        # Compute new position and velocity
        accel = - G_OVER_L * math.sin(theta)
        if numerical_method == EULER:
            theta_old = theta
            theta = theta + tau * omega
            omega = omega + tau * accel

        elif numerical_method == VERLET:   
            theta_new = 2 * theta - theta_old + tau**2 * accel
            omega = (theta_new - theta_old) / (2 * tau)
            theta_old = theta
            theta = theta_new

        elif numerical_method == EULER_CROMER:
            omega = omega + tau * accel
            theta_old = theta
            theta = theta + tau * omega

        elif numerical_method == MIDPOINT:
            omega_new = omega + tau *accel
            theta_old = theta
            theta = theta + tau * (omega_new + omega) / 2
            omega = omega_new

        else:   # Leapfrog
            omega = omega + 2 * tau * accel
            time = time + tau
            theta_old = theta
            theta = theta + 2 * tau * omega

        # Record angle and time for plotting
        theta_list.append(theta * 180 / math.pi) # in degrees
        omega_list.append(omega * 180 / math.pi)

        # See if the pendulum has passed through theta = 0
        # and estimate period.  Exit the loop if the pendulum
        # has finished one complete period.
        if theta * theta_old < 0:
            print "Turning point at time t =", round(time, 4), " [s]."
            if reversals == 0:
                time_old = time
            elif reversals == 1:
                period_list = [2*(time - time_old)]
                time_old = time
            else:
                period_list.append(2*(time - time_old))
                time_old = time
                break
            reversals = reversals + 1

        time = time + tau
        i = i + 1

    # If we are using the Verlet method we need one additional value
    # of omega to make our lists match up.
    if numerical_method == VERLET:
        theta_new = 2 * theta - theta_old + tau**2 * accel
        omega = (theta_new - theta_old) / (2 * tau)
        omega_list.pop(0)
        omega_list.append(omega * 180 / math.pi)


    # Estimate period of oscillation with error
    if reversals > 1:   # Make sure there is at least one data point
        period_sum = 0
        average_period = np.average(period_list)
        error = np.std(period_list)
        print "Average period ", round(average_period, 5), " +/- ", \
              round(error, 5), " [s]"

    # Plot the data
    # Plot the datasets
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(theta_list, omega_list, 'b')
    ax.set_ylim([-theta0 - 5, theta0 + 5])
    ax.set_xlim([-theta0 - 5, theta0 + 5])
    ax.grid(False)
    ax.set_xlabel('Angle [degrees]')
    ax.set_ylabel('Angular velocity [deg/s]')
    plt.show()
    
main()
