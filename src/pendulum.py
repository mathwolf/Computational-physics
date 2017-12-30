# pendulum.py
# Python 2.7.2

# Steven Wray
# Physics 4620
# Project 1: Simple pendulum

# Use Euler's method and Verlet's method to approximate the motion
# of a simple pendulum

import math
import numpy as np
import matplotlib.pyplot as plt

EULER = 1
VERLET = 2

def main():

    # Get parameters from user
    print("Numerical methods")
    print("\t1. Euler")
    print("\t2. Verlet")
    s = input("Enter a numerical method: ")
    numerical_method = int(s)
    s = input("Initial angle [degrees]: ")
    theta0 = float(s)
    s = input("Enter time step [seconds]: ")
    tau = float(s)
    s = input("Enter number of time steps: ")
    nstep = int(s)
    
    # Set physical constants and variables
    G_OVER_L = 1
    time = 0        # Initial time
    reversals = 0   # Count reversals of pendulum

    # Set initial position of pendulum
    theta = theta0 * math.pi / 180  # radians
    omega = 0                       # angular velocity
    
    # Take a backward step to start Verlet method
    accel = - G_OVER_L * math.sin(theta)
    theta_old = theta - omega * tau + 0.5 * tau**2 * accel

    # Lists to store data for plotting
    theta_list = [theta]
    time_list = [0]
    time = time + tau

    for i in range(nstep):
        # Record angle and time for plotting

        # Compute new position
        accel = - G_OVER_L * math.sin(theta)
        if numerical_method == EULER:
            theta_old = theta
            theta = theta + tau * omega
            omega = omega + tau * accel

        else:   # Verlet method
            theta_new = 2 * theta - theta_old + tau**2 * accel
            theta_old = theta
            theta = theta_new

        # Record angle and time for plotting
        theta_list.append(theta * 180 / math.pi) # in degrees
        time_list.append(time)

        # See if the pendulum has passed through theta = 0
        # and estimate period
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
            reversals = reversals + 1

        time = time + tau

    # Estimate period of oscillation with error
    if reversals > 1:
        period_sum = 0
        for i in range(reversals - 1):
            print period_list[i]
        average_period = np.average(period_list)
        error = np.std(period_list)
        print "Average period ", round(average_period, 3), " +/- ", \
              round(error, 3), " [s]"

    # Plot the data
    line1 = plt.plot(time_list, theta_list, 'b')
    plt.ylabel('Angle [degrees]')
    plt.xlabel('Time [s]')
    plt.show()
    
main()
