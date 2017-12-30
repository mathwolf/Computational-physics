# pendulum.py
# Python 2.7.2

# Steven Wray
# Physics 4620
# Project 1: Simple pendulum

# 2.20: Improve the estimate of the pendulum's period by using
# interpolation.

import math
import numpy as np
import matplotlib.pyplot as plt

# Function to interpolate between three data points
# using a quadratic Lagrange polynomial.

def interp(xi, x, y):
    # x and y are arrays that contain the coordinates for the three
    # data points
    # xi: interpolate the function at this value
    yi = (xi-x[1]) * (xi-x[2]) / ((x[0]-x[1]) * (x[0]-x[2])) * y[0] \
         + (xi-x[0]) * (xi-x[2]) / ((x[1]-x[0]) * (x[1]-x[2])) * y[1] \
         + (xi-x[0]) * (xi-x[1]) / ((x[2]-x[0]) * (x[2]-x[1])) * y[2] 
    return yi

EULER = 1
VERLET = 2
EULER_CROMER = 3
MIDPOINT = 4
LEAPFROG = 5

def main():

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

        elif numerical_method == VERLET:   
            theta_new = 2 * theta - theta_old + tau**2 * accel
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
        time_list.append(time)

        # See if the pendulum has passed through theta = 0
        # and estimate period
        if theta * theta_old < 0:
            # Make sure we have at least three data points for interpolation
            # Otherwise just use the current time as turning point
            if (len(theta_list) >= 3):
                # Estimate turning point using interpolation
                turning_point = interp(0.0, \
                                   [theta_list[-3], theta_list[-2], theta_list[-1]],
                                   [time_list[-3], time_list[-2], time_list[-1]])
            else:
                turning_point = time
            
            print "Turning point at time t =", round(turning_point, 4), " [s]."
            if reversals == 0:
                turning_point_old = turning_point
            elif reversals == 1:
                period_list = [2*(turning_point - turning_point_old)]
                turning_point_old = turning_point
            else:
                period_list.append(2*(turning_point - turning_point_old))
                turning_point_old = turning_point
            reversals = reversals + 1

        time = time + tau

    # Estimate period of oscillation with error
    if reversals > 1:
        period_sum = 0
        average_period = np.average(period_list)
        error = np.std(period_list)
        print "Average period ", round(average_period, 5), " +/- ", \
              round(error, 7), " [s]"

    # Plot the data
    line1 = plt.plot(time_list, theta_list, 'b')
    plt.ylabel('Angle [degrees]')
    plt.xlabel('Time [s]')
    plt.show()
    
main()
