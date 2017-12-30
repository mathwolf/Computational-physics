# pendulum.py
# Python 2.7.2

# Steven Wray
# Physics 4620
# Project 1: Simple pendulum

# 2.17: Use Verlet's method to find the period of a pendulum for a range
# of initial angles.  Compare this to the period given by different
# approximations.

import math
import numpy as np
import matplotlib.pyplot as plt


def main():

    # Set physical constants and parameters
    G_OVER_L = 1
    nstep = 300     # Number of timestep for Verlet method
    tau = 0.05      # Timestep

    # Lists to store data for plotting
    theta_list = [0]
    average_period_list = [0]

    for theta0 in range(2, 170):

        # Set initial conditions for pendulum
        time = 0                        # initial time
        reversals = 0                   # count number of reversals
        theta = theta0 * math.pi / 180  # radians
        omega = 0                       # angular velocity
    
        # Take a backward step to start Verlet method
        accel = - G_OVER_L * math.sin(theta)
        theta_old = theta - omega * tau + 0.5 * tau**2 * accel

        for i in range(nstep):
            # Update position with Verlet method
            accel = - G_OVER_L * math.sin(theta)
            theta_new = 2 * theta - theta_old + tau**2 * accel
            theta_old = theta
            theta = theta_new

            # See if the pendulum has passed through theta = 0
            # and estimate period
            if theta * theta_old < 0:
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
            average_period = np.average(period_list)
            error = np.std(period_list)
            if error / average_period > 0.01:
                print ("WARNING: Error in period is greater than 1%")
                print ("\tfor initial angle", theta0)
            # Record data for plotting
            if theta_list[-1] == 0:
                theta_list = [theta0]
                average_period_list = [average_period]
            else:
                theta_list.append(theta0)
                average_period_list.append(average_period)
                
    # Record dataset for plotting small angle approximation of period
    average_period = 2 * math.pi * math.sqrt(1/G_OVER_L)
    small_angle_period = [average_period]
    for theta0 in range(3, 170):
        small_angle_period.append(average_period)

    # Find theta where error in small angle approximation exceeds 10%
    for i in range(len(theta_list)):
        error = math.fabs(average_period_list[i] - small_angle_period[i]) \
                / average_period_list[i]
        if error > 0.1:
            theta = theta_list[i]
            s = 'Error in small angle approximation exceeds 10% for theta = ' \
                + str(theta) + ' degrees.'
            print(s)
            break

    # Record dataset for plotting quadratic approximation of period
    elliptic_period_list = [0]
    for theta0 in range (2,170):
        theta = math.pi * theta0 / 180
        period = 2 * math.pi * math.sqrt(1/G_OVER_L) * (1 + theta**2 / 16)
        if elliptic_period_list[-1] == 0:
            elliptic_period_list = [period]
        else:
            elliptic_period_list.append(period)

    # Find theta where error in quadratic approximation exceeds 10%
    for i in range(len(theta_list)):
        error = math.fabs(average_period_list[i] - elliptic_period_list[i]) \
                / average_period_list[i]
        if error > 0.1:
            theta = theta_list[i]
            s = 'Error in elliptic integral approximation exceeds 10% for theta = ' \
                + str(theta) + ' degrees.'
            print(s)
            break

    
        
    # Plot the datasets
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(theta_list, average_period_list, 'b', \
            theta_list, small_angle_period, 'b--', \
            theta_list, elliptic_period_list, 'b:')
    leg = ax.legend(("Verlet method", "Small angle (constant) approximation", \
                     "Quadratic approximation"), 'upper left', shadow=True)
    ax.set_ylim([0,16])
    ax.grid(False)
    ax.set_xlabel('Initial Angle [degrees]')
    ax.set_ylabel('Period [s]')

    plt.show()
    
main()
