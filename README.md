# Computational Physics Project One
 ## Introduction
This project models a double pendulum. This type of pendulum consists of two pendulums attached to each other-- one hangs from the bob of the other. This double pendulum, in particular, has a sturdy rod connecting the masses. The double pendulum is a chaotic system that has unpredictable motion. The motion of the two arms is described by four coupled, nonlinear differential equations describing its angular positions and angular velocities.  

To do this, a fourth-order Runge-Kutta algorithm is utilized. This is a standard numerical technique for solving differential equations by stepping forward in small increments of time, estimating the rate of change at each step. The smaller the time step, the more accurate the result is.  

One way to check the accuracy of the numerical solution is to track the total mechanical energy of the system. The total mechanical energy is the sum of the kinetic and potential energies of both masses. The kinetic energy can vary because the velocity of the second mass depends on the motion of both arms, which introduces a coupling term. The potential energy is gravitational potential energy of each mass relative to the pivot point. Because this is a conservative system, the total energy remains constant throughout the motion of the system.  

One of the key features of this system is its sensitivity to initial conditions. For small angles, the nonlinear coupling terms become negligible, and the system behaves like two independent linear pendulums with regular, predictable motion. However, for the larger initial angles, the nonlinear terms dominate and the system exhibits chaotic behavior, such as trajectories diverging exponentially in time which makes long-term prediction impossible.  

This project investigates how sensitive the double pendulum is to its initial conditions utilizing the fourth-order Runge-Kutta algorithm and total mechanical energy.  
