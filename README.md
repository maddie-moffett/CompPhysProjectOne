# Computational Physics Project One
 ## Introduction
This project models a double pendulum. This type of pendulum consists of two pendulums attached to each other-- one hangs from the bob of the other. This double pendulum, in particular, has a sturdy rod connecting the masses. The double pendulum is a chaotic system that has unpredictable motion. The motion of the two arms is described by four coupled, nonlinear differential equations describing its angular positions and angular velocities.  

To do this, a fourth-order Runge-Kutta algorithm is utilized. This is a standard numerical technique for solving differential equations by stepping forward in small increments of time, estimating the rate of change at each step. The smaller the time step, the more accurate the result is.  

One way to check the accuracy of the numerical solution is to track the total mechanical energy of the system. The total mechanical energy is the sum of the kinetic and potential energies of both masses. The kinetic energy can vary because the velocity of the second mass depends on the motion of both arms, which introduces a coupling term. The potential energy is gravitational potential energy of each mass relative to the pivot point. Because this is a conservative system, the total energy remains constant throughout the motion of the system.  

One of the key features of this system is its sensitivity to initial conditions. For small angles, the nonlinear coupling terms become negligible, and the system behaves like two independent linear pendulums with regular, predictable motion. However, for the larger initial angles, the nonlinear terms dominate and the system exhibits chaotic behavior, such as trajectories diverging exponentially in time which makes long-term prediction impossible.  

This project investigates how sensitive the double pendulum is to its initial conditions utilizing the fourth-order Runge-Kutta algorithm and total mechanical energy.  
## Numerical Methods
To solve the four coupled equations of motion, a fourth-order Runge-Kutta algorithm was implemented using Python code. The algorithm was coded manually to maintain full control over the integration process and to better understand the numeric. 

The state of the system at any given time is represented as a vector of four quantities: θ₁, θ₂, ω₁, and ω₂. A function was written to evaluate the time derivatives of this given state at any moment. This derivative function takes the current state and returns the instantaneous rates of change of each variable.   

The Runge-Kutta algorithm advances the state forward by one time step, h, by computing four estimates of the slope (k₁, k₂, k₃, and k₄). The first estimate, k₁, uses the slope at the beginning of the interval. The second and third estimates, k₂ and k₃, use the slope at the midpoint of the interval, each informed by the previous estimate. The fourth estimate, k₄, uses the slope at the end of the interval. The new state is then computed as a weighted average of these four estimates, with the midpoint sloped recieving twice the weight of the endpoints. This weighting gives the Runge-Kutta algorithm its high accuracy. To support this, two small helper functions were written which lists a scalar and adding an element list that allows the k values to be computed and combined.  

The integration process proceeds by stepping forward from t = 0 to t = 100 seconds in increments of h, storing the values of θ₁, θ₂, ω₁, and ω₂ at each time step. From these angles and angular velocities, three quantities are computed and plotted: the total mechanical energy, x position of the masses, and the y position of the masses. The Carestian coordinates are obtained from the geometric relations for a double pendulum, where the first position only depends on θ₁ and the position of the second mass depends on both θ₁ and θ₂. For the plots, only every tenth data point is used when plotting.  

The total energy is computed at each sampled time step by evaluating the kinetic and potential energy expressions. Because the system is conservative, any change in energy over time is a result of numerical error caused during integration. The difference between the initial and final energy values is used a direct measure of accuracy of the numerical solution.  
## Results
