import pylab
from math import cos, sin

# Standards
# l = l1 = l2
# m = m1 = m2

# COORDINATES
# x1 = l*sin(theta1)                  y1 = -l*cos(theta1)
# x2 = l*sin(theta1) + l*sin(theta2)  y2 = -l*cos(theta1) - l*cos(theta2)

# Potential Energy
# V = m*g*y1 + m*g*y2 
# = -m*g*l*(2*cos(theta1) + cos(theta2))

# VELOCITIES
# dx1 = l * omega1 * cos(theta1)
# dy1 = l * omega1* sin(theta1)
# dx2 = l * omega1 * cos(theta1) + l * omega2 * cos(theta2)
# dy2 = l * omega1* sin(theta1) + l * omega2* sin(theta2)

# Kinetic Energy
# T = (1/2) * m * (dx1**2 + dy1**2) + (1/2) * m * (dx2**2 + dy2**2)
# T = (1/2) * m * (l**2) * (2 * detheta1**2 + omega2**2 + 2*omega1*omega2(cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2)))
# T = (1/2) * m * (l**2) * (2*(omega1**2) + omega2**2 + 2*omega1*omega2*cos(theta1-theta2))

# EQUATIONS OF MOTION
# domega1 = -(omega1**2 *sin(2*theta1 - 2*theta2) + 2 *omega2**2 *sin(theta1 - theta2) + (g/l)*(sin(theta1- 2*theta2) + 3*sin(theta1)))) / (3-cos(2*theta1 - 2*theta2))
# domega2 = (4*omega1**2 *sin(theta1 - theta2) + omega2**2 *sin(2*theta1 - 2*theta2) + (g/l)*(sin(2*theta1- theta2) - sin(theta2)))) / (3-cos(2*theta1 - 2*theta2))

def PartA(theta1, theta2, omega1, omega2, m, g, l):
    V = -m*g*l*(2*cos(theta1) + cos(theta2))
    T = (1/2)*m*(l**2) * (2*(omega1**2) + omega2**2 + 2*omega1*omega2*cos(theta1-theta2))
    E = V + T
    return E