import pylab
from math import cos, sin

# Standards
# l = l1 = l2
# m = m1 = m2

# COORDINATES
# x1 = l*sin(theta1)                  y1 = -l*cos(theta1)
# x2 = l*sin(theta1) + l*sin(theta2)  y2 = -l*cos(theta1) - l*cos(theta2)

# Potential Energy
# V = m*g*y1 + m*g*y2 = -m*g*l*(2*cos(theta1) + cos(theta2))

# VELOCITIES
# dx1 = l * dtheta1 * cos(theta1)
# dy1 = l * dtheta1* sin(theta1)
# dx2 = l * dtheta1 * cos(theta1) + l * dtheta2 * cos(theta2)
# dy2 = l * dtheta1* sin(theta1) + l * dtheta2* sin(theta2)

# Kinetic Energy
# T = (1/2) * m * (dx1**2 + dy1**2) + (1/2) * m * (dx2**2 + dy2**2)
# T = (1/2) * m * (l**2) * (2 * detheta1**2 + dtheta2**2 + 2*dtheta1*dtheta2(cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2)))
# T = (1/2) * m * (l**2) * (2*(dtheta1**2) + dtheta2**2 + 2*dtheta1*dtheta2*cos(theta1-theta2))

# EQUATIONS OF MOTION
# 