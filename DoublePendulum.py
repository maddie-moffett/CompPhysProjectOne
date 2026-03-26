import pylab
from math import cos, sin
from copy import copy

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

def PartA(theta1, theta2, omega1, omega2, m, l, g = 9.8):
    V = -m*g*l*(2*cos(theta1) + cos(theta2))
    T = (1/2)*m*(l**2) * (2*(omega1**2) + omega2**2 + 2*omega1*omega2*cos(theta1-theta2))
    E = V + T
    return E

def f(r, m, l, g = 9.8):

    theta1 = r[0]
    theta2 = r[1]
    omega1 = r[2]
    omega2 = r[3]

    dtheta1 = omega1
    dtheta2 = omega2

    domega1 = -(omega1**2 *sin(2*theta1 - 2*theta2) + 2 *omega2**2 *sin(theta1 - theta2) + (g/l)*(sin(theta1- 2*theta2) + 3*sin(theta1))) / (3-cos(2*theta1 - 2*theta2))
    domega2 = (4*omega1**2 *sin(theta1 - theta2) + omega2**2 *sin(2*theta1 - 2*theta2) + (g/l)*(sin(2*theta1- theta2) - sin(theta2))) / (3-cos(2*theta1 - 2*theta2))

    return [dtheta1, dtheta2, domega1, domega2]

def multi(h, r):
    j = copy(r)
    for num in range(len(j)):
        j[num] *= h
    return j

def addi(lone, ltwo):
    retti = []
    for i in range(len(lone)):
        retti.append(lone[i] + ltwo[i])
    return retti

def rungkut(r, t, h, m, l, g):
    kone = multi(h, f(r, t, m, l, g))
    ktwo = multi(h, f(addi(r, multi(0.5, kone)), t + 0.5*h, m, l, g))
    kthree = multi(h, f(addi(r, multi(0.5, ktwo)), t + 0.5*h, m, l, g))
    kfour = multi(h, f(addi(r, kthree), t + h, m, l, g))
    k = addi(r, multi((1/6), addi(addi(kone, multi(2, ktwo)), addi(multi(2, kthree), kfour))))
    return k

def calculate(theta1_og, theta2_og, omega1_og, omega2_og, m, l, g, tmin, tmax, h):
    theta1s = [theta1_og]
    theta2s = [theta2_og]
    omega1s = [omega1_og]
    omega2s = [omega2_og]
    ts = [tmin]
    t = tmin

    while t <= tmax:
        theta1_curr = theta1s[-1]
        theta2_curr = theta2s[-1]
        omega1_curr = omega1s[-1]
        omega2_curr = omega2s[-1]
        r_old = [theta1_curr, theta2_curr, omega1_curr, omega2_curr]
        r_new = rungkut(r_old, t, h, m, l, g)
        theta1s.append(r_new[0])
        theta2s.append(r_new[1])
        omega1s.append(r_new[2])
        omega2s.append(r_new[3])
        ts.append(t)
        t += h