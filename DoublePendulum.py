import pylab
from math import cos, sin, pi
from copy import copy
from matplotlib.animation import FuncAnimation

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

def f(r, l, g = 9.8):

    theta1 = r[0]
    theta2 = r[1]
    omega1 = r[2]
    omega2 = r[3]

    dtheta1 = omega1
    dtheta2 = omega2

    domega1topp1 = (omega1**2) * sin(2*theta1 - 2*theta2)
    domega1topp2 = 2 * (omega2**2) * sin(theta1 - theta2)
    domega1topp3 = (g/l) * (sin(theta1 - (2 * theta2)) + 3 * sin(theta1))
    domega1top = domega1topp1 + domega1topp2 + domega1topp3
    domega1bot = 3 - cos(2*theta1 - 2*theta2)

    domega2topp1 = 4 * (omega1**2) * sin(theta1 - theta2)
    domega2topp2 = (omega2**2) * sin((2 * theta1) - (2 * theta2))
    domega2topp3 = 2 * (g/l) * (sin((2 * theta1) - theta2) - sin(theta2))
    domega2top = domega2topp1 + domega2topp2 + domega2topp3
    domega2bot = 3 - cos((2 * theta1) - (2 * theta2))

    domega1 = -domega1top / domega1bot
    domega2 = domega2top / domega2bot

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

def rungkut(r, h, l, g):
    kone = multi(h, f(r, l, g))
    ktwo = multi(h, f(addi(r, multi(0.5, kone)), l, g))
    kthree = multi(h, f(addi(r, multi(0.5, ktwo)), l, g))
    kfour = multi(h, f(addi(r, kthree), l, g))
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
        r_new = rungkut(r_old, h, l, g)
        theta1s.append(r_new[0])
        theta2s.append(r_new[1])
        omega1s.append(r_new[2])
        omega2s.append(r_new[3])
        ts.append(t)
        t += h
    
    return theta1s, theta2s, omega1s, omega2s, ts

def totalEnergy(theta1s, theta2s, omega1s, omega2s, ts, m, l, g, labelop = ""):
    Es =[]
    new_ts = []
    for i in range(0, len(theta1s), 10):
        Es.append(PartA(theta1s[i], theta2s[i], omega1s[i], omega2s[i], m, l, g))
        new_ts.append(ts[i])
    pylab.plot(new_ts, Es, label = labelop)
    pylab.xlabel("Time in seconds")
    pylab.ylabel("Total Energy Change of the System in Joules")
    pylab.title("Conservation of Energy")
    return Es

def getX1(l, theta1):
    return l*sin(theta1)
def getY1(l, theta1):
    return -l*cos(theta1)
def getX2(l, theta1, theta2):
    return l*sin(theta1) + l*sin(theta2)
def getY2(l, theta1, theta2):
    return -l*cos(theta1) - l*cos(theta2)

def GraphXs(theta1s, theta2s, ts, l, labelop1 = "First Mass", labelop2= "Second Mass"):
    x1s = []
    x2s = []
    new_ts = []
    
    for i in range(0, len(theta1s), 10):
        x1s.append(getX1(l, theta1s[i]))
        x2s.append(getX2(l, theta1s[i], theta2s[i]))
        new_ts.append(ts[i])

    pylab.plot(new_ts, x1s, label = labelop1)
    pylab.plot(new_ts, x2s, label = labelop2)
    pylab.xlabel("Time in seconds")
    pylab.ylabel("X-position")
    pylab.title("X-position Over Time")
    pylab.legend()
    return x1s, x2s

def GraphYs(theta1s, theta2s, ts, l, labelop1 = "First Mass", labelop2= "Second Mass"):
    y1s = []
    y2s = []
    new_ts = []
    
    for i in range(0, len(theta1s), 10):
        y1s.append(getY1(l, theta1s[i]))
        y2s.append(getY2(l, theta1s[i], theta2s[i]))
        new_ts.append(ts[i])

    pylab.plot(new_ts, y1s, label = labelop1)
    pylab.plot(new_ts, y2s, label = labelop2)
    pylab.xlabel("Time in seconds")
    pylab.ylabel("Y-position")
    pylab.title("Y-position Over Time")
    pylab.legend()
    return y1s, y2s

def PartB(tmin, tmax, h, theta1, theta2, omega1, omega2, m, l, g = 9.8):
    theta1s, theta2s, omega1s, omega2s, ts = calculate(theta1, theta2, omega1, omega2, m, l, g, tmin, tmax, h)
    totalEnergy(theta1s, theta2s, omega1s, omega2s, ts, m, l, g)
    pylab.show()
    GraphXs(theta1s, theta2s, ts, l)
    pylab.show()
    GraphYs(theta1s, theta2s, ts, l)
    pylab.show()

def PartC():
    tmin = 0
    tmax = 100
    h = 0.002
    l = 0.4
    m = 1          # NOT GIVEN
    theta1 = pi/2
    theta2 = pi/2
    omega1 = 0     # NOT GIVEN
    omega2 = 0     # NOT GIVEN
    PartB(tmin, tmax, h, theta1, theta2, omega1, omega2, m, l)
    # h = 0.001 energy loss 3 * 10^(-7)

def PartDHelper(theta1, theta2, theta3, theta4, theta5, omega1, omega2, m, l, g, tmin, tmax, h):
    thetaops = ["pi/6", "pi/4", "pi/3", "pi/2", "pi"]
    location1 = "C:/Users/maddi/Desktop/CompPhys/ProjectOne/FiguresForPartD/"
    location2ops = ["pio6", "pio4", "pio3", "pio2", "pi"]
    thetas = [theta1, theta2, theta3, theta4, theta5]
    for A in range(len(thetas)):
        for B in range(len(thetas)):
            thetaAs, thetaBs, omega1s, omega2s, ts = calculate(thetas[A], thetas[B], omega1, omega2, m, l, g, tmin, tmax, h)
            GraphXs(thetaAs, thetaBs, ts, l, labelop1 = "First x of " + thetaops[A] + " and " + thetaops[B], labelop2 = "Second x of " + thetaops[A] + " and " + thetaops[B])
            pylab.legend()
            pylab.savefig(location1 + "X" + location2ops[A] + location2ops[B])
            pylab.clf()
            GraphYs(thetaAs, thetaBs, ts, l, labelop1 = "First y of " + thetaops[A] + " and " + thetaops[B], labelop2 = "Second y of " + thetaops[A] + " and " + thetaops[B])
            pylab.legend()
            pylab.savefig(location1 + "Y" + location2ops[A] + location2ops[B])
            pylab.clf()
    

def PartD():
    tmin = 0
    tmax = 100
    h = 0.001
    l = 0.4
    m = 1          # NOT GIVEN
    g = 9.8
    theta1 = pi/6
    theta2 = pi/4
    theta3 = pi/3
    theta4 = pi/2
    theta5 = pi
    omega1 = 0     # NOT GIVEN
    omega2 = 0     # NOT GIVEN
    PartDHelper(theta1, theta2, theta3, theta4, theta5, omega1, omega2, m, l, g, tmin, tmax, h)



def roundem(r):
    t = []
    for val in r:
        t.append(int(round(val, 0)))
    return t

def animate(location2, theta1 = pi/2, theta2 = pi/2, omega1 = 0, omega2 = 0, m = 1, l = 0.4, g = 9.8, tmin = 0, tmax = 100, h = 0.01):
    theta1s, theta2s, omega1s, omega2s, ts = calculate(theta1, theta2, omega1, omega2, m, l, g, tmin, tmax, h)
    ts = multi(100, ts)
    ts = roundem(ts)

    def init():
        ax.plot([0], [0], "ko", ms = 5)

    def update(frame):
        ax.clear()
        x1 = getX1(l, theta1s[frame])
        y1 = getY1(l, theta1s[frame])
        x2 = getX2(l, theta1s[frame], theta2s[frame])
        y2 = getY2(l, theta1s[frame], theta2s[frame])

        ax.plot([0, x1, x2], [0, y1, y2], "k-", lw = 2)
        ax.plot([0], [0], "ko", lw = 5)
        ax.plot([x1, x2], [y1, y2], "ro", ms = 10)
        ax.set_xlim(-0.9, 0.9)
        ax.set_ylim(-0.9, 1.5)
        ax.axis("off")

    location1 = "Animations/"

    fig, ax = pylab.subplots()
    fig.set_facecolor("lightgray")
    ax.set_xlim(-0.9, 0.9)
    ax.set_ylim(-0.9, 1.5)
    ax.axis("off")
    ani = FuncAnimation(fig, update, init_func = init, frames = ts, interval = 10, blit = False, repeat = False)
    ani.save(location1 + "ANIMATION" + location2 + ".gif")
    return ani

def AniHelper(theta1, theta2, theta3, theta4, theta5):
    location2ops = ["pio6", "pio4", "pio3", "pio2", "pi"]
    thetas = [theta1, theta2, theta3, theta4, theta5]
    for A in range(len(thetas)):
        for B in range(len(thetas)):
            animate(location2ops[A] + location2ops[B], thetas[A], thetas[B])

if __name__ == "__main__":
    PartC()