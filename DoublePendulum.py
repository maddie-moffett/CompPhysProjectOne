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

def PartA(theta1, theta2, omega1, omega2, m, l, g = 9.8): # PART A
    V = -m*g*l*(2*cos(theta1) + cos(theta2)) # potential energy (given)
    T = (1/2)*m*(l**2) * (2*(omega1**2) + omega2**2 + 2*omega1*omega2*cos(theta1-theta2)) # kinetic energy (given)
    E = V + T # total energy = potential + kinetic
    return E

def f(r, l, g = 9.8): # the function of calculation

    # unpack r
    theta1 = r[0]
    theta2 = r[1]
    omega1 = r[2]
    omega2 = r[3]

    dtheta1 = omega1 # equation one and two; derivative of theta is omega
    dtheta2 = omega2

    domega1topp1 = (omega1**2) * sin(2*theta1 - 2*theta2)                   # broke it all into pieces because it was a lot.
    domega1topp2 = 2 * (omega2**2) * sin(theta1 - theta2)
    domega1topp3 = (g/l) * (sin(theta1 - (2 * theta2)) + 3 * sin(theta1))
    domega1top = domega1topp1 + domega1topp2 + domega1topp3                 # numerator of domega 1 fraction
    domega1bot = 3 - cos(2*theta1 - 2*theta2)                               # denomenator of domega 1 fraction

    domega2topp1 = 4 * (omega1**2) * sin(theta1 - theta2)
    domega2topp2 = (omega2**2) * sin((2 * theta1) - (2 * theta2))
    domega2topp3 = 2 * (g/l) * (sin((2 * theta1) - theta2) - sin(theta2))
    domega2top = domega2topp1 + domega2topp2 + domega2topp3                  # numerator of domega 2 fraction
    domega2bot = 3 - cos((2 * theta1) - (2 * theta2))                        # denomenator of domega 2 fraction

    domega1 = -domega1top / domega1bot                                       # domega 1 = numerator / denomenator
    domega2 = domega2top / domega2bot                                        # domega 2 = numerator / denomenator

    return [dtheta1, dtheta2, domega1, domega2] # return the derivatives in the same order received

def multi(h, r):               # helper function for scalar multiplication of a vector
    j = copy(r)                # copy the original vector
    for num in range(len(j)):  # iterate through
        j[num] *= h            # multiply by the scalar
    return j                   # return

def addi(lone, ltwo):          # helper function to add two vectors pairwise
    retti = []                 # initialize empty list
    for i in range(len(lone)): # iterate through and append to the empty list
        retti.append(lone[i] + ltwo[i]) #  the sum of the values in the two vectors at that position
    return retti               # return

def rungkut(r, h, l, g): # rungekutta one step function;
                         # essentially given but uses helper functions for multiplication and division of the arrays
    kone = multi(h, f(r, l, g))
    ktwo = multi(h, f(addi(r, multi(0.5, kone)), l, g))
    kthree = multi(h, f(addi(r, multi(0.5, ktwo)), l, g))
    kfour = multi(h, f(addi(r, kthree), l, g))
    k = addi(r, multi((1/6), addi(addi(kone, multi(2, ktwo)), addi(multi(2, kthree), kfour))))
    return k

def calculate(theta1_og, theta2_og, omega1_og, omega2_og, m, l, g, tmin, tmax, h):
    theta1s = [theta1_og] # initialize list of values with given as the first entry
    theta2s = [theta2_og] # initialize list of values with given as the first entry
    omega1s = [omega1_og] # initialize list of values with given as the first entry
    omega2s = [omega2_og] # initialize list of values with given as the first entry
    ts = [tmin]           # initialize list of values with given as the first entry
    t = tmin              # start iterator value at the minimum t given

    while t <= tmax:              # while loop run from tmin to tmax
        theta1_curr = theta1s[-1] # grab the current positions from the arrays
        theta2_curr = theta2s[-1] # grab the current positions from the arrays
        omega1_curr = omega1s[-1] # grab the current positions from the arrays
        omega2_curr = omega2s[-1] # grab the current positions from the arrays
        r_old = [theta1_curr, theta2_curr, omega1_curr, omega2_curr] # assign current vals to own array
        r_new = rungkut(r_old, h, l, g) # calculate step forward using rungekutta function
        theta1s.append(r_new[0])  # append new values to lists of values
        theta2s.append(r_new[1])  # append new values to lists of values
        omega1s.append(r_new[2])  # append new values to lists of values
        omega2s.append(r_new[3])  # append new values to lists of values
        ts.append(t)              # store current time
        t += h                    # step forward in time
    
    return theta1s, theta2s, omega1s, omega2s, ts # return all pertinent information

def totalEnergy(theta1s, theta2s, omega1s, omega2s, ts, m, l, g, labelop = ""):   # calculate and graph total energy over time; option for label
    Es =[]                                                                        # initialize empty array for Energy vals
    new_ts = []                                                                   # initialize empty array for t vals
                                                                                  # shrink bc rungekutta gives too many t vals to graph
    for i in range(0, len(theta1s), 10):                                          # iterate through og t vals, only every 10 (^^)
        Es.append(PartA(theta1s[i], theta2s[i], omega1s[i], omega2s[i], m, l, g)) # use part A to calculate E, append to E array
        new_ts.append(ts[i])                                                      # append the new t val as well
    pylab.plot(new_ts, Es, label = labelop)                                       # plot all E vals over time
    pylab.xlabel("Time in seconds")                                               # label x axis
    pylab.ylabel("Total Energy Change of the System in Joules")                   # label y axis
    pylab.title("Conservation of Energy")                                         # label graph
    return Es                                                                     # return Es

def getX1(l, theta1):            # calculate x1 position for given l and theta
    return l*sin(theta1)
def getY1(l, theta1):            # calculate y1 position for given l and theta
    return -l*cos(theta1)
def getX2(l, theta1, theta2):    # calculate x2 position for given l and thetas
    return l*sin(theta1) + l*sin(theta2)
def getY2(l, theta1, theta2):    # calculate y2 position for given l and thetas
    return -l*cos(theta1) - l*cos(theta2)

def GraphXs(theta1s, theta2s, ts, l, labelop1 = "First Mass", labelop2= "Second Mass"): # calculate and graph x positions over time; option to change label
    x1s = []    # initialize empty array for x1 vals
    x2s = []    # initialize empty array for x2 vals
    new_ts = [] # initialize empty array for new t vals
    
    for i in range(0, len(theta1s), 10):             # iterate through one every ten ts
        x1s.append(getX1(l, theta1s[i]))             # append the x1 pos to its array at that t
        x2s.append(getX2(l, theta1s[i], theta2s[i])) # append the x2 pos to its array at that t
        new_ts.append(ts[i])                         # append the t val to its array

    pylab.plot(new_ts, x1s, label = labelop1)        # plot all of the x1s over time
    pylab.plot(new_ts, x2s, label = labelop2)        # plot all of the x2s over time
    pylab.xlabel("Time in seconds")                  # label the x axis
    pylab.ylabel("X-position")                       # label the y axis
    pylab.title("X-position Over Time")              # title the graph
    pylab.legend()                                   # create the legend for the graph
    return x1s, x2s                                  # return the calculated values

def GraphYs(theta1s, theta2s, ts, l, labelop1 = "First Mass", labelop2= "Second Mass"): # calculate and graph y positions over time; option to change label
    y1s = []    # initialize empty array for y1 vals
    y2s = []    # initialize empty array for y2 vals
    new_ts = [] # initialize empty array for new t vals
    
    for i in range(0, len(theta1s), 10):             # iterate through one every ten ts
        y1s.append(getY1(l, theta1s[i]))             # append the y1 pos to its array at that t
        y2s.append(getY2(l, theta1s[i], theta2s[i])) # append the y2 pos to its array at that t
        new_ts.append(ts[i])                         # append the t val to its array

    pylab.plot(new_ts, y1s, label = labelop1)        # plot all of the y1s over time
    pylab.plot(new_ts, y2s, label = labelop2)        # plot all of the y2s over time
    pylab.xlabel("Time in seconds")                  # label the x axis
    pylab.ylabel("Y-position")                       # label the y axis
    pylab.title("Y-position Over Time")              # title the graph
    pylab.legend()                                   # create the legend for the graph
    return y1s, y2s                                  # return the calculated values

def PartB(tmin, tmax, h, theta1, theta2, omega1, omega2, m, l, g = 9.8): # function for part b
    # use the calculate function to collect the thetas and omegas over time from given min to given max t
    theta1s, theta2s, omega1s, omega2s, ts = calculate(theta1, theta2, omega1, omega2, m, l, g, tmin, tmax, h)
    totalEnergy(theta1s, theta2s, omega1s, omega2s, ts, m, l, g) # calculate the total energy over t, graph it
    pylab.show()                                                 # show above graph
    GraphXs(theta1s, theta2s, ts, l)                             # calculate the x-positions over t, graph it
    pylab.show()                                                 # show above graph
    GraphYs(theta1s, theta2s, ts, l)                             # calculate the y-positions over t, graph it
    pylab.show()                                                 # show above graph

def PartC():       # part C
    tmin = 0       # minimum t = 0
    tmax = 100     # maximum t = 100
    h = 0.001      # h = 0.001
    l = 0.4        # length1 = length2 = 0.4
    m = 1          # NOT EXPLICITLY GIVEN; assumed from the book
    theta1 = pi/2  # theta1 = theta2 = pi/2
    theta2 = pi/2  # theta1 = theta2 = pi/2
    omega1 = 0     # NOT EXPLICITLY GIVEN; assumed from part D
    omega2 = 0     # NOT EXPLICITLY GIVEN; assumed from part D
    PartB(tmin, tmax, h, theta1, theta2, omega1, omega2, m, l)

 # helper function iterates through angles, makes the x and y graphs, and saves them to a folder within the git directory
def PartDHelper(theta1, theta2, theta3, theta4, theta5, omega1, omega2, m, l, g, tmin, tmax, h):
    thetaops = ["pi/6", "pi/4", "pi/3", "pi/2", "pi"]      # theta options spelled out for legend and naming purposes
    location1 = "FiguresForPartD/"                         # location within the git directory
    location2ops = ["pio6", "pio4", "pio3", "pio2", "pi"]  # location 2 options; theta options spelled out for file saving
    thetas = [theta1, theta2, theta3, theta4, theta5]      # theta options in usable form
    for A in range(len(thetas)):                           # iterate through thetas for theta 1
        for B in range(len(thetas)):                       # iterate through thetas for theta 2
            # calculate to get the thetas over time
            thetaAs, thetaBs, omega1s, omega2s, ts = calculate(thetas[A], thetas[B], omega1, omega2, m, l, g, tmin, tmax, h)
            # graph the xs, give the alternative label option to make it easier to see what graph it is
            GraphXs(thetaAs, thetaBs, ts, l, labelop1 = "First x of " + thetaops[A] + " and " + thetaops[B], labelop2 = "Second x of " + thetaops[A] + " and " + thetaops[B])
            pylab.legend()                                                      # show the legend
            pylab.savefig(location1 + "X" + location2ops[A] + location2ops[B])  # save the figure at the desired location
            pylab.clf()                                                         # clear the figure
            # graph the ys, give the alternative label option to make it easier to see what graph it is
            GraphYs(thetaAs, thetaBs, ts, l, labelop1 = "First y of " + thetaops[A] + " and " + thetaops[B], labelop2 = "Second y of " + thetaops[A] + " and " + thetaops[B])
            pylab.legend()                                                      # show the legend
            pylab.savefig(location1 + "Y" + location2ops[A] + location2ops[B])  # save the figure at the desired location
            pylab.clf()                                                         # clear the figure
    

def PartD():       # Part D
    tmin = 0       # minimum t = 0
    tmax = 100     # maximum t = 100
    h = 0.001      # h = 0.001 (accuracy from part C)
    l = 0.4        # length1 = length2 = 0.4
    m = 1          # NOT GIVEN; assumed from the book
    g = 9.8        # known constant
    theta1 = pi/6  # theta1 option one is pi/6
    theta2 = pi/4  # theta1 option two is pi/4
    theta3 = pi/3  # theta1 option three is pi/3
    theta4 = pi/2  # theta1 option four is pi/2
    theta5 = pi    # theta1 option five is pi
    omega1 = 0     # GIVEN omega1 = omega2 = 0
    omega2 = 0     # GIVEN omega1 = omega2 = 0
    PartDHelper(theta1, theta2, theta3, theta4, theta5, omega1, omega2, m, l, g, tmin, tmax, h) # use above in partD helper



def roundem(r):                       # round vals in a vector HELPER function
    t = []                            # initialize empty vector for new vals
    for val in r:                     # iterate through vals in old vector
        t.append(int(round(val, 0)))  # append to new vector the rounded version of the one in the old vector
    return t                          # return

# FUNction to make all of our fancy animations!!
def animate(location2, theta1 = pi/2, theta2 = pi/2, omega1 = 0, omega2 = 0, m = 1, l = 0.4, g = 9.8, tmin = 0, tmax = 100, h = 0.01):
    # calculate the values over time
    theta1s, theta2s, omega1s, omega2s, ts = calculate(theta1, theta2, omega1, omega2, m, l, g, tmin, tmax, h)
    ts = multi(100, ts) # multiply the t values by 100 so they're all unique
    ts = roundem(ts)    # round the t values so they can be used

    def init():                         # nested function: initiate the animation
        ax.plot([0], [0], "ko", ms = 5) # plot the anchor point

    def update(frame):                  # nested function: what happens each frame
        ax.clear()                      # clear old frame
        x1 = getX1(l, theta1s[frame])   # calculate x val
        y1 = getY1(l, theta1s[frame])   # calculate y val
        x2 = getX2(l, theta1s[frame], theta2s[frame])   # calculate second x val
        y2 = getY2(l, theta1s[frame], theta2s[frame])   # calculate second y val

        ax.plot([0, x1, x2], [0, y1, y2], "k-", lw = 2) # plot the lines between our two masses (the rods holding them together)
        ax.plot([0], [0], "ko", lw = 5)                 # plot the anchor point as a black dot
        ax.plot([x1, x2], [y1, y2], "ro", ms = 10)      # plot the two masses as red circles
        ax.set_xlim(-0.9, 0.9)                          # domain of x vals
        ax.set_ylim(-0.9, 0.6)                          # range of y vals (note this was actually +1.5 ONLY for the pi/pi animation)
        ax.axis("off")                                  # do not show the axes

    location1 = "Animations/"                           # saving to the animations folder

    fig, ax = pylab.subplots()                          # create the plot background
    fig.set_facecolor("lightgray")                      # set background color
    ax.set_xlim(-0.9, 0.9)                              # domain of x vals
    ax.set_ylim(-0.9, 0.6)                              # range of y vals (note this was actually +1.5 ONLY for the pi/pi animation)
    ax.axis("off")                                      # don't show the axes
    ani = FuncAnimation(fig, update, init_func = init, frames = ts, interval = 10, blit = False, repeat = False) # ANIMATE
    ani.save(location1 + "ANIMATION" + location2 + ".gif")                                                       # save the animation at location
    return ani                                          # return just in casekies

def AniHelper(theta1, theta2, theta3, theta4, theta5):                        # helper function for animate
    location2ops = ["pio6", "pio4", "pio3", "pio2", "pi"]                     # provide the options for thetas for the save file name
    thetas = [theta1, theta2, theta3, theta4, theta5]                         # provide the options for the thetas for calculations
    for A in range(len(thetas)):                                              # iterate through theta1s
        for B in range(len(thetas)):                                          # iterate through theta2s
            animate(location2ops[A] + location2ops[B], thetas[A], thetas[B])  # animate that pair of thetas, save at the location w/ their names