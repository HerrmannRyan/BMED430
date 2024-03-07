import pandas
import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import odeint # <--this may not have to be used since we are doing euler

"""For the following process,
We are using a simpified equation given assumptions
u is the cell density and we will be solving it through an euler method
needs a start and end condition

We need to note the initial conditions:
u(x,t0) is 0 when x = omega where omega is the center of the pore
u(x,t0) is u0 <-- this value needs to be assumed and that is at x=domega and that is the boundary of the pore

The values should be run from time 0 which is given at day four to day 28 -- we declare the time step of 1 hour,
this means that it is 28-4 days 24 days time 24 hours gives us the number of steps

We also need to probably deal with different parts of pore, but for now we only are looking that the conditions
along the pore boundary (this may be a terrible assumption to make)"""

#Define constants
fiberDiameter = 1.94 #microns
fiberDiameterSecond = 1.74 #microns

#using endothelial cells
prolifRate = .5268 #1/day
maxDensity = 4e-3 #cells/um <-- this is our K defined in the def direchBC
dt = 1 #hour
dt = dt/24 #days


#initial conditions
sigfigs = 4
timestart = 0
timeend = 24 #days assuming that time 0 is technically day 4

initDensity = 2e-3 #cells/um^2 <-- from the paper can be changed. This is u0
#initDensity = 0 #cells/um^2 <-- this is at x = 0

#pore diam --> equation with time step --> with both fiber diameters

#governing equation to be iterated by dt
def direchBC (u):
    lamb = prolifRate
    K = maxDensity
    dudt = lamb*u*(1-(u/K))
    return dudt

#declare lists to be tracked
densityList = []
densityListFormatted = []
timeList = []

#set up and start loop
density = initDensity
timecount = 0

while(timestart<timeend):
    u = direchBC(density)
    newDensity = density + u*dt

    densityList.append(density)
    density = newDensity
    
    timecount += 1
    timestart = timecount*dt
    timeList.append(timestart)

#takes values in the density list and formats them to a new list with correct significant figures
for i in range(len(densityList)):
    densityListFormatted.append('%.*g' % (sigfigs,densityList[i]))

#print(densityListFormatted)

#this needs to be formatted
fig = plt.figure(figsize = (10,8))
plt.plot(timeList,densityList, color = "tab:red")
plt.grid(True)
plt.show()

#need the coverage function
tau = 0.5 #<-- this was fixed in the paper

#A_void function needed here