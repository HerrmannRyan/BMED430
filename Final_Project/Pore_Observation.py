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
we will be simplifying this to having an inital density over the total fiber

The values should be run from time 0 which is given at day four to day 28 -- we declare the time step of 1 hour,
this means that it is 28-4 days 24 days time 24 hours gives us the number of steps

We also need to probably deal with different parts of pore, but for now we only are looking that the conditions
along the pore boundary (this may be a terrible assumption to make)"""

#Define constants
fiberDiameter = 1.70 #microns

#using endothelial cells
prolifRate = .5268 #1/day
maxDensity = 4e-3 #cells/um <-- this is our K defined in the def direchBC
dt = 1 #hour
dt = dt/24 #days


#initial conditions
sigfigs = 4
timestart = 0
timeend = 24 #days assuming that time 0 is technically day 4
daycount = []

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
densityDayList = []
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

    if(timecount%24 == 0):
        densityDayList.append(density)
        daycount.append(timecount/24)

#takes values in the density list and formats them to a new list with correct significant figures
for i in range(len(densityDayList)):
    densityListFormatted.append('%.*g' % (sigfigs,densityDayList[i]))

#this needs to be formatted
fig = plt.figure(figsize = (10,8))
plt.plot(timeList,densityList, color = "tab:red")
plt.title("Cell Density in the Pore vs. Time", fontsize = 20)
plt.xlabel("Time (Days)",fontsize = 16.5)
plt.ylabel("Density (cell/um^2)", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
#plt.show()

#save figure
fig.savefig('Final_Project/densityGraph.jpeg',dpi = 300,bbox_inches = 'tight')

"""The coverage function takes use of the area of the void and then
uses that tpo calculate the coverage of the cells. There needs to be a 'for
loop' that will look each value of u(t) and then calculate coverage

For this equation we are using L which is the pore size length and that
is set from the experiment. For this we are using 400 microns"""

#Coverage values need to be specified
length = 84.1 #microns

#A_void function
def a_void(ut):
    diam = fiberDiameter
    poreLength = length

    poreArea = poreLength**2
    pi = np.pi
    num_cells = ut*poreArea
    r = diam/2

    avoid = poreArea - ((pi*(r**2))*num_cells)

    return avoid

#coverage function
def coverageFunct(ut):
    avoid = a_void(ut) #calculate a_void
    cov = 1 - (avoid/length**2)
    
    return cov

#create list for the coverage
coverageList = []
coverageListFormatted = []

#get coverage for each density
for i in range(len(densityList)):
    coverage = coverageFunct(densityList[i])

    #normalizing the function for the coverage and then multiplying by 100 for %
    coverage = coverage/.01*100

    #add to list 
    coverageList.append(coverage)

#format to list with correct sigfigs using the previous list 
for i in range(len(densityList)):
    if(i % 24 == 0):
        formattedValue = '%.*g%%' % (sigfigs,coverageList[i])
        coverageListFormatted.append(formattedValue)

fig2 = plt.figure(figsize = (10,8))
plt.plot(timeList,coverageList, color = "tab:green")
plt.title("Cell Coverage vs Time", fontsize = 20)
plt.xlabel("Time (Days)",fontsize = 16.5)
plt.ylabel("Coverage %", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
#plt.show()

#save fig
fig.savefig('Final_Project/coverageGraph.jpeg',dpi = 300,bbox_inches = 'tight')


#create a dataframe for the coverage and density
results = {"Time (day)":daycount,"Cell Density (cells/um^2)": densityListFormatted,"Coverage":coverageListFormatted}
df1 = pandas.DataFrame(results)
df1.set_index("Time (day)",inplace = True)
print(df1)
print("")

#save files
L_dfs = [df1]
with open('Final_Project/DensityAndCoverage.csv','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")