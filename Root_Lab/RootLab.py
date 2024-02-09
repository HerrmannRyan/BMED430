import numpy as np
import pandas
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#import variables
a = 100 #nm
sigmar = 0.2 #fraction reflected
epi = 1e-5 #convergence criteria
error = 10.0 #initialize the error
sigfigs = 4
guess = 0.9 #Lambda = a/rpore
x0 = guess
rpore = a/guess
icount = 0 #initial count

#define lists for del^2
L_error = ['%.*g' % (sigfigs,error)]
L_lam = ['%.*g' % (sigfigs,guess)]
L_rpore = ['%.*g' % (sigfigs,rpore)]
L_icount = [icount]

#list for fsolve
L_fsolve = []

def g(x,sigr):
    kpore = (1-x)**2
    term2 = 2-kpore
    term3 = 1- (2*x**2)/3 - 0.163*x**3
    term4 = kpore*term2*term3
    f = sigr-(1-term4)
    return f

#test = g(x0,sigmar)
lam1 = fsolve(g,x0,sigmar)
L_fsolve.append('%.*g' % (sigfigs,lam1[0]))
print(L_fsolve)

#del squared functions
def g1(x,sigr):
    kpore = (1-x)**2
    term2 = 2-kpore
    term3 = 1 - (2*x**2)/3 - 0.163*x**3
    f1 = 1.0 - math.sqrt((1-sigr)/(term2*term3))
    return f1

def gprime(y0,y1,y2):
    f2 = (y2-y1)/(y1-y0)
    return f2

def g3(y1,y2,gdev):
    f3 = y2 + (gdev/(1-gdev))*(y2-y1)
    return f3

while error >= epi:
    x1 = g1(x0,sigmar)
    x2 = g1(x1,sigmar)
    derivative = gprime(x0,x1,x2)
    
    xdel2 = g3(x1,x2,derivative)

    error = np.abs(xdel2 -x0)/xdel2
    x0 = xdel2
    icount = icount + 1
    L_error.append('%.*g' % (sigfigs,error))
    L_icount.append(icount)
    L_lam.append('%.*g' % (sigfigs,x0))
    L_rpore.append('%.*g' % (sigfigs,rpore))

#creation of dictionary
results = {"Iterations":L_icount,"Lambda":L_lam, "Pore Radius (nm)":L_rpore, "Relative Error":L_error}

df1 = pandas.DataFrame(results)
df1.set_index("Iterations",inplace= True)
print(df1)
print("")

given = {"Sigmar $\\sigma_r$":sigmar,"Fsolve":L_fsolve, "Del Squared":L_lam[-1],"Initial Guess":L_lam[0],"Iterations":L_icount[-1],"Rel Error":L_error[-1]}
df2 = pandas.DataFrame(given)
print(df2)

#plotting the data
xnew = np.linspace(0.01,0.99,30)
ynew = g(xnew,sigmar)

fig = plt.figure(figsize = (10,8))
plt.plot(xnew,ynew, color = "tab:orange")
plt.plot(xnew,np.zeros(len(xnew)), "--", color = "tab:gray")
plt.plot(x0,g(x0,sigmar),marker = "p",color = "tab:purple", markersize = 8)
plt.title("Sigmar-f($\\lambda$) vs. $\\lambda$", fontsize = 20)
plt.xlabel("$\\lambda$",fontsize = 16.5)
plt.ylabel("Sigmar-f($\\lambda$)", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
plt.show()

fig.savefig('delSquaredGraph.jpeg',dpi = 300,bbox_inches = 'tight')

L_dfs = [df2,df1]
with open('DelsquaredCSV.csv','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")