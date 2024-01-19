import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import pandas as pd

#
# Read in excel_file
#
excel_df = pd.read_excel('Lab2/viscosity_data.xlsx', sheet_name ='data_file')
#excel_df

tempD = excel_df.iloc[:,0].values
visd = excel_df.iloc[:,1].values

print(tempD)
print(visd)
#target temps
temps = [15.0,27.0,37.0,42.0,56.0,64.0] #deg C
L_exact = [1.139,0.852,0.692,0.629,0.496,0.440] #cP

sigfigs = 4

#linear interpolation
vis1_temp = np.interp(temps,tempD,visd)
vis1 = vis1_temp*1000.0
print(vis1)
print(L_exact)

#Spline
f=CubicSpline(tempD,visd,bc_type='natural')
vis2=f(temps)*1000.0
print(vis2)

#Curvefit
def fit_func(x,a,b):
    fit=a*np.log(x) + b
    return fit

#results from curvefitting
params = curve_fit(fit_func,tempD,visd)
a_fit = params[0][0]
b_fit = params[0][1]

#load lists
L_vis1 = [] #linear interpolation
L_vis2 = [] #cubic spine
L_vis3 = [] #curve fit
L_vis3r = [] #curve fit unformatted
L_temp = [] #temperature formatted

for i in temps:
    v1 = np.interp(i,tempD,visd)*1000.0
    v2 = f(i)*1000.0
    v3 = fit_func(i,a_fit,b_fit)*1000.0

    L_vis1.append('%.*g' % (sigfigs,v1))
    L_vis2.append('%.*g' % (sigfigs,v2))
    L_vis3.append('%.*g' % (sigfigs,v3))
    L_temp.append('%.*g' % (sigfigs,i))
    L_vis3r.append(v3) #for plotting curve-fit results

vis_dictionary = {'Temperature $^\\circ$%C':L_temp, 'Linear $cP':L_vis1, 'Cubic Spline $cP': L_vis2, 'Curve Fit $cP': L_vis3, 'Exact $cP':L_exact}
print(vis_dictionary)
df1 = pd.DataFrame(vis_dictionary)
df1.set_index('Temperature $^\\circ$%C',inplace=True)
#df1

xnew = np.linspace(5,85,30)
ynew = f(xnew)*1000.0

#New plot lines
fig = plt.figure(figsize=(10,8))
plt.subplot(3,1,1)
plt.rcParams.update({'font.size':14})
plt.plot(tempD,visd*1000,'bs')
plt.plot(temps,vis1,'ro')
plt.plot(xnew,ynew,'--k')
plt.title("Linear Interpolation",fontsize = 22)
plt.ylabel('Viscosity $cP$',fontsize = 18)
plt.xlabel('Temperature $^\\circ$C', fontsize = 18)
plt.grid(True)


plt.subplot(3,1,2)
plt.rcParams.update({'font.size':14})
plt.plot(tempD,visd*1000,'bs')
plt.plot(temps,vis2,'ro')
plt.plot(xnew,ynew,'--k')
plt.title("Cubic Fit",fontsize = 22)
plt.ylabel('Viscosity $cP$',fontsize = 18)
plt.xlabel('Temperature $^\\circ$C', fontsize = 18)
plt.grid(True)


plt.subplot(3,1,3)
plt.rcParams.update({'font.size':14})
plt.plot(tempD,visd*1000,'bs')
plt.plot(temps,L_vis3r,'ro')
plt.plot(xnew,ynew,'--k')
plt.title("Curve-Fit",fontsize = 22)
plt.ylabel('Viscosity $cP$',fontsize = 18)
plt.xlabel('Temperature $^\\circ$C', fontsize = 18)
plt.grid(True)

plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9, hspace = 1)

fig.savefig('Lab2/interpolationGraphs.jpeg',dpi = 300,bbox_inches = 'tight')

L_dfs = [excel_df,df1]
with open('Lab2/interpolationCSV.csv','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")