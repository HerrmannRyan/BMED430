#Import for python
import pandas
import numpy as np
import scipy.constants
import math
import matplotlib.pyplot as plt

#Import table
vol = 100 #liters
vol = vol*1000 #switch to mL
sigfigs = 4
Num_part = 22736
p_diam = 0.2 #cm
p_rad = p_diam/2
temp = 300 #K
vis = 0.852 #cP
molec_weight = 600 #g/mol
rhoA_s = 0.05 #g/cm^3
rho_drug = 1 #g/cc
rhoA_f = 0 #initial concentration
massA_f = rhoA_f*vol #g

#Timing for the loop for the interations and methods etc
timei = 0 #initial time
timef = 10 #min
timef = timef*60 #seconds
dt = 0.05 #seconds
snp_time = 30 #look at every thirty seconds

#Scipy constants
Kb = scipy.constants.Boltzmann #Boltzman constant
Av = scipy.constants.Avogadro

#conversion factors
#convert cp to pa*s
cP_Pas = 1e-3 

vis = vis*cP_Pas

p_surf = 4*np.pi*p_rad**2

#Definitions of variables given the input
def a_sol (MW):
    return ((3*MW)/(4*scipy.constants.pi*Av))**(1/3)

def Jaf_funct (K_m,Paf,Pas):
    return K_m*(Pas - Paf)

def Km_funct(D0, Diamp):
    return 2*D0/Diamp

def D0_funct(TempD,viscos,a):
    x = (Kb*TempD)/(6*scipy.constants.pi*viscos*a)
    return x

#using the function to make as
As = a_sol(molec_weight)/100

#moving that to D0
D0 = D0_funct(temp, vis, As)
D0 = D0*100**2

#Km calculation
Km = Km_funct(D0,p_diam)

#calculate for jaf
Jaf = Jaf_funct(Km,rhoA_f,rhoA_s)

#print(As, " ", D0, " ", Km, " ", Jaf) #print check

#initial table formatting
L_time = [timei]
L_timet = ['%.*g' % (sigfigs,timei)]
L_conc = [rhoA_f]
L_conct = ['%.*g' % (sigfigs,rhoA_f)]
L_mass = [massA_f]
L_masst = ['%.*g' % (sigfigs,massA_f)]

#set up loop
icount = 0
jcount = 0 #snapshot in time

while timei<timef:
    J_af = Jaf_funct(Km,rhoA_f,rhoA_s)
    drhoA = J_af*p_surf*Num_part/vol
    rhoA_n = rhoA_f + drhoA*dt

    rhoA_f = rhoA_n
    massA_f = rhoA_f * vol

    icount += 1
    jcount += 1
    timei = icount*dt
    stime = jcount*dt

    #snapshot
    if stime == snp_time:
        L_time.append(timei)
        L_timet.append('%.*g' % (sigfigs,timei))
        L_conc.append(rhoA_f)
        L_conct.append('%.*g' % (sigfigs,rhoA_f))
        L_mass.append(massA_f)
        L_masst.append('%.*g' % (sigfigs,massA_f))
        jcount = 0

print(f'concentration after ten minutes: {rhoA_f:.4g} g/cc')
print(f'concentration after ten minutes: {massA_f:.4g} g')

#create dictionary
results = {"Time (s)":L_time,"Drug Concentration (g/cc)": L_conct,"Mass Concentration (g)":L_masst}
df1 = pandas.DataFrame(results)
df1.set_index("Time (s)",inplace = True)
print(df1)
print("")

#make the theory graph
phi = Km*p_surf*Num_part/vol
xt = np.linspace(0,600,30)
y_massconc = rhoA_s*(1-np.exp(-phi*xt))*vol

fig = plt.figure(figsize = (10,8))
plt.plot(L_time,L_conc, color = "tab:red")
plt.title("Concentration Versus Time", fontsize = 20)
plt.xlabel("Time (s)",fontsize = 16.5)
plt.ylabel("Concentration (g/cc)", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
plt.show()
fig.savefig('Euler_Lab/Concentration.jpeg',dpi = 300,bbox_inches = 'tight')

fig = plt.figure(figsize = (10,8))
plt.plot(L_time,L_mass, color = "tab:red")
plt.title("Mass Concentration Versus Time", fontsize = 20)
plt.xlabel("Time (s)",fontsize = 16.5)
plt.ylabel("Mass Concentration (g)", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
plt.show()
fig.savefig('Euler_Lab/MassConcentration.jpeg',dpi = 300,bbox_inches = 'tight')

fig = plt.figure(figsize = (10,8))
plt.plot(xt,y_massconc, color = "tab:red")
plt.title("Mass Concentration Versus Time", fontsize = 20)
plt.xlabel("Time (s)",fontsize = 16.5)
plt.ylabel("Mass Concentration (g)", fontsize = 16.5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)
plt.show()
fig.savefig('Euler_Lab/MassConcentrationTheory.jpeg',dpi = 300,bbox_inches = 'tight')

L_dfs = [df1]
with open('Euler_Lab/ConcentrationTable.csv','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")