#Import for python
import pandas
import numpy as np
import scipy.constants
import math
import matplotlib.pyplot as plt

#Basic Formatting
sigfigs = 4
N_samples = 100

#Import table
vol = 100 #liters
Num_part = 22736
p_diam = 0.2 #cm
p_rad = p_diam/2
temp = 300 #K
vis = 0.852 #cP
molec_weight = 600 #g/mol
rho_drug = 1 #g/cc
rhoA_f = 0 #initial concentration
massA_f = rhoA_f*vol #g
unc = 0.05

#Timing for the loop for the interations and methods etc
timei = 0 #initial time
timef = 10 #min
timef = timef*60 #seconds
dt = 0.05 #seconds

#Scipy constants
Kb = scipy.constants.Boltzmann #Boltzman constant
Av = scipy.constants.Avogadro

#Separate Constants

#conversion factors
#convert cp to pa*s
cP_Pas = 1e-3
vis = vis*cP_Pas  
#print(rhoAs_a)

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

vol = vol*1000 #switch to mL
rhoA_s = 0.05 #g/cm^3

#lists for results
stDev_List = []
stDev_List_M= []
concMean_list = []
concMean_listM = []

stDev_ListF = []
stDev_List_MF= []
concMean_listF = []
concMean_listMF = []

#lists for number of samples 
num_samples = []

#lists of +/- rho and vol these will be formatted
rhoA_uncertainty = []
volf_uncertainty = []

#all values for the plot with number of samples as the first column
whole_list = []
plot_samples = np.arange(0,N_samples, dtype=int)
whole_list.append(plot_samples)

#hardcoded list
end_Status = 4
j=0

#poor use of a while loop when a for loop should have been used but it works the same. 
#this is to make the differnt starts. Graphs will be just based on samples
while j<end_Status:

    num_samples.append(N_samples)

    unc_volf = unc*vol # +/- g/cc
    unc_rhoA_S = unc*rhoA_s # +/- g/cc

    #change via cases
    #holds rhoAs constant
    if j==0:
        unc_rhoA_S = 0
    #holds volf constant
    elif j==1:
        unc_volf = 0
    #holds both constant
    elif j==3:
        unc_volf = 0
        unc_rhoA_S = 0

    rhoA_uncertainty.append('%.*g' % (sigfigs,unc_rhoA_S))
    volf_uncertainty.append('%.*g' % (sigfigs,unc_volf))

    volf_a = np.random.normal(vol, unc_volf, N_samples)
    rhoAs_a = np.random.normal(rhoA_s, unc_rhoA_S, N_samples)

    L_conc = []
    L_conct = []
    L_mass = []
    L_masst = []

    for i in range(N_samples):
        rhoA_f = 0 #initial concentration
        massA_f = rhoA_f*vol #g
        #set up loop
        icount = 0
        time = timei

        while time<timef:

            rhoA_s = rhoAs_a[i]
            J_af = Jaf_funct(Km,rhoA_f,rhoA_s)
            vol = volf_a[i]
            
            drhoA = J_af*p_surf*Num_part/vol
            rhoA_n = rhoA_f + drhoA*dt

            rhoA_f = rhoA_n
            massA_f = rhoA_f * vol

            icount += 1
            time = icount*dt
        # print(f'concentration after ten minutes: {rhoA_f:.4g} g/cc')
        # print(f'concentration after ten minutes: {massA_f:.4g} g')        
        L_conc.append(rhoA_f)
        L_conct.append('%.*g' % (sigfigs,rhoA_f))
        L_mass.append(massA_f)
        L_masst.append('%.*g' % (sigfigs,massA_f))
    
    whole_list.append(L_mass)

    #print(L_masst)
    mean_conc = np.mean(L_conc)
    stdev_conc = np.std(L_conc)
    # print(mean_conc)
    # print(stdev_conc)

    stDev_List.append(stdev_conc)
    concMean_list.append(mean_conc)

    #formatted lists for the table
    stDev_ListF.append('%.*g' % (sigfigs,stdev_conc))
    concMean_listF.append('%.*g' % (sigfigs,mean_conc))

    mean_mass_conc = np.mean(L_mass)
    stdev_mass_conc = np.std(L_mass)

    stDev_List_M.append(stdev_mass_conc)
    concMean_listM.append(mean_mass_conc)

    #formatted list for table
    stDev_List_MF.append('%.*g' % (sigfigs,stdev_mass_conc))
    concMean_listMF.append('%.*g' % (sigfigs,mean_mass_conc))

    j+=1

# print(stDev_List_MF)
# print(concMean_listMF)
# print(volf_uncertainty)
# print(rhoA_uncertainty)
# print(num_samples)

#finally hard print a list to lable
listNames = ["N/A","Volf","Rho As","Both"]

#create dictionary
results = {"Held Constant":listNames,"Samples":num_samples,"Rho As (+/-)":rhoA_uncertainty,"Volf (+/-)":volf_uncertainty,"Result (g)":concMean_listMF,"Stdev (+/-)":stDev_List_MF}

#create df
df1 = pandas.DataFrame(results)
df1.set_index("Held Constant",inplace = True)
print(df1)
print("")
#print(whole_list)

#chat gpt with evan to make a subplot of 4
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

axes[0, 0].plot(whole_list[0], whole_list[1], 'ro')
axes[0, 0].set_title("Mass Concentration with neither constant", fontsize=16)
axes[0, 0].set_xlabel("Sample #", fontsize=14)
axes[0, 0].set_ylabel("Mass Concentration (g)", fontsize=14)
axes[0, 0].tick_params(labelsize=12)
axes[0, 0].grid(True)

axes[0, 1].plot(whole_list[0], whole_list[2], 'bo')
axes[0, 1].set_title("Mass Concentration with volf constant", fontsize=16)
axes[0, 1].set_xlabel("Sample #", fontsize=14)
axes[0, 1].set_ylabel("Mass Concentration (g)", fontsize=14)
axes[0, 1].tick_params(labelsize=12)
axes[0, 1].grid(True)

axes[1, 0].plot(whole_list[0], whole_list[3], 'go')
axes[1, 0].set_title("Mass Concentration with rho As constant", fontsize=16)
axes[1, 0].set_xlabel("Sample #", fontsize=14)
axes[1, 0].set_ylabel("Mass Concentration (g)", fontsize=14)
axes[1, 0].tick_params(labelsize=12)
axes[1, 0].grid(True)

axes[1, 1].plot(whole_list[0], whole_list[4], 'ro')
axes[1, 1].set_title("Mass Concentration with both constant", fontsize=16)
axes[1, 1].set_xlabel("Sample #", fontsize=14)
axes[1, 1].set_ylabel("Mass Concentration (g)", fontsize=14)
axes[1, 1].tick_params(labelsize=12)
axes[1, 1].grid(True)

plt.tight_layout()
plt.show()

# Make sure to change the write path in your own code
fig.savefig('Monte_Carlo/MonteCarloGraphsNotIn.jpeg',dpi = 300,bbox_inches = 'tight')

L_dfs = [df1]
with open('Monte_Carlo/MCResultsTableNotIn.csv','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")