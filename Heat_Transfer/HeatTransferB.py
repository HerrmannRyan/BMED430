import numpy as np
import matplotlib.pyplot as plt
import pandas

#setup m matrix finite diff metrix
#setup c matrix solution matrix
L = 1 #m
sigfigs = 4
epi = 1e-5 #convergence criteria
m_err = 10.0 #max error
LeftHandSide = 0.0
RightHandSide = 100.0

n = 38
n1 = n+1
n2 = n1 + 1
dx = L/(n2-1)

mMat = np.zeros((n2,n2))
cMat = np.zeros(n2)
u1 = np.zeros(n2)
u1n = np.zeros(n2)

L_xp = [0]
L_xpf = ['%.*f' % (sigfigs,LeftHandSide)]

u1f = []

#append format so that I dont have to keep writing the same thing over and over
def L_xpfAppend(n):
    L_xpf.append('%.*f' % (sigfigs,n))

for i in range(1,n1):
    mMat[i,i] = -2
    mMat[i, i-1] = 1
    mMat[i, i+1] = 1
    cMat[i] = 0
    L_xp.append(i*dx)
    L_xpfAppend(i*dx)

mMat[n1,n1] = 1
mMat[0,0] = 1
cMat[-1] = 100.00
u1[0] = 0.0
u1[-1] = 100

#test print matrix
#print(mMat)
#print(cMat)

L_xp.append(L)
L_xpfAppend(L)

linsolve_solve = np.linalg.solve(mMat,cMat)

icount = 0
L_merr = []
L_err = []
L_count = []

#gauss sidell method use u1 and u1n
while m_err > epi:
    icount += 1

    for j in range(1,n1):
        u1n[j] = (1/mMat[j,j])*(cMat[j] - mMat[j,j-1]*u1[j-1] - mMat[j,j+1]*u1[j+1])
        err = np.abs(u1n[j]-u1[j]) #absolute error
        L_err.append(err)
        #print(u1[j])
        u1[j] = u1n[j]
    
    m_err = max(L_err)
    L_err = []
    L_merr.append('%.*g' % (sigfigs,m_err))
    L_count.append(icount)


#test print
print(linsolve_solve)
print(u1)
print(icount)

L_aberr = [] #absolute error
for i in range(0,n2):
    u1f.append('%.*f' % (sigfigs-1,u1[i]))
    L_aberr.append(np.abs(linsolve_solve[i]-u1[i]))


fig = plt.figure(figsize = (10,8))

plt.subplot(2,1,1)
plt.plot(L_xp, linsolve_solve, 'o', color = "tab:red")
plt.plot(L_xp, linsolve_solve, '--', color = "tab:blue")
plt.title('LU Decomp Temperature vs Position', fontsize = 12)
plt.ylabel('Temperature $^\\circ$C', fontsize = 12)
plt.xlabel('Position (m)', fontsize = 12)
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(L_xp, u1, 'o', color = "tab:red")
plt.plot(L_xp, u1, '--', color = "tab:blue")
plt.title('GS Temperature vs Position', fontsize = 12)
plt.ylabel('Temperature $^\\circ$C', fontsize = 12)
plt.xlabel('Position (m)', fontsize = 12)
plt.grid(True)

plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9, hspace = .4)
#plt.show()


fig.savefig('Heat_Transfer/Gaussidell6b.jpeg',dpi = 300,bbox_inches = 'tight')

#Data Frame
results = {'Position (m)': L_xpf, 'LU Decomp Temperature (C)': linsolve_solve, 'GS Temperature (C)':u1f, 'Abs Error':L_aberr}

#checks for value and then adds percentage symbol on top
if 'Abs Error' in results and all(isinstance(val, (int, float)) for val in results['Abs Error']):
    results['Abs Error'] = [f"{val*100:.4f}%" for val in results['Abs Error']]

if 'LU Decomp Temperature (C)' in results and all(isinstance(val, (int, float)) for val in results['LU Decomp Temperature (C)']):
    results['LU Decomp Temperature (C)'] = [f"{val:.3f}" for val in results['LU Decomp Temperature (C)']]

df1 = pandas.DataFrame(results)
df1.set_index("Position (m)",inplace=True)
print(df1)

#Iterations DataFrame
stats_dict = {'Iterations':L_count, 'Max Error':L_merr}
df_stats = pandas.DataFrame(stats_dict)
df_stats.set_index('Iterations',inplace=True)
print(df_stats)

L_header = ["Iterations", "Absolute Error (degC)", "Convergence Criteria"]
L_stats = [L_count[-1], L_merr[-1], epi]
stats_dict2 = {'Results': L_stats}
df_stats2 = pandas.DataFrame(stats_dict2)
df_stats2.index = L_header
print(df_stats2)

L_dfs = [df1,df_stats,df_stats2]
with open('Heat_Transfer/TempTable6b','w',newline='') as f:
    for df in L_dfs:
        df.to_csv(f)
        f.write("\n")