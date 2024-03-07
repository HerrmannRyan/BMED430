import numpy as np
import matplotlib.pyplot as plt
import pandas

#setup m matrix finite diff metrix
#setup c matrix solution matrix
L = 1 #m
sigfigs = 4
n = 3
n1 = n+1
n2 = n1 + 1
dx = L/n1

mMat = np.zeros((n2,n2))
cMat = np.zeros((n2))

mMat[n1,n1] = 1
mMat[0,0] = 1
cMat[-1] = 100.00

LoHS = 0.0
L_xp = [0]
L_xpf = ['%.*g'%(sigfigs,LoHS)]

#append format so that I dont have to keep writing the same thing over and over
def L_xpfAppend(n):
    L_xpf.append('%.*g' % (sigfigs,n))

for i in range(1,n1):
    L_xp.append(i*dx)
    L_xpfAppend(i*dx)

    mMat[i,i] = -2
    mMat[i, i-1] = 1
    mMat[i, i+1] = 1
    cMat[i] = 0

#test print matrix
print(mMat)
print(cMat)

L_xp.append(L)
L_xpfAppend(L)

solve = np.linalg.solve(mMat,cMat)

#test print
print(solve)
print(L_xp)

fig = plt.figure(figsize = (10,8))
plt.plot(L_xp, solve, 'o', color = "tab:red")
plt.plot(L_xp, solve, '-', color = "tab:blue")
plt.title('Temperature vs Position', fontsize = 12)
plt.ylabel('Temperature $^\\circ$C', fontsize = 12)
plt.xlabel('Position (m)', fontsize = 12)
plt.grid(True)
plt.show()
fig.savefig('Heat_Transfer/ResultFig.jpeg',dpi = 300,bbox_inches = 'tight')

results = {'Position (m)': L_xpf, 'Temperature (C)': solve}
df1 = pandas.DataFrame(results)
print(df1)

with open('Heat_Transfer/TempTable.csv','w',newline='') as f:
    df1.to_csv(f)
    f.write("\n")