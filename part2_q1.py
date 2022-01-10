import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii
from GLS import GLS as GLS
from scipy import optimize

#Data for RV of system %i circuit
#t time in days (JD-2,400,000)
#RV radial velocity in m/s
#e_RV uncertainty on RV in m/s

#question 1

data = ascii.read('part2_data1.dat') 

t = np.array(data['t'])
# t = t*24*60*60
RV = np.array(data['RV'])
RVerror = np.array(data['e_RV'])

# plt.plot(t,RV)
# plt.show()

# f = np.fft.fftfreq(t.shape[0])

P = np.logspace(-1.3,np.log10(30./2.),int(1e5))

p_om,a_om,b_om=GLS(t,RV,RVerror,P) 

# plt.plot(P,p_om)
# plt.show()

#question 2
index = np.argmax(p_om)
strongest = P[index]
# print(strongest)

#question 3
phi = (t/strongest) %1
# plt.plot(phi,RV,'o')
# plt.show()

steps = 50

#get K and D values from looking at the graph
K = 55
P = strongest
D = 3.8

intervalK = 0.2
intervalP = 4.05e-6
intervalD = 0.015

gridK = np.arange(K-intervalK*steps,K+(intervalK)*(steps-1.5),intervalK)
gridP = np.arange(P-intervalP*steps,P+(intervalP)*(steps-1.5),intervalP)
gridD = np.arange(D-intervalD*steps,D+(intervalD)*(steps-1.5),intervalD)

#calculate chi squared for each point on the grid (you will need a nested loop)

chi = []

def radialvelocity(K,t,P,D):
  return K*np.sin(2*np.pi*((t/P)-D))

for i in range(99):
  chi.append([])
  for j in range(99):
     chi[i].append([])
     for k in range(99):
      chi[i][j].append(np.sum((((RV-radialvelocity(gridK[i],t,gridP[j],gridD[k],))**2))/RVerror**2))

chi = np.array(chi)

minchi = np.amin(chi)

index = np.unravel_index(np.argmin(chi,axis=None),chi.shape)

plt.plot(gridK,chi[:,index[1],index[2]])
plt.axhline(y=(minchi+3.53), color='r', linestyle='-')
sigma = np.argwhere(np.diff(np.sign((np.min(chi[i:,index[1],index[2]]) + 3.53) - chi[:,index[1],index[2]]))).flatten()
print('K', gridK[index[0]], '+/-',np.abs(gridK[index[0]]-(((gridK[sigma[0]])))))
plt.show()

print(sigma)

plt.plot(gridP,chi[index[0],:,index[2]])
plt.axhline(y=(minchi+3.53), color='r', linestyle='-')
sigma = np.argwhere(np.diff(np.sign((np.min(chi[index[0], :, index[2]]) + 3.53) - chi[index[0], :, index[2]]))).flatten()
print('P', gridP[index[1]], '+/-', np.abs(gridP[index[1]]-((gridP[sigma[0]]+gridP[sigma[1]])/2)))
plt.show()

print(sigma)

plt.plot(gridD,chi[index[0],index[1],:])
plt.axhline(y=(minchi+3.53), color='r', linestyle='-')
sigma = np.argwhere(np.diff(np.sign((np.min(chi[index[0],index[1],:]) + 3.53) - chi[index[0],index[1],:]))).flatten()
print('delta phi', gridD[index[2]], '+/-', np.abs(gridD[index[2]]-((gridD[sigma[0]]+gridD[sigma[1]])/2)))
plt.show()

print(sigma)

#for 3 its 3.53
#for 4 its 4.72

#question 4

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A1,sigma=A1error)
print('A1',popt)
