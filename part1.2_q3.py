import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii
import scipy
from scipy import optimize

def radioactivedecay(t,decayconstant1,A01,decayconstant2,A02):
    return np.array(A01*np.exp(-decayconstant1*t)+A02*np.exp(-decayconstant2*t))

# here's how to read in some data: 
data = ascii.read('part1_data2.dat') 

t = np.array(data['Time'])

A1 = np.array(data['R1'])/t[1]
A2 = np.array(data['R2'])/t[1]
A3 = np.array(data['R3'])/t[1]
A4 = np.array(data['R4'])/t[1]

A1error = np.array(data['e_R1'])/t[1]
A2error = np.array(data['e_R2'])/t[1]
A3error = np.array(data['e_R3'])/t[1]
A4error =np.array(data['e_R4'])/t[1]

fig, axs = plt.subplots(2,2)

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A1,sigma=A1error)
print('A1',popt)
axs[0,0].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A2,sigma=A2error)
print('A2',popt)
axs[1,0].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A3,sigma=A3error)
print('A3',popt)
axs[0,1].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A4,sigma=A4error)
print('A4',popt)
axs[1,1].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")
plt.show()

#the values are not similar to those found in q1
#this proves it is different from a single isotope and a mixture of isotopes


