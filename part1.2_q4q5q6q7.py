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
m11 = popt[0]
A011 = popt[1]
m21 = popt[2]
A021 = popt[3]
axs[0,0].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A2,sigma=A2error)
print('A2',popt)
m12 = popt[0]
A012 = popt[1]
m22 = popt[2]
A022 = popt[3]
axs[1,0].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A3,sigma=A3error)
print('A3',popt)
m13 = popt[0]
A013 = popt[1]
m23 = popt[2]
A023 = popt[3]
axs[0,1].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")

popt,pcov = scipy.optimize.curve_fit(radioactivedecay,t,A4,sigma=A4error)
print('A4',popt)
m14 = popt[0]
A014 = popt[1]
m24 = popt[2]
A024 = popt[3]
axs[1,1].plot(t,radioactivedecay(t,*popt),'--b',label=f"scipy fit")
plt.show()

#the values are not similar to those found in q1
#this proves it is different from a single isotope and a mixture of isotopes

steps = 5

intervalm1 = 1e-5
intervalA01 = 2
intervalm2 = 1e-5
intervalA02 = 1 
intervalm3 = 1e-5
intervalA03 = 2
intervalm4 = 1e-5
intervalA04 = 1.5

gridm11 = np.arange(m11-intervalm1*steps,m11+(intervalm1)*(steps-0.5),intervalm1)
gridm21 = np.arange(m21-intervalm1*steps,m21+(intervalm1)*(steps-0.5),intervalm1)

gridm12 = np.arange(m12-intervalm2*steps,m12+(intervalm2)*(steps-0.5),intervalm2)
gridm22 =np.arange(m22-intervalm2*steps,m22+(intervalm2)*(steps-0.5),intervalm2)

gridm13 = np.arange(m13-intervalm3*steps,m13+(intervalm3)*(steps-0.5),intervalm3)
gridm23 = np.arange(m23-intervalm3*steps,m23+(intervalm3)*(steps-0.5),intervalm3)

gridm14 = np.arange(m14-intervalm4*steps,m14+(intervalm4)*(steps-0.5),intervalm4)
gridm24 = np.arange(m24-intervalm4*steps,m24+(intervalm4)*(steps-0.5),intervalm4)

gridA011 = np.arange(A011-intervalA01*steps,A011+(intervalA01)*steps,intervalA01)
gridA021 = np.arange(A021-intervalA01*steps,A021+(intervalA01)*steps,intervalA01)

gridA012 = np.arange(A012-intervalA02*steps,A012+(intervalA02)*steps,intervalA02)
gridA022 =np.arange(A022-intervalA02*steps,A022+(intervalA02)*(steps-0.5),intervalA02)

gridA013 = np.arange(A013-intervalA03*steps,A013+(intervalA03)*steps,intervalA03)
gridA023 = np.arange(A023-intervalA03*steps,A023+(intervalA03)*steps,intervalA03)

gridA014 = np.arange(A014-intervalA04*steps,A014+(intervalA04)*steps,intervalA04)
gridA024 = np.arange(A024-intervalA04*steps,A024+(intervalA04)*steps,intervalA04)

chi1 = []
chi2 = []
chi3 = []
chi4 = []

def decay(t,gridm1,gridA01,gridm2,gridA02):
  return gridA01*np.exp(-gridm1*t)+gridA02*np.exp(-gridm2*t)

#4d grid means that you need to append at certain points

for i in range(10):
  chi1.append([])
  for j in range(10):
     chi1[i].append([])
     for k in range(10):
      chi1[i][j].append([])
      for l in range(10):
        chi1[i][j][k].append([np.sum((((A1-decay(t,gridm11[i],gridA011[j],gridm21[k],gridA021[l]))**2))/A1error**2)])

chi1 = np.array(chi1)

print(np.shape(chi1))

for i in range(10):
  chi2.append([])
  for j in range(10):
     chi2[i].append([])
     for k in range(10):
      chi2[i][j].append([])
      for l in range(10):
        chi2[i][j][k].append([np.sum((((A2-decay(t,gridm12[i],gridA012[j],gridm22[k],gridA022[l]))**2))/A2error**2)])

chi2 = np.array(chi2)

for i in range(10):
  chi3.append([])
  for j in range(10):
     chi3[i].append([])
     for k in range(10):
      chi3[i][j].append([])
      for l in range(10):
        chi3[i][j][k].append([np.sum((((A3-decay(t,gridm13[i],gridA013[j],gridm23[k],gridA023[l]))**2))/A3error**2)])

chi3 = np.array(chi3)

for i in range(10):
  chi4.append([])
  for j in range(10):
     chi4[i].append([])
     for k in range(10):
      chi4[i][j].append([])
      for l in range(10):
        chi4[i][j][k].append([np.sum((((A4-decay(t,gridm14[i],gridA014[j],gridm24[k],gridA024[l]))**2))/A4error**2)])

chi4 = np.array(chi4)

minchi1 = np.amin(chi1)
minchi2 = np.amin(chi2)
minchi3 = np.amin(chi3)
minchi4 = np.amin(chi4)

# print(minchi1,minchi2,minchi3,minchi4)

#finding uncertainties

index1 = np.unravel_index(np.argmin(chi1,axis=None),chi1.shape)
index2 = np.unravel_index(np.argmin(chi2,axis=None),chi2.shape)
index3 = np.unravel_index(np.argmin(chi3,axis=None),chi3.shape)
index4 = np.unravel_index(np.argmin(chi4,axis=None),chi4.shape)

#m11
# plt.plot(gridm11,chi1[index1[0],index1[1],index1[2],:])
# plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
sigma1 = np.around((chi1[index1[0],index1[1],index1[2],:] - (minchi1+2.3) ),decimals=0)
sigma1 = np.where(sigma1 == 0)[0]
# print('decayconstant11', gridm11[index1[0]], '+/-', np.abs(gridm11[index1[0]]-gridm11[sigma1[0]]))
# plt.show()

#m11
# plt.plot(gridm11,chi1[index1[0],index1[1],index1[2],:])
# plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
# print('decayconstant11', gridm11[index1[0]], '+/-', np.abs(gridm11[index1[0]]-gridm11[sigma1[0]]))
# plt.show()

#m21
# plt.plot(gridm21,chi1[index1[0],index1[1],index1[2],:])
# plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
# print('decayconstant21', gridm21[index1[0]], '+/-', np.abs(gridm21[index1[0]]-gridm21[sigma1[0]]))
# plt.show()

#A011
# plt.plot(gridA011,chi1[index1[0],index1[1],index1[2],:])
# plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
# print('A011', gridA011[index1[0]], '+/-', np.abs(gridA011[index1[0]]-gridA011[sigma1[0]]))
# plt.show()

#A021
# plt.plot(gridA021,chi1[index1[0],index1[1],index1[2],:])
# plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
# print('A021', gridA021[index1[0]], '+/-', np.abs(gridA021[index1[0]]-gridA021[sigma1[0]]))
# plt.show()

#m12
# plt.plot(gridm12,chi2[index2[0],index2[1],index2[2],:])
# plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
sigma2 = np.around((chi2[index2[0],index2[1],index2[2],:] - (minchi2+2.3) - 0.5),decimals=0)
sigma2 = np.where(sigma2 == 0)[0]
# print('decayconstant12', gridm12[index2[0]], '+/-', np.abs(gridm12[index2[0]]-gridm12[sigma2[0]]))
# plt.show()

#m22
# plt.plot(gridm22,chi2[index2[0],index2[1],index2[2],:])
# plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
# print('decayconstant22', gridm22[index2[0]], '+/-', np.abs(gridm22[index2[0]]-gridm22[sigma2[0]]))
# plt.show()

#A012
# plt.plot(gridA022,chi2[index2[0],index2[1],index2[2],:])
# plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
# print('A012', gridA012[index2[0]], '+/-', np.abs(gridA012[index2[0]]-gridA012[sigma2[0]]))
# plt.show()

#A022
# plt.plot(gridA022,chi2[index2[0],index2[1],index2[2],:])
# plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
# print('A022', gridA022[index2[0]], '+/-', np.abs(gridA022[index2[0]]-gridA022[sigma2[0]]))
# plt.show()

#m13
# plt.plot(gridm13,chi3[index3[0],index3[1],index3[2],:])
# plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
sigma3 = np.around((chi3[index3[0],index3[1],index3[2],:] - (minchi3+2.3) ),decimals=0)
sigma3 = np.where(sigma3 == 0)[0]
# print('decayconstant13', gridm13[index3[0]], '+/-', np.abs(gridm13[index3[0]]-gridm13[sigma3[0]]))
# plt.show()

#m23
# plt.plot(gridm23,chi3[index3[0],index3[1],index3[2],:])
# plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
# print('decayconstant23', gridm23[index3[0]], '+/-', np.abs(gridm23[index3[0]]-gridm23[sigma3[0]]))
# plt.show()

#A013
# plt.plot(gridA013,chi3[index3[0],index3[1],index3[2],:])
# plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
# print('A013', gridA013[index3[0]], '+/-', np.abs(gridA013[index3[0]]-gridA013[sigma3[0]]))
# plt.show()

#A023
# plt.plot(gridA023,chi3[index3[0],index3[1],index3[2],:])
# plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
# print('A023', gridA023[index3[0]], '+/-', np.abs(gridA023[index3[0]]-gridA023[sigma3[0]]))
# plt.show()

# m14
# plt.plot(gridm14,chi4[index4[0],index4[1],index4[2],:])
# plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
sigma4 = np.around((chi4[index4[0],index4[1],index4[2],:] - (minchi4+2.3) ),decimals=0)
sigma4 = np.where(sigma4 == 0)[0]
# print('decayconstant14', gridm14[index4[0]], '+/-', np.abs(gridm14[index4[0]]-gridm14[sigma4[0]]))
# plt.show()

#m24
# plt.plot(gridm24,chi4[index4[0],index4[1],index4[2],:])
# plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
# print('decayconstant24', gridm24[index4[0]], '+/-', np.abs(gridm24[index4[0]]-gridm24[sigma4[0]]))
# plt.show()

#A014
# plt.plot(gridA014,chi4[index4[0],index4[1],index4[2],:])
# plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
# print('A014', gridA014[index4[0]], '+/-', np.abs(gridA014[index4[0]]-gridA014[sigma4[0]]))
# plt.show()

#A024
# plt.plot(gridA024,chi4[index4[0],index4[1],index4[2],:])
# plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
# print('A024', gridA024[index4[0]], '+/-', np.abs(gridA024[index4[0]]-gridA024[sigma4[0]]))
# plt.show()

#AIC and BIC

# p = 4 #p is the number of parameters - the parameters are A01, A02, decayconstant1 and decayconstant2
# N1 = len(A1) #N is number of data points of R1/ R2/ R3/ R4
# N2 = len(A2)
# N3 = len(A3)
# N4 = len(A4)

# AIC1 = minchi1 + 2*p
# BIC1 = minchi1 + np.log(N1)*p

# AIC2 = minchi2 + 2*p
# BIC2 = minchi2 + np.log(N2)*p

# AIC3 = minchi3 + 2*p
# BIC3 = minchi3 + np.log(N3)*p

# AIC4 = minchi4 + 2*p
# BIC4 = minchi4 + np.log(N4)*p

# print(AIC1,BIC1)
# print(AIC2,BIC2)
# print(AIC3,BIC3)
# print(AIC4,BIC4)

#question 5

#these AIC and BIC values are much smaller than the values from q2
#does this mean the two isoptope model is better?

#question 6

#need to find half life then compare to table
#halflife = 0.693/decayconstant

halflife11 = 0.693/gridm11[index1[0]]
halflife21 = 0.693/gridm21[index1[0]]
halflife12 = 0.693/gridm12[index2[0]]
halflife22 = 0.693/gridm22[index2[0]]
halflife13 = 0.693/gridm13[index3[0]]
halflife23 = 0.693/gridm23[index3[0]]
halflife14 = 0.693/gridm14[index4[0]]
halflife24 = 0.693/gridm24[index4[0]]

#how do i find the uncertainty im dumb

# halferror11 = 0.693/(np.abs(gridm11[index1[0]]-gridm11[sigma1[0]]))
# halferror21 = 0.693/(np.abs(gridm12[index1[0]]-gridm12[sigma1[0]]))
# halferror12 = 0.693/(np.abs(gridm21[index2[0]]-gridm21[sigma2[0]]))
# halferror22 = 0.693/(np.abs(gridm22[index2[0]]-gridm22[sigma2[0]]))
# halferror13 = 0.693/(np.abs(gridm13[index3[0]]-gridm13[sigma3[0]]))
# halferror23 = 0.693/(np.abs(gridm23[index3[0]]-gridm23[sigma3[0]]))
# halferror14 = 0.693/(np.abs(gridm14[index4[0]]-gridm14[sigma4[0]]))
# halferror24 = 0.693/(np.abs(gridm24[index4[0]]-gridm24[sigma4[0]]))


# print('halflife11', halflife11, '+/-', halferror11)
# print('halflife21', halflife12, '+/-', halferror21)
# print('halflife12', halflife21, '+/-', halferror12)
# print('halflife22', halflife22, '+/-', halferror22)
# print('halflife13', halflife13, '+/-', halferror13)
# print('halflife23', halflife23, '+/-', halferror23)
# print('halflife14', halflife14, '+/-', halferror14)
# print('halflife24', halflife24, '+/-', halferror24)

#they are annoyingly different values

#question 7
#idk how to change the equation to fit this, is the way i did it okay?

# all distance in metres
d = np.array([0.0500, 0.0500, 0.1000, 0.1000, 0.1800, 0.1800, 0.3000, 0.3000])
derror = np.array([0.0005,0.0005,0.0005,0.0005])

A0 = np.array([gridA011[index1[0]],gridA021[index1[0]],gridA012[index2[0]],gridA022[index2[0]],gridA013[index3[0]],gridA023[index3[0]],gridA014[index4[0]],gridA024[index4[0]]])
A0error = np.array([np.abs(gridA011[index1[0]]-gridA011[sigma1[0]]),np.abs(gridA021[index1[0]]-gridA021[sigma1[0]]),np.abs(gridA012[index2[0]]-gridA012[sigma2[0]]),np.abs(gridA022[index2[0]]-gridA022[sigma2[0]]),np.abs(gridA013[index3[0]]-gridA013[sigma3[0]]),np.abs(gridA023[index3[0]]-gridA023[sigma3[0]]),np.abs(gridA014[index4[0]]-gridA014[sigma4[0]]),np.abs(gridA024[index4[0]]-gridA024[sigma4[0]])])

#plot ln(A0) on y and ln(d) on x, ln(c) is intercept
y = np.log(A0)

yerror = A0error/A0

x = np.log(d)
# xerror = derror/d 

print(y)
print(x)

def leastsquares(x,y,error):
  x = np.array(x)
  y = np.array(y)
  error = np.array(error)

  s = np.sum(1/(error**2))
  sx = np.sum(x/(error**2))
  sy = np.sum(y/(error**2))
  sxy = np.sum((x*y)/(error**2))
  sxx = np.sum((x**2)/(error**2))

  m = (s*sxy-sx*sy)/(s*sxx-sx**2)
  c = (sxx*sy-sx*sxy)/(s*sxx-sx**2)

  return m,c

#what about the uncertainty in d
m,c = leastsquares(x, y, yerror) 

#why are my error bars so small?
y2 = m * x + c
plt.plot(x,y2,color='red')
plt.plot(x,y,color='red',ls='none')
plt.errorbar(x,y,yerr=yerror,fmt='ro')
plt.show()

print(m)
print(np.exp(c))

#need to find uncertainties - using chi squared instead of bootstrap

steps = 50

intervalm = 0.1
intervalc = 0.1

gridm = np.arange(m-intervalm*steps,m+(intervalm)*(steps-1),intervalm)
gridc = np.arange(m-intervalm*steps,m+(intervalm)*(steps-1),intervalm)

chi = []

def isotropic(t,gridm,gridA0):
  return gridc * d**gridm

for i in range(len(gridm)):
  chi.append([])
  for j in range(len(gridc)):
    chi[i].append(np.sum((((A0-isotropic(d,gridm[i],gridc[j]))**2))/A0error**2))

chi = np.array(chi)

minchi = np.amin(chi)

index = np.unravel_index(np.argmin(chi,axis=None),chi.shape)

plt.plot(gridm,chi[index[0],:])
plt.axhline(y=(minchi+2.3), color='r', linestyle='-')
sigma = np.argwhere(np.diff(np.sign(chi[index[0],:] - (minchi+2.3)))).flatten()
print('powerlaw index', gridm[index1[0]], '+/-',np.abs(gridm[index1[0]]-((gridm[sigma1[0]]+gridm[sigma1[1]])/2)))
plt.show()

#need to find how many particles are in the sample





