import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii

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

#plot ln(A) on y and t on x, ln(A0) is intercept and lamda is gradient
y1 = np.log(A1)
y1error = A1error/A1 #log error only works when delta x is much smaller than x

y2 = np.log(A2)
y2error = A2error/A2

y3 = np.log(A3)
y3error = A3error/A3

y4 = np.log(A4)
y4error = A4error/A4

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

m,c = leastsquares(t, y1, y1error) 
ynew1 = m * t + c
m1 = -m #note i have change it to -m because wanted to make m positive, easier to work with i guess
A01 = np.exp(c)

m,c = leastsquares(t, y2, y2error) 
ynew2 = m * t + c
m2 = -m
A02 = np.exp(c)

m,c = leastsquares(t, y3, y3error) 
ynew3 = m * t + c
m3 = -m
A03 = np.exp(c)

m,c = leastsquares(t, y4, y4error) 
ynew4 = m * t + c
m4 = -m
A04 = np.exp(c)

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(t, y1,ls='none')
axs[0,0].errorbar(t,y1,yerr=y1error,fmt='ro')
axs[0,0].plot(t,ynew1,color='black')
axs[0, 0].set_title('Axis [0, 0]')

axs[0, 1].plot(t, y2, 'tab:orange',ls='none')
axs[0,1].errorbar(t,y2,yerr=y2error,fmt='ro')
axs[0,1].plot(t,ynew2,color='black')
axs[0, 1].set_title('Axis [0, 1]')

axs[1, 0].plot(t, y3, 'tab:green',ls='none')
axs[1,0].errorbar(t,y3,yerr=y3error,fmt='ro')
axs[1,0].plot(t,ynew3,color='black')
axs[1, 0].set_title('Axis [1, 0]')

axs[1, 1].plot(t, y4, 'tab:red',ls='none')
axs[1,1].errorbar(t,y4,yerr=y4error,fmt='ro')
axs[1,1].plot(t,ynew4,color='black')
axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()

# print(m1,m2,m3,m4)
# print(A01,A02,A03,A04)



#how to make a chi^2 grid



#two lists - one decay constants and one A0
#make the values centred around the value found in previous part

steps = 50

intervalm1 = 1e-5
intervalA01 = 0.7
intervalm2 = 1e-5
intervalA02 = 0.2
intervalm3 = 1e-5
intervalA03 = 0.15
intervalm4 = 1e-5
intervalA04 = 0.1

gridm1 = np.arange(m1-intervalm1*steps,m1+(intervalm1)*(steps-1),intervalm1)
gridm2 = np.arange(m2-intervalm2*steps,m2+(intervalm2)*(steps-1),intervalm2)
gridm3 = np.arange(m3-intervalm3*steps,m3+(intervalm3)*(steps-1),intervalm3)
gridm4 =np.arange(m4-intervalm4*steps,m4+(intervalm4)*(steps),intervalm4)

# print(len(gridm1),len(gridm2),len(gridm3),len(gridm4))

gridA01 = np.arange(A01-intervalA01*steps,A01+(intervalA01)*(steps),intervalA01)
gridA02 = np.arange(A02-intervalA02*(steps-1),A02+(intervalA02)*(steps+1),intervalA02)
gridA03 = np.arange(A03-intervalA03*(steps-1),A03+(intervalA03)*(steps),intervalA03)
gridA04 = np.arange(A04-intervalA04*steps,A04+(intervalA04)*(steps),intervalA04)

# print(len(gridA01),len(gridA02),len(gridA03),len(gridA04))

#calculate chi squared for each point on the grid (you will need a nested loop)

chi1 = []
chi2 = []
chi3 = []
chi4 = []

def decay(t,gridm,gridA0):
  return gridA0 * np.exp(-gridm*t)

for i in range(len(gridm1)):
  chi1.append([])
  for j in range(len(gridA01)):
    chi1[i].append(np.sum((((A1-decay(t,gridm1[i],gridA01[j]))**2))/A1error**2))

chi1 = np.array(chi1)

for i in range(len(gridm2)):
  chi2.append([])
  for j in range(len(gridA02)):
    chi2[i].append(np.sum(((A2-decay(t,gridm2[i],gridA02[j]))**2)/A2error**2))

chi2 = np.array(chi2)

for i in range(len(gridm3)):
  chi3.append([])
  for j in range(len(gridA03)):
    chi3[i].append(np.sum(((A3-decay(t,gridm3[i],gridA03[j]))**2)/A3error**2))

chi3 = np.array(chi3)

for i in range(len(gridm4)):
  chi4.append([])
  for j in range(len(gridA04)):
    chi4[i].append(np.sum(((A4-decay(t,gridm4[i],gridA04[j]))**2)/A4error**2))

chi4 = np.array(chi4)

minchi1 = np.amin(chi1)
minchi2 = np.amin(chi2)
minchi3 = np.amin(chi3)
minchi4 = np.amin(chi4)

index1 = np.unravel_index(np.argmin(chi1,axis=None),chi1.shape)
index2 = np.unravel_index(np.argmin(chi2,axis=None),chi2.shape)
index3 = np.unravel_index(np.argmin(chi3,axis=None),chi3.shape)
index4 = np.unravel_index(np.argmin(chi4,axis=None),chi4.shape)

plt.plot(gridm1,chi1[index1[0],:])
plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
sigma1 = np.argwhere(np.diff(np.sign(chi1[index1[0],:] - (minchi1+2.3)))).flatten()
print('decayconstant1', gridm1[index1[0]], '+/-',np.abs(gridm1[index1[0]]-((gridm1[sigma1[0]]+gridm1[sigma1[1]])/2)))
plt.show()

plt.plot(gridA01,chi1[index1[0],:])
plt.axhline(y=(minchi1+2.3), color='r', linestyle='-')
print('A01', gridA01[index1[0]], '+/-', np.abs(gridA01[index1[0]]-((gridA01[sigma1[0]]+gridA01[sigma1[1]])/2)))
plt.show()

plt.plot(gridm2,chi2[index2[0],:])
plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
sigma2 = np.argwhere(np.diff(np.sign(chi2[index2[0],:] - (minchi2+2.3)))).flatten()
print('decayconstant2', gridm2[index2[0]], '+/-', np.abs(gridm2[index2[0]]-((gridm2[sigma2[0]]+gridm2[sigma2[1]])/2)))
plt.show()

plt.plot(gridA02,chi2[index2[0],:])
plt.axhline(y=(minchi2+2.3), color='r', linestyle='-')
print('A02', gridA02[index2[0]], '+/-', np.abs(gridA02[index2[0]]-((gridA02[sigma2[0]]+gridA02[sigma2[1]])/2)))
plt.show()

plt.plot(gridm3,chi3[index3[0],:])
plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
sigma3 = np.argwhere(np.diff(np.sign(chi3[index3[0],:] - (minchi3+2.3)))).flatten()
print('decayconstant3', gridm3[index3[0]], '+/-', np.abs(gridm3[index3[0]]-((gridm3[sigma3[0]]+gridm3[sigma3[1]])/2)))
plt.show()

plt.plot(gridA03,chi3[index3[0],:])
plt.axhline(y=(minchi3+2.3), color='r', linestyle='-')
print('A03', gridA03[index3[0]], '+/-', np.abs(gridA03[index3[0]]-((gridA03[sigma3[0]]+gridA03[sigma3[1]])/2)))
plt.show()

plt.plot(gridm4,chi4[index4[0],:])
plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
sigma4 = np.argwhere(np.diff(np.sign(chi4[index4[0],:] - (minchi4+2.3)))).flatten()
sigma4 = sigma4[1]
print('decayconstant4', gridm4[index4[0]], '+/-', np.abs(gridm4[index4[0]]-gridm4[sigma4]))
plt.show()

plt.plot(gridA04,chi4[index4[0],:])
plt.axhline(y=(minchi4+2.3), color='r', linestyle='-')
print('A04', gridA04[index4[0]], '+/-', np.abs(gridA04[index4[0]]-gridA04[sigma4]))
plt.show()

#need to change graphs so that i can find both intersection points
#then the final uncertainty will be the average of those two uncertainties from the two points

#point where its minimum chi - intersection point = uncertainty

#explainvideo
#so you have the grid with all the values of chi squared
#the place in the grid where the minimum chi squared is, take its row and column
#plot them i guess
#the chi squared value should go up by 2.3 for 68% (this is the 1 sigma confidence interval so this should be our error)
#whatever value of decay constant which makes the chi squared value go up by 2.3, the difference between that and minimum chi is the error

#finding AIC and BIC

p = 2 #p is the number of parameters - the parameters are A0 and decay constant
N1 = len(A1) #N is number of data points of R1/ R2/ R3/ R4
N2 = len(A2)
N3 = len(A3)
N4 = len(A4)


AIC1 = minchi1 + 2*p
BIC1 = minchi1 + np.log(N1)*p

AIC2 = minchi2 + 2*p
BIC2 = minchi2 + np.log(N2)*p

AIC3 = minchi3 + 2*p
BIC3 = minchi3 + np.log(N3)*p

AIC4 = minchi4 + 2*p
BIC4 = minchi4 + np.log(N4)*p

# print(AIC1,BIC1)
# print(AIC2,BIC2)
# print(AIC3,BIC3)
# print(AIC4,BIC4)

#no idea if AIC and BIC values are correct or if they make sense


#question 6

# halflife1 = 0.693/m1
# halflife2 = 0.693/m2
# halflife3 = 0.693/m3
# halflife4 = 0.693/m4

# print(halflife1,halflife2,halflife3,halflife4)
#much more similar to each other therefore should i go with this model?
