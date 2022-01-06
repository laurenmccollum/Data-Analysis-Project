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

#plot ln(A) on y and t on x, A0 is intercept
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
m1 = m
A01 = np.exp(c)

m,c = leastsquares(t, y2, y2error) 
ynew2 = m * t + c
m2 = m
A02 = np.exp(c)

m,c = leastsquares(t, y3, y3error) 
ynew3 = m * t + c
m3 = m
A03 = np.exp(c)

m,c = leastsquares(t, y4, y4error) 
ynew4 = m * t + c
m4 = m
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

print(m1,m2,m3,m4)
print(A01,A02,A03,A04)
