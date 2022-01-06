#detector cross section = 3.142 +/- 0.001 cm^2

import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii

# here's how to read in some data: 
data = ascii.read('part1_data1.dat') 

A1 = np.array(data['R1'])
A2 = np.array(data['R2'])
A3 = np.array(data['R3'])
A4 = np.array(data['R4'])

t = np.array(data['Time'])

A1error = np.array(data['e_R1'])
A2error = np.array(data['e_R2'])
A3error = np.array(data['e_R3'])
A4error =np.array(data['e_R4'])

#plot ln(A) on y and t on x, A0 is intercept
y1 = np.log(A1)
y1error = A1error/A1

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

m,c = leastsquares(t, y2, y2error) 
ynew2 = m * t + c

m,c = leastsquares(t, y3, y3error) 
ynew3 = m * t + c

m,c = leastsquares(t, y4, y4error) 
ynew4 = m * t + c

# plt.plot(t,y1,color='red',ls='none')
# plt.errorbar(t,y1,yerr=y1error,fmt='ro')
# plt.plot(t,ynew1,color='black')
# plt.show()

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

allm1 = []
allA01 = []

allm2 = []
allA02 = []

allm3 = []
allA03 = []

allm4 = []
allA04 = []

randomres = []
randomerror = []

def bootstrap(y,y2,t,yerror,allm,allA0):
    #bootstrap method steps:
    #1) Calculate the best fit model and calculate the residuals (You will need to use these residuals for steps 2-5).  
    #2) Label each data point  
    #3) Draw at random a sample of N data with replacement from the residuals (use a uniform random number generator) and add the best fit model to these points (note: N is the number of data points you have)  
    #4) Recalculate the the best fit slope and intercept & store the results  
    #5) Repeat this process as many times as possible 

    #step 1 - residuals - not sure if the errors of residuals are right or not
    for j in range(3000):
        residual = y-y2
        plt.plot(t,residual,'o')
        residualerror = yerror
        plt.errorbar(t,residual,yerr=residualerror, ls = 'none')

        #step 2 - label points - already done from eg y[0] is first data point

        #step 3 - draw at random a sample of N data with replacement from the residuals and add the best fit model to these points
        #need to generate random numbers then change [number] to [new number]
        randomres = []
        randomerror = []

        for i in np.random.randint(0,len(residual), size=len(residual)):
            randomres.append(residual[i])
            randomerror.append(residualerror[i])

        data = randomres + y2
        dataerror = randomerror

        #step 4 = recalculate the best fit slope and intercept and store the results (step 5 - repeat as many times as possible)
        m,c = leastsquares(t,data,dataerror)

        A0 = np.exp(c)

        allm.append(m)
        allA0.append(A0)

    return allm,allA0

allm, allA0 = bootstrap(y1,ynew1,t,y1error,allm1,allA01)

# the first histogram is always messed up, not sure why
y,x,_=plt.hist(allm, bins = 50)

fig, axs = plt.subplots(2, 4)
axs[0, 0].hist(allm, bins = 50)
axs[1,0].hist(allA0, bins = 50)
axs[0, 0].set_title('Axis [0, 0]')

allm, allA0 = bootstrap(y2,ynew2,t,y2error,allm2,allA02)

axs[0, 1].hist(allm, bins = 50)
axs[1,1].hist(allA0, bins = 50)
axs[0, 1].set_title('Axis [0, 1]')

allm, allA0 = bootstrap(y3,ynew3,t,y3error,allm3,allA03)

axs[0, 2].hist(allm, bins = 50)
axs[1,2].hist(allA0, bins = 50)
axs[1, 0].set_title('Axis [1, 0]')

allm, allA0 = bootstrap(y4,ynew4,t,y4error,allm4,allA04)

axs[0, 3].hist(allm, bins = 50)
#last histogram also always messed up now for some reason
axs[1,3].hist(allA0, bins = 50)
axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()

#mean is value and standard deviation is error

#detector 1 calcs
meanm1 = np.mean(allm1,0)
standm1 = np.std(allm1,0)

meanA01 = np.mean(allA01,0)
standA01 = np.std(allA01,0)

print('detector 1')
print(meanm1,standm1)
print(meanA01,standA01)

halflife1 = 0.693/meanm1
halflifeerror1 = standm1/meanm1

print(halflife1,halflifeerror1)

#detector 2 calcs
meanm2 = np.mean(allm2,0)
standm2 = np.std(allm2,0)

meanA02 = np.mean(allA02,0)
standA02 = np.std(allA02,0)

print('detector 2')
print(meanm2,standm2)
print(meanA02,standA02)

halflife2 = 0.693/meanm2
halflifeerror2 = standm2/meanm2

print(halflife2,halflifeerror2)

#detector 3 calcs
meanm3 = np.mean(allm3,0)
standm3 = np.std(allm3,0)

meanA03 = np.mean(allA03,0)
standA03 = np.std(allA03,0)

print('detector 3')
print(meanm3,standm3)
print(meanA03,standA03)

halflife3 = 0.693/meanm3
halflifeerror3 = standm3/meanm3

print(halflife3,halflifeerror3)

#detector 4 calcs - these values are messed up ands i dont know why
meanm4 = np.mean(allm4,0)
standm4 = np.std(allm4,0)

meanA04 = np.mean(allA04,0)
standA04 = np.std(allA04,0)

print('detector 4')
print(meanm4,standm4)
print(meanA04,standA04)

halflife4 = 0.693/meanm4
halflifeerror4 = standm4/meanm4

print(halflife4,halflifeerror4)

#might have to somehow take into account large error bars and disregard this values
#have to write some code to do this?

#log error only works when delta x is much smaller than x

#Using the determination of 𝐴𝐴0 for each of the 4 detectors determine whether the 
#radiation is isotropic, by using the least-squares method (and bootstrap) to 
#determine the powerlaw index for 𝐴𝐴0 vs the distance of the detector from the 
#source (i.e. 𝐴𝐴0(𝑑𝑑)   =  𝐶𝐶 𝑑𝑑𝛼𝛼 with C a constant and 𝛼𝛼 the powerlaw index).

#A0 = C * d**alpha (C is a constant and alpha is the powerlaw index)
#logA0 = logd**alpha + logC
#logA0 = alpha*logd + logC

#need to make a graph with 4 points, one for each of the detectors

#detector 1 placed at 0.0500 +/- 0.0005 m
#detector 2 placed at 0.1000+/-0.0005 m
#detector 3 placed at 0.1800+/-0.0005 m
#detector 4 placed at 0.3000+/-0.0005 m

# all distance in metres
d = np.array([0.0500, 0.1000, 0.1800,0.3000])
derror = np.array([0.0005,0.0005,0.0005,0.0005])

A0 = np.array([meanA01,meanA02,meanA03,meanA04])
A0error = np.array([standA01,standA02,standA03,standA04])

#plot ln(A0) on y and ln(d) on x, ln(c) is intercept
y = np.log(A0)

yerror = A0error/A0

x = np.log(d)
xerror = derror/d 

print(y)
print(x)

#what about the uncertainty in d
m,c = leastsquares(x, y, yerror) 

#why are my uncertainties so small?
y2 = m * x + c
plt.plot(x,y2,color='red')
plt.plot(x,y,color='red',ls='none')
plt.errorbar(x,y,yerr=yerror,fmt='ro')
plt.show()

index = m
constant = np.exp(c)

allindex = []
allconstant = []

allm, allA0 = bootstrap(y,y2,x,yerror,allindex,allconstant)

# the first histogram is always messed up, not sure why
y,x,_=plt.hist(allm, bins = 50)
plt.show()
y,x,_=plt.hist(allindex, bins = 50)
plt.show()
y,x,_=plt.hist(allconstant, bins = 50)
plt.show()

meanindex = np.mean(allindex,0)
standindex = np.std(allindex,0)

meanconstant = np.mean(allconstant,0)
standconstant = np.std(allconstant,0)

print(meanindex,standindex)
print(meanconstant,standconstant)

#things to fix:
#why do my histograms plot extra points?