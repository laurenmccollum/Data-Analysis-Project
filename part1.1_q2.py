#detector cross section = 3.142 +/- 0.001 cm^2

#detector 1 placed at 0.0500 +/- 0.0005 m
#detector 2 placed at 0.1000+/-0.0005 m
#detector 3 placed at 0.1800+/-0.0005 m
#detector 4 placed at 0.3000+/-0.0005 m

import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii

# here's how to read in some data: 
data = ascii.read('part1_data1.dat') 

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
y = np.log(A1)

yerror = A1error/A1

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

m,c = leastsquares(t, y, yerror) 

y2 = m * t + c
plt.plot(t,y,color='red',ls='none')
plt.errorbar(t,y,yerr=yerror,fmt='ro')
plt.plot(t,y2,color='black')
plt.show()

allm = []
allA0 = []

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
        random = np.random.randint(0,len(residual), size=len(residual))

        for i in range(len(residual)):
            residual[i] = residual[random[i]]
            residualerror[i] = residualerror[random[i]]

        y = residual + y2
        yerror = residualerror

        #step 4 = recalculate the best fit slope and intercept and store the results (step 5 - repeat as many times as possible)
        m,c = leastsquares(t,y,yerror)
        y2 = m * t + c

        A0 = np.exp(c)

        allm.append(m)
        allA0.append(A0)

    return allm,allA0

allm, allA0 = bootstrap(y,y2,t,yerror,allm,allA0)

# the first histogram is always messed up, not sure why
y,x,_=plt.hist(allm, bins = 50)
plt.show()
y,x,_=plt.hist(allm, bins = 50)
plt.show()
y,x,_=plt.hist(allA0, bins = 50)
plt.show()

#mean is value and standard deviation is error

meanm = np.mean(allm,0)
standm = np.std(allm,0)

meanA0 = np.mean(allA0,0)
standA0 = np.std(allA0,0)

print(meanm,standm)
print(meanA0,standA0)

halflife = 0.693/meanm
halflifeerror = standm/meanm

print(halflife,halflifeerror)

#why do I get different values every time? should this happen or not?

#might have to somehow take into account large error bars and disregard this values
#have to write some code to do this?

#log error only works when delta x is much smaller than x