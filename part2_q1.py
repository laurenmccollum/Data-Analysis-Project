import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii
from GLS import GLS as GLS

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

#you also want to oversample the periods to get a nice smooth periodogram
#mine is not smooth so whats up
P = np.logspace(-1.3,np.log10(30./2.),int(1e4)) 

p_om,a_om,b_om=GLS(t,RV,RVerror,P) 

plt.plot(p_om)
plt.show()

#question 2










#things needed to plot periodogram
#t,mdl_noise,sig,P

#mdl noise = mdl+np.random.normal(0.,5.,t.shape[0])
#mdl = np.sin(2.*np.pi*t/Per+np.pi/4) 
#need period, K and phase offset
#K = amplitude

#sig = t*0.+5.
#dont really understand multiplying time by 0

#period. no idea how to find

# t0 = t[1] #or could i just use t0 = 0???
# P = #period of orbit
# K = #radial velocity semi amplitude
# phase = #phase offset, a number between 0 and 1
# v = K * np.sin(2*np.pi*(t/P-phase))

#questions:
#how to find the period
#how to find radial velocity semi amplitude
#do I just choose a phase offset

# t=np.arange(0,30.,0.01) 
# Per=4.5 
# mdl=np.sin(2.*np.pi*t/Per+np.pi/4) 
# mdl_noise1=mdl+np.random.normal(0.,0.1,t.shape[0]) 
# mdl_noise2=mdl+np.random.normal(0.,5.,t.shape[0]) 
# sig1=t*0.+0.1 
# sig2=t*0.+5. 
# plt.plot(t,mdl_noise2,'k.') #this is the noise i think
# plt.plot(t,mdl,linewidth=3) 
# plt.show()

# plt.plot(t,mdl_noise1,'k.') #this is the noise i think
# plt.plot(t,mdl,linewidth=3) 
# plt.show()

# P=np.logspace(-1.3,np.log10(30./2.),int(1e5)) 
# p_om,a_om,b_om=GLS(t,mdl_noise2,sig2,P) 

# plt.plot(p_om) #is this just how you're meant to plot it?
# plt.show()