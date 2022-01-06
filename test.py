import numpy as np 
import matplotlib.pyplot as plt 
import astropy.io.ascii as ascii
from GLS import GLS as GLS

t=np.arange(0,30.,0.01) 
Per=4.5 
mdl=np.sin(2.*np.pi*t/Per+np.pi/4) 
mdl_noise1=mdl+np.random.normal(0.,0.1,t.shape[0]) 
mdl_noise2=mdl+np.random.normal(0.,5.,t.shape[0]) 
sig1=t*0.+0.1 
sig2=t*0.+5. 

#Generate a random sampling of 1/5th the points and sort the array 
time_sampling=np.random.randint(0,t.shape[0],int((t.shape[0]/5))) 
time_sampling.sort() 
 
#Discard any duplicate points 
idx=np.unique(time_sampling) 
 
#Cut out two larger chunks, in this case from 70 to 110 and 220 to 240 
idx=np.hstack( (idx[0:70],idx[110:220],idx[240:])) 
 
#define the new time and data & sigma arrays. 
t_sample=t[idx].copy() 
mdl_noise2_sample=mdl_noise2[idx].copy() 
sig2_sample=sig2[idx].copy() 
 
plt.plot(t_sample,mdl_noise2_sample,'k.') 
plt.plot(t,mdl,linewidth=3) 
plt.show()

P=np.logspace(-1.3,np.log10(30./2.),int(1e5)) 
p_om,a_om,b_om=GLS(t,mdl_noise2,sig2,P) 
plt.plot(t,p_om)
plt.show()
 