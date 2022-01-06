# -*- coding: utf-8 -*-

import numpy as np

#Routine to calculate the Generalised Lomb Scargle (see Zechmeister & Kuster 2009).
#Inputs:
# t - time 
# y - data values
# sig - uncertainties on y
# period - array containing the periods on which the periodogram should be calculated.
#
#Outputs:
# p_om - the periodogram
# a_om - the amplitude of the cosine term 
# b_om - the amplitude of the sine term
# Note: total amplitude is sqrt(a_om+b_om)
def GLS(t,y,sig,period=np.logspace(-1.,1.4,int(1e5))):
    W=np.sum(1./sig**2)
    w=1./W * 1./sig**2

    Y=np.sum(w*y)
    YYh=np.sum(w*y*y)
    p_om=np.zeros(period.shape[0])
    a_om=np.zeros(period.shape[0])
    b_om=np.zeros(period.shape[0])
    for i,P in enumerate(period):
        omt=2.*np.pi*t/P
        cos_omt=np.cos(omt)
        sin_omt=np.sin(omt)
        C=np.sum(w*cos_omt)
        S=np.sum(w*sin_omt)
        YCh=np.sum(w*y*cos_omt)
        YSh=np.sum(w*y*sin_omt)
        CCh=np.sum(w*cos_omt*cos_omt)
        SSh=np.sum(w*sin_omt*sin_omt)
        CSh=np.sum(w*cos_omt*sin_omt)
        YY=YYh-Y*Y
        YC=YCh-Y*C
        YS=YSh-Y*S
        CC=CCh-C*C
        SS=SSh-S*S
        CS=CSh-C*S
        D=CC*SS-CS**2
        p_om[i]=(SS*YC**2+CC*YS**2-2*CS*YC*YS)/(YY*D)
        a_om[i]=(YC*SS-YS*CS)/D
        b_om[i]=(YS*CC-YC*CS)/D
    return p_om,a_om,b_om