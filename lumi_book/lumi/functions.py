import numpy as np
import math 
import matplotlib.pyplot as plt
from math import factorial
from matplotlib.pyplot import cm
from scipy import special

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

def psi(Q,n,alpha):
    
    N = np.sqrt(1/(2**n*factorial(n))) * (alpha/np.pi)**(1/4)
    psi= N * special.eval_hermite(n,np.sqrt(alpha)*Q) * np.exp(-(alpha*Q**2)/2)
    
    return psi

def energy(n,omega):
    E=(n+0.5)*omega
    return E

def potential(omega,Q):
    return (1/2)*omega**2*Q**2
    
def inverse_potential(E,omega):
    return np.sqrt((2*E))/omega

def get_S(omega,Q_offset):
    return omega*Q_offset**2/2

def get_omega(mu,k):
    return np.sqrt(k/mu)

def lineshape(n,Energy,E_zpl,omega,Q_offset,broadening=0.1):
    sigma=broadening*omega
    A=np.zeros(len(Energy))
    S=get_S(omega,Q_offset)
    for i in range(n):
        A+= np.exp(-S)*(S**i)/(factorial(i)) * (1/sigma*np.sqrt(2*np.pi))*np.exp(-0.5*(((Energy-(E_zpl-omega*i))/sigma)**2))

    return A/max(A)

def lineshape_envelope(n,Energy,E_zpl,omega,Q_offset,broadening=0.7):
    
    sigma_envelop=omega*broadening
    A_envelop=np.zeros(len(Energy))
    S=get_S(omega,Q_offset)        
    for i in range(n):
        A_envelop+=np.exp(-S)*(S**i)/(factorial(i)) * (1/sigma_envelop*np.sqrt(2*np.pi))*np.exp(-0.5*(((Energy-(E_zpl-omega*i))/sigma_envelop)**2))
    
    return A_envelop/max(A_envelop)



