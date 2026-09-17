#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:46:13 2023

@author: michele_mac
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import tkinter as tk
import random
from scipy.optimize import curve_fit


def Tc_finder(magna_medi_T, T_step, T_i):
    for i in range(len(magna_medi_T)-2):
        if magna_medi_T[i+2]+magna_medi_T[i]-2*magna_medi_T[i+1]>0.05 and T_i+i*T_step>1:
           return T_i+(1+i)*T_step
    return "buuu"
    

def Onsa_fit(T):
    if T<2.269:
       return (1-(math.sinh(2/T))**(-4))**(1/8)
    return 0
    

def E_per_site(pg,J):
    energy_mat=-J*(np.roll(pg,1,0)+np.roll(pg,-1,0)+np.roll(pg,1,1)+np.roll(pg,-1,1))
    E_flatt=np.ravel(energy_mat)*np.ravel(pg)
    return np.average(E_flatt)

def total_E(pg,J):
    energy_mat=-J*(np.roll(pg,1,0)+np.roll(pg,-1,0)+np.roll(pg,1,1)+np.roll(pg,-1,1))
    E_flatt=np.ravel(energy_mat)*np.ravel(pg)
    return np.sum(E_flatt)

def mean_magnetization_per_site(t_eq, L, record):
    
    data=record[::3*(L**2)]
    return abs(np.average(record))



def pittura (pg, L,C):
    
    wid=10
    x0=100
    y0=10
    a=x0, y0
    b=x0+wid, y0
    c=x0+wid, y0+wid
    d=x0, y0+wid
    sh=a,b,c,d
    col=["orange","","purple"]
    for i in range(L):
        for j in range(L):
            site_c=col[pg[i,j]+1]
            
            x0+=wid
            a=x0, y0
            b=x0+wid, y0
            c=x0+wid, y0+wid
            d=x0, y0+wid
            sh=a,b,c,d
            C.create_polygon(sh, fill=site_c)
        y0+=wid
        x0=100
        
        
#create the matrix of spin configurations
def create_play_ground(L):
    pg=np.random.randint(0,2,size=(L,L))*2-1
    return pg


#calculate energy of flip, site is the list with pos in the matrix
def site_flip_energy(site,pg,J,L):
    y=site[0]
    x=site[1]
    flipped=-pg[y,x]
    
    #each site interacts with its neighbors, first column interacts with last column, 
    #same for rows, so the matrix it's like the surface of a toroid
    
    E_f=-J*(pg[y,x-1]+pg[y,(x+1)%L]+pg[y-1,x]+pg[(y+1)%L,x])*flipped
    E_0=-J*(pg[y,x-1]+pg[y,(x+1)%L]+pg[y-1,x]+pg[(y+1)%L,x])*pg[y,x]
    return E_f-E_0



#function that flips the spin of a site, decision made with metropolis algorithm
def flip(site, pg,J,L, T):
    d_E=site_flip_energy(site, pg, J, L)
    s_E=np.sign(d_E)
    bu=np.array([s_E,-d_E/T])
    k=np.min(bu)
    if math.exp(k)>s_E*random.random():
        pg[site[0],site[1]]=-pg[site[0],site[1]]
    #if d_E<=0:
     #   pg[site[0],site[1]]=-pg[site[0],site[1]]
    #else:
     #   if math.e**(-d_E/T)>random.random():
      #       pg[site[0],site[1]]=-pg[site[0],site[1]]
    return pg



def road_to_nowhere(pg,t_eq,L, J, T,C):
    k=10000
    record=[]
    E_store=[]
    E_tot=[]
    for i in range(int(0.95*t_eq)):
        
        site=[np.random.randint(0,L), np.random.randint(0,L)]
        pg=flip(site, pg, J, L, T)
    
    
    """
    for i in range(int(0.05*t_eq)):
        record.append(np.average(pg))
        E_store.append(E_per_site(pg, J))
        E_tot.append(total_E(pg, J))
        site=[np.random.randint(0,L), np.random.randint(0,L)]
        pg=flip(site, pg, J, L, T)
        
    
    return [mean_magnetization_per_site(t_eq, L, record),np.average(E_store), np.average(E_tot)]
"""

    pittura(pg, L, C)    
        #if i%k==0 :
            #print(magnetization_per_site(pg))
    
    
    #plt.plot(np.arange(len(record)),record)



L=40
T=0.5
J=1
eq_distance=100
t_eq=(L**2)*eq_distance
pg=create_play_ground(L)


D=tk.Canvas(tk.Tk(), bg="white", height=1500, width=1500)

magna_medi_T=[]
mean_E_T=[]
E_tot=[]
magn_Onsa=[]
T_f=2
T_i=0.1
T_step=0.1
#for i in range(int((T_f-T_i)/T_step)):   
#    T=T_i+i*T_step
#    data_T=road_to_nowhere(pg, t_eq, L, J, T, D)
#    magna_medi_T.append(data_T[0])
#    mean_E_T.append(data_T[1])
#    E_tot.append(data_T[2])
#    magn_Onsa.append(Onsa_fit(T))

#m_E_T=np.array(E_tot)
#cap_T=(np.roll(m_E_T,-1)[:m_E_T.size-1]-m_E_T[:m_E_T.size-1])/T_step
#log_C=[]
#for i in range(int((3-1.5)/T_step)-1):
#    if cap_T[int(1.5/T_step)+i]<0 :
#        cap_T[int(1.5/T_step)+i]=1
#    log_C.append(math.log(cap_T[int(1.5/T_step)+i],L))


"""
plt.scatter(np.arange(start=T_i,stop=T_f, step=T_step),magna_medi_T)
plt.plot(np.arange(start=T_i,stop=T_f, step=T_step),magn_Onsa, color="red")
plt.xlabel("Temperature")
plt.ylabel("Mean Magnetization per site")

plt.show()
plt.plot(np.arange(start=T_i,stop=T_f, step=T_step),mean_E_T)
plt.xlabel("Temperature")
plt.ylabel("Mean Energy per site")
plt.show()
plt.plot(np.linspace(1.5, 3-T_step,int((3-1.5)/T_step) -1),log_C)
plt.xlabel("Temperature")
plt.ylabel("Heat Capacity (log)")
plt.show()


print(np.linspace(T_i, T_f-2*T_step,int((T_f-T_i)/T_step) -1))
"""
T=1

road_to_nowhere(pg, t_eq, L, J, T, D)

D.pack()
tk.Tk().mainloop()

#popt, pcov= curve_fit(Onsa_fit,np.arange(start=T_i,stop=T_f, step=T_step), magna_medi_T, p0=None)

#print('m',popt, pcov)

#plt.axline([0,0],xy2=None, slope=popt[0])

    
    
print(Tc_finder(magna_medi_T, T_step, T_i))



