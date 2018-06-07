# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 09:20:54 2018

@author: sp4009
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib import rc
from matplotlib import rcParams
rc('font',**{'family':'serif','serif':['STIXGeneral'], 'weight': 'normal', 'size':'12'})
rcParams['mathtext.fontset'] = 'stix'
rcParams['mathtext.default'] = 'regular'
rcParams['lines.linewidth'] = 1

import numpy, scipy.interpolate, scipy.integrate, urllib.request, io, tarfile
import numericalunits as nu

#%%
'''

My way....

'''
h = 4.135667662e-15 #eV.s Planks constant
c = 299792458 #m/s speed of light
k = 8.6173303e-5 #eV/K Boltzmann constant
T = 300 #K
e = 1.60217662e-19
JeV = 1.6e-19
eVnm = 1239.84193

AM = pd.read_excel(r'C:\Users\sp4009\OneDrive - Imperial College London\PhDThesis\C_Intro\AM\ASTMG173\ASTMG173.xlsx', header = 1, index_col = 0)
'''
AM.index                wavelength                  nm
AM['eV']                electron volt               eV
AM['Global tilt']       Solar irrandiance           J / s / m**2 / nm
AM['GTphoton']          Photon intensity            photon / s / m**2 
AM['GTphotonTrapz']     Intergrated photon int.     photons[E-E_infinitiy] / s / m**2
AM['RR0']               Spectral radiance density   done per area
#AM['RR0Trapz']
AM['RRint']             Radiative recombination     ???
AM['Jsc']               current density             mA / cm**2
AM['Voc']               OC voltage                  V
AM['Pmax']              Max power                   ???
'''

# convert nm to eV
AM['eV'] = [eVnm/i for i in AM.index]
# convert power to photons with E={hc}{\lambda }
AM['GTphoton'] = AM['Global tilt'] / (AM['eV']*JeV) 
# intergrate photon count
AM['GTphotonTrapz'] = [np.trapz(AM['GTphoton'][:i], x = AM['GTphoton'][:i].index) for i in AM.index]

# my intergration didn't work so had to use Steve's little code below
#AM['RR'] = (2 * np.pi) / (c**2 * h**3) * ((AM['eV'])**2 / (np.exp(AM['eV'] / (k*T)) - 1) )
#AM['RRTrapz'] = [np.trapz(AM['RR'][:i], x = AM['RR'][:i].index) for i in AM.index]#intergrating from bandgap to max energy so needs to be reverse order of wavelength

# calculate the amount of black body radiative recombination, with a given voltage spillting fermi level into quasi fermi level
def RR02(Egap):
    integrand = lambda E : E**2 / (np.exp(E / (k * T)) - 1)
    integral = scipy.integrate.quad(integrand, Egap, AM['eV'].max(), full_output=1)[0]
    return ((2 * np.pi) / (c**2 * h**3)) * integral
AM['RRint'] = [RR02(i) for i in AM['eV']]

# calculate net current at voltage, V, from all photons absorbed - radiative recombination
def J(V):
    '''
    photons - recombination
    where recbomination = radiative recombination at 0V * exp(eV/kT) due to extra charges away from fermi level
    and convert from A / m2 to mA / cm2  
    '''
    #Wgap = int(1239.8/Egap)
    return e*(AM['GTphotonTrapz'] - AM['RRint'] * np.exp(e * V / (k * T))) *1000 / (100*100) 

# current at 0 V
AM['Jsc'] = J(0)
# voltage when photons absorbed = RR0, radiative recombination
AM['Voc'] = e * (k * T / e) * np.log(AM['GTphotonTrapz'] / AM['RRint'])  

# calculating the JV cureves for all wavelengths from Vrange
Vrange = np.linspace(-1,4, num=1000)
P = pd.DataFrame()
JV = pd.DataFrame()
for i in Vrange:
    df = J(i/e)
    df.columns = i
    JV  = pd.concat([JV,df], axis=1)
    P  = pd.concat([P,df.multiply(i)], axis=1)
P.columns =  Vrange
P = P.transpose()
JV.columns =  Vrange
JV = JV.transpose()

# add maximum power point to dataframe
AM['Pmax'] = P.max()

AM['FF'] = AM['Pmax']/(AM['Jsc']*AM['Voc'])

#%%
'''
Thesis plots...
'''
#%%
fig, ax = plt.subplots(figsize=(3,3))

ax1 = ax.twiny()
AM['Etr'].plot.area(ax=ax1, color='1', alpha=0)
#AM['Etr'].plot.area(ax=ax, color='0.5')
AM['Global tilt'][:2100].plot.area(ax=ax, color='0')


ax2=ax.twinx()
AM['Pmax'][:2100].plot(ax=ax2, color = 'r')
[axs.spines['right'].set_color('r') for axs in [ax2,ax1,ax]]
#[axs.tick_params(axis='y', colors='red') for axs in [ax2,ax1,ax]]
ax2.tick_params(axis='y', colors='red')
ax2.set_ylabel('$PCE_{max}$ (%)', color='r')

xlim = [200,2200]
ax.set_xlim(xlim)
ax1.set_xlim(xlim) #[1239.84/x for x in xlim])
ax.set_ylim(-0.12,2.4)
ax2.set_ylim(-1.8,36)
xticks = [400,800,1200,1600,2000]
ax.set_xticks(xticks)
ax1.set_xticks(xticks)
ax1.set_xticklabels([round(1240/tick,1) for tick in xticks], color='0')
ax1.set_xlabel('Photon energy (eV)', color='0')
ax.set_ylabel('AM1.5 Irradiance (W m$^{-2}$ nm$^{-1}$)')

# plot devices
devices= {'Perovskite: 22.7%' : [1240/1.6, 23],
          'Silicon: 27.6%' : [1240/1.1, 26.8],
          'GaAs: 29.3%' : [1240/1.39, 29.5],
          'OPV: 13.2%' : [1240/1.6, 13.2]}
          #'CdTe: 22.1%' : [1240/1.5, 22]}
for i in devices.keys():
    ax2.text(devices[i][0]+50, devices[i][1], i, color='r', fontsize=9, 
         bbox=dict(facecolor='1', edgecolor='none', pad=0))
    ax2.plot(*devices[i], 'r*')
ax2.text(1240/1.5+50, 21.5, 'CdTe: 22.1%', color='r', fontsize=9, va='top', ha='left',
         bbox=dict(facecolor='1', edgecolor='none', pad=0))
ax2.plot(1240/1.5, 21.5, 'r*')

fig.savefig(r'C:\Users\sp4009\OneDrive - Imperial College London\PhDThesis\C_Intro\AM/AM15.pdf', bbox_inches='tight',dpi=600)

#%% plotting photons
f, ax = plt.subplots(1,2, figsize=(6,3))
AM['Global tilt'].plot(ax=ax[0], xlim=(280,2000))
AM['GTphoton'].plot(ax=ax[1], xlim=(280,2000))
AM['GTphoton'][int(1239.84193/2)]
print('total power = ' + str(AM['Global tilt'].sum()))
print('total photons = ' + str(AM['GTphoton'].sum()))

#%% photons above bandgap
f, ax = plt.subplots(1,2, figsize=(6,3))
AM['GTphotonTrapz'].plot(ax=ax[0])
ax[1].plot(AM['eV'].values, AM['GTphotonTrapz'].values)
AM['GTphotonTrapz'][int(1240/1.1)]
print('total photons above 1.1 eV = ' + str(AM['GTphotonTrapz'][int(1240/1.1)]))

#%% plotting radiation
f, ax = plt.subplots(1,2, figsize=(6,3))
AM['RRint'].plot(ax=ax[0])
ax[1].plot(AM['eV'], AM['RRint'].values)
[ax.set_yscale('log') for ax in ax]
[ax.set_yscale('log') for ax in ax]

#%% plotting Jsc and Voc

f, ax = plt.subplots(1,2, figsize=(6,3))
AM['Voc'].plot(ax=ax[0])
ax0 = ax[0].twinx()
AM['Jsc'].plot(ax=ax0)
ax0.set_xlim(200,4500)

ax[1].plot(AM['eV'].values, AM['Voc'].values)
ax[1].plot(AM['eV'].values, AM['eV'].values, 'C0--')
[h(-0.5,4.5) for h in [ax[1].set_xlim, ax[1].set_ylim]]
plt.tight_layout()

print(JSC(1.1 * nu.eV) / (nu.mA / nu.cm**2)) # Steve's calcs
print(VOC(1.1 * nu.eV) / nu.V)
print(AM['Jsc'][int(eVnm/1.1)]) #1127
print(AM['Voc'][int(eVnm/1.1)])
#%% plotting Jsc and Voc and FF

clrs = ['C0','C1','C2']

f, ax0 = plt.subplots(figsize=(3,3))
AM['Voc'].plot(ax=ax0, c = clrs[0])
ax1 = ax0.twinx()
AM['Jsc'].plot(ax=ax1, c = clrs[1])
ax2 = ax0.twinx()
AM['FF'].plot(ax=ax2, c = clrs[2])

ax0.set_ylabel('V (V)', color=clrs[0])
ax1.set_ylabel('J (mA $cm^{-2}$)', color=clrs[1])
ax2.set_ylabel('FF', color=clrs[2])
ax0.set_ylim(-0.4,4)
ax1.set_ylim(-7,75)
ax2.set_ylim(-0.1,1.1)
ax0.spines['left'].set_color(clrs[0])
ax0.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['right'].set_color(clrs[1])
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_color(clrs[2])
ax0.tick_params(axis='y', colors=clrs[0])
ax1.tick_params(axis='y', colors=clrs[1])
ax2.tick_params(axis='y', colors=clrs[2])

ax2.spines['right'].set_position(('outward', 50))  


#%% plotting example JV curves
f, ax = plt.subplots(1,2, figsize=(6,3))
a = [300,350,400,450,500,550,600,650,700,800,900,950,1000,1100,1200,1300,1500,1700,1900,2300,2500,3000,3500]
JV[a[::-1]].plot(ax = ax[0], ylim=(-5,80), colormap = 'coolwarm', legend=None)
P[a[::-1]].plot(ax = ax[1], ylim=(-5,40), colormap = 'coolwarm', legend=None)

#%% Comparing to lines device

#lin JV
LinDevice = pd.read_csv(r'C:\Users\sp4009\OneDrive - Imperial College London\PhDThesis\C_Intro\AM\LinBestDevice\A3-MAPI-NMAI-age--251-rev_100_5.txt', header = 0, index_col = 0, sep='	')
LinDevice['Current Density /mAcm-2'] = LinDevice['Current Density /mAcm-2']
LinDevice['P'] = LinDevice['Current Density /mAcm-2'].multiply(LinDevice.index)

f, ax = plt.subplots(1,2, figsize=(6,3))
a = [775]
JV[a].plot(ax = ax[0], legend=False)
LinDevice['Current Density /mAcm-2'].multiply(-1).plot(ax=ax[0])
P[a].plot(ax = ax[1], legend=False)
LinDevice['P'].multiply(-1).plot(ax=ax[1])

[ax.set_xlim(-0.2,1.6) for ax in ax]
[ax.set_ylim(-5,32) for ax in ax]
ax[0].set_ylim(-2.5,28)
ax[1].set_ylim(-4,56)
ax[0].set_ylabel('J (mA $cm^{-2}$)')
ax[1].set_ylabel('P (mW $cm^{-2}$)')
ax[1].legend(['SQ limit w/ 1.6eV', 'Lin champion dev'])

plt.tight_layout()
f.savefig(r'C:\Users\sp4009\OneDrive - Imperial College London\PhDThesis\C_Intro\AM/Lin', bbox_inches='tight',dpi=600)


#%% calculate SQ efficency
'''
Steves normal
'''
import numpy, scipy.interpolate, scipy.integrate, urllib.request, io, tarfile
import numericalunits as nu

#%%
data_url = 'http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/compressed/ASTMG173.csv.tar'
download_as_bytes = urllib.request.urlopen(data_url).read()
download_as_file = io.BytesIO(download_as_bytes)
download_as_tarfile_object = tarfile.open(fileobj=download_as_file)
my_csv_file = download_as_tarfile_object.extractfile('ASTMG173.csv')
downloaded_array = numpy.genfromtxt(my_csv_file, delimiter=",", skip_header=2)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[:,[0,2]]

AM15[:,0] *= nu.nm
AM15[:,1] *= nu.W * nu.m**-2 * nu.nm**-1

wavelength_min = 280 * nu.nm
wavelength_max = 4000 * nu.nm
E_min = nu.hPlanck * nu.c0 / wavelength_max
E_max = nu.hPlanck * nu.c0 / wavelength_min   
AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])  


Tcell = 300 * nu.K

def SPhotonsPerTEA(Ephoton):
    wavelength = nu.hPlanck * nu.c0 / Ephoton
    return AM15interp(wavelength) * (1 / Ephoton) * (nu.hPlanck * nu.c0 / Ephoton**2)
print(SPhotonsPerTEA(2 * nu.eV) * (1 * nu.meV) * (1 * nu.m**2) * (1 * nu.s))

PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
print(solar_constant / (nu.W/nu.m**2))

def solar_photons_above_gap(Egap):
    return scipy.integrate.quad(SPhotonsPerTEA, Egap, E_max, full_output=1)[0]
print(solar_photons_above_gap(1.1 * nu.eV) * (1 * nu.m**2) * (1 * nu.s))


#%%

def RR0(Egap):
    integrand = lambda E : E**2 / (np.exp(E / (k * ( nu.eV / nu.K) * T * nu.K)) - 1)
    integral = scipy.integrate.quad(integrand, Egap, E_max, full_output=1)[0]
    return ((2 * np.pi) / ((c * nu.m / nu.s)**2 * (h* (nu.eV * nu.s))**3)) * integral


f, ax = plt.subplots(1,2, figsize=(6,3))

#AM['RRTrapz'].multiply(1e2).plot(ax=ax[0])

ax0=ax[0].twinx()
Egap_list = AM['eV'] * nu.eV
Egap_list2 =  numpy.linspace(0.31 * nu.eV, 4.42 * nu.eV, num=1000)
RR0_list = numpy.array([RR0(E) for E in Egap_list])
ax[0].plot(1240/(Egap_list / nu.eV) , RR0_list / nu.s, '-b')
RR0_list2 = numpy.array([RR0(E) for E in Egap_list2])
ax[0].plot(1240/(Egap_list2 / nu.eV) , RR0_list2 / nu.V, 'b.')


def RR02(Egap):
    integrand = lambda E : E**2 / (np.exp(E / (k * T)) - 1)
    integral = scipy.integrate.quad(integrand, Egap, AM['eV'].max(), full_output=1)[0]
    return ((2 * np.pi) / (c**2 * h**3)) * integral
AM['RRint'] = [RR02(i) for i in AM['eV']]



AM['RRint'] = [RR02(i) for i in AM['eV']]

ax[1].plot(1240/(Egap_list / nu.eV), AM['RRint'], color='r')
AM['RRint'].plot(ax=ax[1], color='b')
#AM['bb'].multiply(1e2).plot(ax=ax[0], color='r')
#[ax.set_ylimyscale('log') for ax in [ax[0],ax0]]
#[ax.set_ylim(2e-48,3e28) for ax in [ax[0],ax0]]
plt.tight_layout()
#%%

def current_density(V, Egap):
    return nu.e * (solar_photons_above_gap(Egap) - RR0(Egap) * np.exp(nu.e * V / (nu.kB * Tcell)))
def JSC(Egap):
    return current_density(0, Egap)
def VOC(Egap):
    return (nu.kB * Tcell / nu.e) * np.log(solar_photons_above_gap(Egap) / RR0(Egap))
print(JSC(1.1 * nu.eV) / (nu.mA / nu.cm**2))
print(VOC(1.1 * nu.eV) / nu.V)

def J(V):
    '''
    photons - recombination
    where recbomination = radiative recombination at 0V * exp(eV/kT) due to extra charges away from fermi level
    and convert from A / m2 to mA / cm2  
    '''
    #Wgap = int(1239.8/Egap)
    return e*(AM['GTphotonTrapz'] - AM['RRint'] * np.exp(e * V / (k * T))) *1000 / (100*100) 
'''
def current_density(V, Egap):
    return nu.e * (solar_photons_above_gap(Egap) - RR0(Egap) * np.exp(nu.e * V / (nu.kB * Tcell)))
print(JSC(1.1 * nu.eV) / (nu.mA / nu.cm**2))
'''
AM['Jsc'] = J(0)
AM['Voc'] = e * (k * T / e) * np.log(AM['GTphotonTrapz'] / AM['RRint'])  

print(AM['Jsc'][int(eVnm/1.1)])
print(AM['Voc'][int(eVnm/1.1)])




