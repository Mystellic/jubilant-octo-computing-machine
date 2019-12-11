import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#from Inertia import Ixx, Iyy, centx, centy


span = np.linspace(0, 10.8695, 600)
rho = 0.287368
Vel = 250
input_cl = 0.787
input_cd = 0.07648
input_cm = -0.11
df = pd.read_csv('AoA0.csv')


Y = np.array(df['y-span'].to_list())
Cl = np.array(df['Cl'].to_list())
Cd = np.array(df['ICd'].to_list())
Cm = np.array(df['CmAirf@chord/4'].to_list())
c = np.array(df['Chord'].to_list())


CL_interpolate = sp.interpolate.interp1d(Y,Cl,kind = 'cubic', fill_value = 'extrapolate')
CD_interpolate = sp.interpolate.interp1d(Y,Cd,kind = 'cubic', fill_value = 'extrapolate')
CM_interpolate = sp.interpolate.interp1d(Y,Cm,kind = 'cubic', fill_value = 'extrapolate')
c_interpolate = sp.interpolate.interp1d(Y,c, kind = 'linear', fill_value = 'extrapolate')
#Angke of attack 10

df10 = pd.read_csv('AoA10.csv')

Y10 = np.array(df10['y-span'].to_list())
Cl10 = np.array(df10['Cl'].to_list())
Cd10 = np.array(df10['ICd'].to_list())
Cm10 = np.array(df10['CmAirf@chord/4'].to_list())



CL10_interpolate = sp.interpolate.interp1d(Y,Cl10,kind = 'cubic', fill_value = 'extrapolate')
CD10_interpolate = sp.interpolate.interp1d(Y,Cd10,kind = 'cubic', fill_value = 'extrapolate')
CM10_interpolate = sp.interpolate.interp1d(Y,Cm10,kind = 'cubic', fill_value = 'extrapolate')


BIG_CL10 = 1.253876
BIG_CL = 0.414486

BIG_CD10 = 0.054363
BIG_CD = 0.006049

BIG_CM10 = -1.20412
BIG_CM = -0.478928


def Coeff(y,cl,cd,cm):

    #Preliminary FUnctions
    cd_profile = 0.00060458 - 0.00003337 * (y)
    CD0 = cd_profile + CD_interpolate(y)
    CD10 = cd_profile + CD10_interpolate(y)

    #Function for CL
    CL = CL_interpolate(y) + ((cl - BIG_CL)/(BIG_CL10-BIG_CL))*(CL10_interpolate(y)-CL_interpolate(y))
    CD = CD0 + ((cd - BIG_CD)/(BIG_CD10-BIG_CD))*(CD10 - CD0)
    CM = CM_interpolate(y) + ((cm - BIG_CM) / (BIG_CM10 - BIG_CM)) * (CM10_interpolate(y) - CM_interpolate(y))

    return(CL, CD, CM)


inpCl2 = -0.45
inpCD2 = 0.018
inpCM2 = 0.5
Cl, CD, CM = Coeff(span, input_cl, input_cd, input_cm)
Cl2, CD2, CM2 = Coeff(span, inpCl2, inpCD2, inpCM2)


L = 0.5*rho*Vel**2*Cl*c_interpolate(span)
D = 0.5*rho*Vel**2*CD*c_interpolate(span)
M = 0.5*rho*Vel**2*CM*c_interpolate(span)*c_interpolate(span)


L2 = 0.5*rho*Vel**2*Cl2*c_interpolate(span)
D2= 0.5*rho*Vel**2*CD2*c_interpolate(span)
M2 = 0.5*rho*Vel**2*CM2*c_interpolate(span)*c_interpolate(span)
#FIgure Out angle of attack

aoa = np.arcsin(((input_cl-BIG_CL)/(BIG_CL10 - BIG_CL))*np.sin(0.174532925))
aoa2 = np.arcsin(((inpCl2-BIG_CL)/(BIG_CL10 - BIG_CL))*np.sin(0.174532925))


N = np.cos(aoa)*L + np.sin(aoa)*D
T = np.cos(aoa)*D + np.sin(aoa)*L

N2 = np.cos(aoa2)*L2 + np.sin(aoa2)*D2
T2 = np.cos(aoa2)*D2 + np.sin(aoa2)*L2

Lfull = (10.8695/600)*sum(L)
Dfull = (10.8695/600)*sum(D)
#print(integrate.trapz(L2,span))

Nfull = integrate.trapz(N, span)
Tfull = integrate.trapz(T, span)
Mfull = integrate.trapz(M, span)

Nfull2 = integrate.trapz(N2, span)
Tfull2 = integrate.trapz(T2, span)
Mfull2 = integrate.trapz(M2, span)

# Now that the forces parallel and perpendicular to the chord of the airplane are known
#it is now time to turn these force diagrams into NVM diagrams for the sar of the wing.

N_v = Nfull - integrate.cumtrapz(N, span, initial = 0)
N_m = -integrate.trapz(N_v, span) + integrate.cumtrapz(N_v, span, initial = 0)

T_v = Tfull - integrate.cumtrapz(T, span, initial = 0)
T_m = integrate.trapz(T_v, span) - integrate.cumtrapz(T_v, span, initial = 0)

IntM = Mfull - integrate.cumtrapz(M, span, initial = 0)

N_v2 = Nfull2 - integrate.cumtrapz(N2, span, initial = 0)
N_m2 = -integrate.trapz(N_v2, span) + integrate.cumtrapz(N_v2, span, initial = 0)

T_v2 = Tfull2 - integrate.cumtrapz(T2, span, initial = 0)
T_m2 = integrate.trapz(T_v2, span) - integrate.cumtrapz(T_v2, span, initial = 0)

IntM2 = Mfull2 - integrate.cumtrapz(M2, span, initial = 0)

#TwistMom = 0.5*rho*Vel**2*CM*c_interpolate(span)*c_interpolate(span)


#Y is span
#Z up positive
#X is along chord
'''
#3rd wingbox
dens = 2810
b = np.linspace(0.349, 0.1396, 600)
bshort = np.linspace(0.3399, 0.1359644, 600)
a = np.linspace(1.55, 0.620, 600)

y = 0.0010
part1 = np.ones(250)*0.019
part2 = np.ones(250)*0.020
part3 = np.ones(100)*0.022

A_design2 = np.concatenate((part1,part2,part3))

Topplate = y * (1.5525 + 0.4 * 1.5525) * 10.8695 / 2
Sideplate = y * (0.349 + 0.4 * 0.349) * 10.8695 / 2
#Stringers = (250/600)*10.8695*0.02+(250/600)*10.8695*0.02+(100/600)*10.8695*0.015
Stringers = np.sum((np.ones(600)*(10.8695/600))*A_design2)
W = dens * (Topplate * 2 + Sideplate * 2 + Stringers)


Ixx_Design2 = 2*(((1/12)*y*b**3)+a*y*(b/2)**2) + A_design2*(b/2)**2

def stress_at_any_point(x,z):
    x = np.linspace(x, 0.4*x,600)
    z = np.linspace(z, 0.4 * z, 600)
    x = x - centx
    z = z- centz
    stress = -((N_m *z)/Ixx) - ((T_m *x)/Izz)
    stress2 = -((N_m *z)/Ixx2) - ((T_m *x)/Izz2)
    stress3 = -((N_m *z)/Ixx_Design2) - ((T_m *x)/Izz)
    return stress, stress2, stress3

def stress_at_any_point2(x,z):
    x = np.linspace(x, 0.4*x,600)
    z = np.linspace(z, 0.4 * z, 600)
    x = x - centx
    z = z- centz
    stress = -((N_m2 *z)/Ixx) - ((T_m2 *x)/Izz)
    stress2 = -((N_m2 *z)/Ixx2) - ((T_m2 *x)/Izz2)
    stress3 = -((N_m2 *z)/Ixx_Design2) - ((T_m2 *x)/Izz)
    return stress, stress2, stress3

s1, s2, s3 = stress_at_any_point(1.55, 0)
as1, as2, as3 = stress_at_any_point(0, 0.348)

bs1, bs2, bs3 = stress_at_any_point2(0, 0)
cs1, cs2, cs3 = stress_at_any_point2(1.55, 0.348)

E1 = 71*10**9
E2 = 107*10**9

#Deflection graph for the lift

theta = integrate.cumtrapz(N_m, span, initial = 0)
theta2 = integrate.cumtrapz(N_m2, span, initial = 0)
v = -(1/(E1 *Ixx))*integrate.cumtrapz(theta, span, initial = 0)
va = -(1/(E2 *Ixx2))*integrate.cumtrapz(theta, span, initial = 0)
vb = -(1/(E1 *Ixx_Design2))*integrate.cumtrapz(theta, span, initial = 0)

vb[250:500] = vb[250:500] - (vb[250]-vb[249])
vb[500:600] = vb[500:600] - (vb[500] - vb[499])

v2 = -(1/(E1 *Ixx))*integrate.cumtrapz(theta2, span, initial = 0)
va2 = -(1/(E2 *Ixx2))*integrate.cumtrapz(theta2, span, initial = 0)
vb2= -(1/(E1 *Ixx_Design2))*integrate.cumtrapz(theta2, span, initial = 0)

vb2[250:500] = vb2[250:500] - (vb2[250]-vb2[249])
vb2[500:600] = vb2[500:600] - (vb2[500] - vb2[499])

#theta2 = integrate.cumtrapz(N_m2, span, initial = 0)
#v2 = -(1/(E *Ixx))*integrate.cumtrapz(theta2, span, initial = 0)

#Deflection drag for graph

#theta_D = integrate.cumtrapz(T_m,span,initial = 0)
#v_d = (1/(E *Izz))*integrate.cumtrapz(theta_D, span, initial = 0)


G1 = 26*10**9
G2 = 42.1*10**9
t1 = 0.001
t2 = 0.001
t3 = 0.001
Am = (0.101*c_interpolate(span)+0.0985*c_interpolate(span))*(0.45*c_interpolate(span)-0)/2
S = 0.101*c_interpolate(span)+0.0985*c_interpolate(span)+2*np.sqrt((0.008625*c_interpolate(span))**(2)+(0.45*c_interpolate(span))**(2))
dtheta1 = IntM/(4*Am**2*G1)*(S/t1)


dtheta2 = IntM/(4*Am**2*G2)*(S/t2)
dtheta3 = IntM/(4*Am**2*G1)*(S/t3)
totaltheta1 = integrate.cumtrapz(dtheta1, span, initial = 0)
totaltheta1 = totaltheta1*180/np.pi
totaltheta2 = integrate.cumtrapz(dtheta2, span, initial = 0)
totaltheta2 = totaltheta2*180/np.pi
totaltheta3 = integrate.cumtrapz(dtheta3, span, initial = 0)
totaltheta3 = totaltheta3*180/np.pi


kv = 1.5
ks = 5.5
vp = 0.33
E = 71*10**9
sf = 1.1
kc = 4

#Skin buckling
Areas = S*t+AS

pitch_stringer = t*((abs(s1)*Areas*12*(1-vp**2))/(np.pi**2*kc*E))**-0.5
#plt.plot(span[0:500], pitch_stringer[0:500])
#plt.show()


#Column buckling
#-316177152.627236

I = (abs(s[0])*0.7246**2)/(0.25*np.pi**2*E)
print(I)
print(Ixx[0])
#Shear stress


tau_ave_shear = N_v/(b*t+bshort*t)

tau_max_shear = kv*tau_ave_shear
tau = IntM/(2*Am*t)
#plt.plot(span, IntM)

tau_max_max = tau + tau_max_shear
#tcr = (np.pi**2*ks*E)/(12(1-vp**2))*(t/pitch_stiffener)**2
pitch_stiffener = t*((12*(1-vp**2)*tau_max_max*sf)/(np.pi**2*ks*E))**-0.5

pitch_stiffener0 = pitch_stiffener[0]   #0 delta 2.5
pitch_stiffener1 = pitch_stiffener[138] #2.5 delta 2.5
pitch_stiffener2 = pitch_stiffener[276] #5.0 delta 2.5
pitch_stiffener3 = pitch_stiffener[414] #7.5delta 0.5
pitch_stiffener4 = pitch_stiffener[442] #8.0 delta 0.5
pitch_stiffener5 = pitch_stiffener[469] #8.5 delta 0.5
pitch_stiffener6 = pitch_stiffener[497] #9.0 delta 0.5
pitch_stiffener7 = pitch_stiffener[524] #9.5 delta 1.3695
stiffeners = 2.5/pitch_stiffener0 + 2.5/pitch_stiffener1 + 2.5/pitch_stiffener2 + 0.5/pitch_stiffener3 + 0.5/pitch_stiffener4 + 0.5/pitch_stiffener5 + 0.5/pitch_stiffener6 + 1.3695/pitch_stiffener7
'''



