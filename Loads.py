import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from MOI import Ixx, Izz, Ixx2, Izz2, centx, centz
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

    #Preliminary Functions
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
inpCM2 = -0.11
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

#3rd wingbox
b = np.linspace(0.349, 0.1396, 600)
a = np.linspace(1.55, 0.620, 600)
y = 0.0012
part1 = np.ones(250)*0.02
part2 = np.ones(250)*0.02
part3 = np.ones(100)*0.016

A_design2 = np.concatenate((part1,part2,part3))

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

'''
plt.subplot(221)
plt.title("Stress Load Case 1")
plt.plot(span, s1, label ="Wing box 1")
plt.plot(span, s2, label ="Wing box 2")
plt.plot(span, s3, label ="Wing box 3")
plt.xlabel("Span location (m)")
plt.ylabel("Stress (Pa)")
plt.grid(True)
plt.legend()


plt.subplot(222)
plt.title("Stress Load Case 1")
plt.plot(span, as1, label ="Wing box 1")
plt.plot(span, as2, label ="Wing box 2")
plt.plot(span, as3, label ="Wing box 3")
plt.xlabel("Span location (m)")
plt.ylabel("Stress (Pa)")
plt.grid(True)
plt.legend()


plt.subplot(223)
plt.title("Stress Load Case 2")
plt.plot(span, bs1, label ="Wing box 1")
plt.plot(span, bs2, label ="Wing box 2")
plt.plot(span, bs3, label ="Wing box 3")
plt.xlabel("Span location (m)")
plt.ylabel("Stress (Pa)")
plt.grid(True)
plt.legend()


plt.subplot(224)
plt.title("Stress Load Case 2")
plt.plot(span, cs1, label ="Wing box 1")
plt.plot(span, cs2, label ="Wing box 2")
plt.plot(span, cs3, label ="Wing box 3")
plt.xlabel("Span location (m)")
plt.ylabel("Stress (Pa)")
plt.grid(True)
plt.legend()
plt.subplots_adjust(bottom = 0.1, top = 0.95, hspace = 0.4)
plt.show()
'''


E = 71*10**9


#Deflection graph for the lift

theta = integrate.cumtrapz(N_m, span, initial = 0)
v = -(1/(E *Ixx))*integrate.cumtrapz(theta, span, initial = 0)
#theta2 = integrate.cumtrapz(N_m2, span, initial = 0)
#v2 = -(1/(E *Ixx))*integrate.cumtrapz(theta2, span, initial = 0)

#Deflection drag for graph

theta_D = integrate.cumtrapz(T_m,span,initial = 0)
v_d = (1/(E *Izz))*integrate.cumtrapz(theta_D, span, initial = 0)


G = 26*10**9
t = 0.001
Am = (0.101*c_interpolate(span)+0.0985*c_interpolate(span))*(0.45*c_interpolate(span)-0)/2
S = 0.101*c_interpolate(span)+0.0985*c_interpolate(span)+2*np.sqrt((0.008625*c_interpolate(span))**(2)+(0.45*c_interpolate(span))**(2))
dtheta = IntM/(4*Am**2*G)*(S/t)
totaltheta = integrate.trapz(dtheta, span)
#print(totaltheta*180/sp.pi)

"""
plt.plot(span, v)
plt.title('Deflection of wing')
plt.axis('equal')
plt.grid(True)
plt.show()

plt.plot(span, v2)
"""
