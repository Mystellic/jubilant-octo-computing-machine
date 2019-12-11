import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from Inertia import Inertia_Calc
from LOADS import span, c_interpolate, N_v, N_m, T_v, T_m, IntM, N_v2, N_m2, T_v2, T_m2, IntM2

t = 0.003
AS = np.ones(600)*0.0142

AS_lower = np.ones(600)*0.0142

b = np.linspace(0.349, 0.1396, 600)
bshort = np.linspace(0.3399, 0.1359644, 600)
a = np.linspace(1.55, 0.620, 600)

def stress_at_any_point(x, y, I_xx, I_yy, centx, centy, Nm, Tm):
    x = np.linspace(x, 0.4*x,600)
    y = np.linspace(y, 0.4 * y, 600)
    x = x - centx
    y = y - centy
    stress = ((Nm *y)/I_xx) - ((Tm *x)/I_yy)
    return stress

E1 = 71*10**9
E2 = 107*10**9
G1 = 26*10**9
G2 = 42.1*10**9


#Deflections
'''
theta = integrate.cumtrapz(N_m, span, initial = 0)
theta2 = integrate.cumtrapz(N_m2, span, initial = 0)
v = -(1/(E1 *Ixx))*integrate.cumtrapz(theta, span, initial = 0)
v2 = -(1/(E1 *Ixx))*integrate.cumtrapz(theta2, span, initial = 0)
'''

#Angle of twist

Am = (0.101*c_interpolate(span)+0.0985*c_interpolate(span))*(0.45*c_interpolate(span)-0)/2

S = 0.101*c_interpolate(span)+0.0985*c_interpolate(span)+2*np.sqrt((0.008625*c_interpolate(span))**(2)+(0.45*c_interpolate(span))**(2))
'''
dtheta1 = IntM/(4*Am**2*G1)*(S/t)
totaltheta1 = (integrate.cumtrapz(dtheta1, span, initial = 0)) *180/np.pi # Degrees
'''



vp = 0.33

#Shear stress
kv = 1.5
sf = 1.1
ks = 5.5

tau_ave_shear = N_v/(b*t+bshort*t)

tau_max_shear = kv*tau_ave_shear
tau = IntM/(2*Am*t)

tau_max_max = tau + tau_max_shear
pitch_stiffener = t*((12*(1-vp**2)*tau_max_max*sf)/(np.pi**2*ks*E1))**-0.5
#Varying pitch
pitch_stiffener0 = pitch_stiffener[0]   #0 delta 2.5
pitch_stiffener1 = pitch_stiffener[138] #2.5 delta 2.5
pitch_stiffener2 = pitch_stiffener[276] #5.0 delta 2.5
pitch_stiffener3 = pitch_stiffener[414] #7.5delta 0.5
pitch_stiffener4 = pitch_stiffener[442] #8.0 delta 0.5
pitch_stiffener5 = pitch_stiffener[469] #8.5 delta 0.5
pitch_stiffener6 = pitch_stiffener[497] #9.0 delta 0.5
pitch_stiffener7 = pitch_stiffener[524] #9.5 delta 1.3695
stiffeners = 2.5/pitch_stiffener0 + 2.5/pitch_stiffener1 + 2.5/pitch_stiffener2 + 0.5/pitch_stiffener3 + 0.5/pitch_stiffener4 + 0.5/pitch_stiffener5 + 0.5/pitch_stiffener6 + 1.3695/pitch_stiffener7

#Rib_locations = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 10.8695]

Rib_locations = np.arange(0, 10.8695, 0.2)

Rib_locations[-1] = 10.8695
Lengths = []
spanlocations = []
for i in range(1, len(Rib_locations)):
    Length = Rib_locations[i] - Rib_locations[i-1]
    Lengths.append(Length)
    spanlocations.append(int(Rib_locations[i-1]/0.018116))
Lengths = np.array(Lengths)
spanlocations.append(600)

#Rib section start, Rib section stop, n_s, A_s, P_s, H_s, b_s, Pitch, Actual pitch
sections = np.array([0,0,0,0,0,0])

Stringer_range = np.arange(1, 40 ,2)
Stringer_range[0] = 2

#Upper skin
for x in range(np.size(Rib_locations)-1):
#for x in range(0,2):
    Length = Lengths[x]
    spanloc_start = spanlocations[x]
    spanloc_end = spanlocations[x+1]
    


    #stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, N_m, T_m)
    for i in Stringer_range:
        A_s = 0.000568
        AS_section = A_s * i
        AS[spanloc_start:spanloc_end] = AS_section
        Areas = S*t+AS
        kc = 4
        t_s = 0.004
        j = 1/4
        H_s = 0.1136
        b_s = j*H_s
        centroid_s = (0.5*(H_s**2)*t_s)/A_s
        Ixx_stringers = (1/12*(H_s**3)*t_s+A_s*(centroid_s)**2)*i

        centx, centy, Ixx, Iyy = Inertia_Calc(t, AS, AS_lower)
        stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, centx, centy, N_m, T_m)
        I = (abs(stress[spanloc_start])*Length**2)/(0.25*np.pi**2*E1)
        pitch_stringer = t*((abs(stress[spanloc_start])*Areas[spanloc_start]*12*(1-vp**2))/(np.pi**2*kc*E1))**-0.5
        Apitch = a[spanloc_start]/i
        if Ixx_stringers > I:
            section = np.array([Rib_locations[x], Rib_locations [x+1], i, AS_section, pitch_stringer, Apitch])
            sections = np.vstack((sections, section))
            break
sections = np.delete(sections, 0,0)

df = pd.DataFrame(sections)
df.columns = ['Section start', 'Section end', 'Amount','AS_section', ' RPitch', 'APitch']
print("Upper skin")
print(df)

stress = stress_at_any_point(-0.775, -0.1745, Ixx, Iyy, centx, centy, N_m, T_m)

#Lower skin
sections_low = np.array([0,0,0,0,0,0])
for x in range(np.size(Rib_locations)-1):

    Length = Lengths[x]
    spanloc_start = spanlocations[x]
    spanloc_end = spanlocations[x+1]

    #stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, N_m, T_m)
    for i in Stringer_range:
        A_s = 0.000568
        AS_lowsection = A_s * i
        AS_lower[spanloc_start:spanloc_end] = AS_lowsection
        Areas = S*t+AS
        kc = 4
        t_s = 0.004
        j = 1/4
        H_s = 0.1136
        b_s = j*H_s
        centroid_s = (0.5*(H_s**2)*t_s)/A_s
        Ixx_stringers = (1/12*(H_s**3)*t_s+A_s*(centroid_s)**2)*i

        centx, centy, Ixx, Iyy = Inertia_Calc(t, AS, AS_lower)
        stress = stress_at_any_point(-0.775, -0.1745, Ixx, Iyy, centx, centy, N_m2, T_m2)
        I = (abs(stress[spanloc_start])*Length**2)/(0.25*np.pi**2*E1)
        pitch_stringer_low = t*((abs(stress[spanloc_start])*Areas[spanloc_start]*12*(1-vp**2))/(np.pi**2*kc*E1))**-0.5
        Apitch = a[spanloc_start]/i
        if Ixx_stringers > I:
            section_low = np.array([Rib_locations[x], Rib_locations [x+1], i, AS_lowsection, pitch_stringer_low, Apitch])
            sections_low = np.vstack((sections_low, section_low))
            break
sections_low = np.delete(sections_low, 0,0)

df_low = pd.DataFrame(sections_low)
df_low.columns = ['Section start', 'Section end', 'Amount','AS_section', ' RPitch', 'APitch']
print("Lower skin")
print(df_low)

#Moment of inertia top plate
'''
n_s = 25
A_s = AS/n_s 
t_s = 0.004
j = 1/4
H = A_s/((1+j)*t_s)
b = j*H

#L_s = (A_s/2)/t_s
#centroid_s = (0.5*(H**2)*t_s)/A_s
#Ixx_stringers = (1/12*(H**3)*t_s+A_s*(centroid_s)**2)*n_s



centx, centy, Ixx, Iyy = Inertia_Calc(0.003, 0.010792, 0)

stress = stress_at_any_point(0.775, 0.1745, Ixx, Iyy, N_m, T_m)



print("Area of stringers = ", A_s)
'''
#Skin buckling
kc = 4
Areas = S*t+AS







#Weight of yo momma
dens = 2810
Topplate = t*(1.5525+0.4*1.5525)*10.8695/2
Sideplate = t*(0.36+0.4*0.36)*10.8695/2
Stringers = integrate.trapz(AS, span)
Stringers_lower = integrate.trapz(AS_lower, span)

TotalArea = 0
for x in range(np.size(Rib_locations)-1):
    Area = Am[spanlocations[x]]*0.001
    TotalArea = TotalArea + Area

Weight = dens*(Topplate*2+Sideplate*2+Stringers+Stringers_lower+TotalArea)
print("Weight =", Weight)

